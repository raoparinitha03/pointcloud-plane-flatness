#pragma once
// RANSAC plane detection from scratch.
//   repeat N times:
//     1. pick 3 random points  -> they define a plane (normal = cross product)
//     2. count points whose distance to that plane < threshold  (the "consensus set")
//     3. keep the plane with the largest consensus set
//   N is adapted:  N = log(1-p) / log(1 - w^3),  w = best inlier ratio so far, p = confidence.
//   Finally refine: least-squares fit on the inliers, re-collect inliers, repeat.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>
#include <random>
#include "planeperc/plane_fit.hpp"

namespace planeperc {

struct RansacParams {
  double threshold = 0.01;  // metres: max point-to-plane distance to count as inlier
  int maxIters = 1000;
  double confidence = 0.99;
  size_t minInliers = 200;
  int refineIters = 3;
  uint64_t seed = 42;
};

struct PlaneResult {
  Plane plane;
  std::vector<size_t> inliers;  // indices into the cloud passed to ransacPlane
  int iterations = 0;
};

inline std::vector<size_t> collectInliers(const Cloud& pts, const Plane& pl, double thr) {
  std::vector<size_t> out;
  for (size_t i = 0; i < pts.size(); ++i)
    if (std::fabs(pl.dist(pts[i])) < thr) out.push_back(i);
  return out;
}

inline std::optional<PlaneResult> ransacPlane(const Cloud& pts, const RansacParams& prm) {
  const size_t n = pts.size();
  if (n < 3 || n < prm.minInliers) return std::nullopt;

  std::mt19937_64 rng(prm.seed);
  std::uniform_int_distribution<size_t> pick(0, n - 1);

  Plane best;
  size_t bestCount = 0;
  int needed = prm.maxIters;
  int it = 0;
  for (; it < needed; ++it) {
    const size_t ia = pick(rng), ib = pick(rng), ic = pick(rng);
    if (ia == ib || ib == ic || ia == ic) continue;
    const Vec3 &a = pts[ia], &b = pts[ib], &c = pts[ic];
    Vec3 nrm = cross(b - a, c - a);
    const double len = norm(nrm);
    if (len < 1e-12) continue;  // collinear sample -> degenerate plane
    nrm = nrm * (1.0 / len);
    const Plane cand{nrm, -dot(nrm, a)};

    size_t count = 0;
    for (const Vec3& p : pts)
      if (std::fabs(cand.dist(p)) < prm.threshold) ++count;

    if (count > bestCount) {
      bestCount = count;
      best = cand;
      // adaptive iteration count
      const double w = static_cast<double>(count) / n;
      const double w3 = w * w * w;
      if (w3 >= 1.0 - 1e-12) {
        needed = it + 1;
      } else if (w3 > 1e-12) {
        const double N = std::log(1.0 - prm.confidence) / std::log(1.0 - w3);
        needed = std::min(prm.maxIters, static_cast<int>(std::ceil(N)));
      }
    }
  }
  if (bestCount < prm.minInliers) return std::nullopt;

  // Refinement: LS fit on inliers -> new inliers -> repeat
  std::vector<size_t> inl = collectInliers(pts, best, prm.threshold);
  for (int r = 0; r < prm.refineIters; ++r) {
    Plane refined;
    if (!fitPlaneLS(pts, inl, refined)) break;
    best = refined;
    inl = collectInliers(pts, best, prm.threshold);
  }
  if (inl.size() < prm.minInliers) return std::nullopt;
  return PlaneResult{best, std::move(inl), it};
}

// Sequential RANSAC: find the biggest plane, remove its inliers, repeat.
inline std::vector<PlaneResult> extractPlanes(const Cloud& cloud, const RansacParams& prm, int maxPlanes) {
  std::vector<PlaneResult> planes;
  std::vector<char> used(cloud.size(), 0);
  for (int k = 0; k < maxPlanes; ++k) {
    Cloud sub;
    std::vector<size_t> map;  // sub index -> original index
    for (size_t i = 0; i < cloud.size(); ++i)
      if (!used[i]) {
        sub.push_back(cloud[i]);
        map.push_back(i);
      }
    RansacParams p = prm;
    p.seed = prm.seed + 7919ull * (k + 1);
    auto res = ransacPlane(sub, p);
    if (!res) break;
    PlaneResult out = *res;
    for (size_t& i : out.inliers) {
      i = map[i];
      used[i] = 1;
    }
    planes.push_back(std::move(out));
  }
  return planes;
}

}  // namespace planeperc
