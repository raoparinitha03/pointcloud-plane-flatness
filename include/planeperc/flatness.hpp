#pragma once
// Surface-flatness analysis: how far do the points of a surface deviate from its fitted plane?
// RANSAC inliers only contain points within `threshold` of the plane, so a bump taller than the
// threshold would be hidden.  We therefore re-collect a wider "band" of points around the plane
// (restricted to the plane's own extent and excluding points owned by other planes).
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include "planeperc/types.hpp"

namespace planeperc {

struct FlatnessStats {
  size_t n = 0;
  double rms = 0;      // root-mean-square deviation (m)
  double maxAbs = 0;   // largest |cell-median deviation| (m)
  double p95Abs = 0;   // 95th percentile of |deviation| (m)
  double maxPos = 0;   // largest deviation toward the normal (m)  (bump toward viewer)
  double maxNeg = 0;   // largest deviation away from the normal (m) (dent)
};

// Build an orthonormal in-plane basis (u, v) from the normal n.
inline void planeBasis(const Vec3& n, Vec3& u, Vec3& v) {
  const Vec3 helper = (std::fabs(n.x) < 0.9) ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
  u = cross(n, helper);
  u = u * (1.0 / norm(u));
  v = cross(n, u);
}

// Indices of points that belong to the surface "band" of `pl`.
inline std::vector<size_t> surfaceBand(const Cloud& cloud, const Plane& pl, const std::vector<size_t>& planeInliers,
                                       const std::vector<int8_t>& ownedByOther, double band) {
  Vec3 u, v;
  planeBasis(pl.n, u, v);
  double umin = 1e300, umax = -1e300, vmin = 1e300, vmax = -1e300;
  for (size_t i : planeInliers) {
    const double a = dot(cloud[i], u), b = dot(cloud[i], v);
    umin = std::min(umin, a);
    umax = std::max(umax, a);
    vmin = std::min(vmin, b);
    vmax = std::max(vmax, b);
  }
  std::vector<size_t> out;
  for (size_t i = 0; i < cloud.size(); ++i) {
    if (ownedByOther[i]) continue;
    if (std::fabs(pl.dist(cloud[i])) > band) continue;
    const double a = dot(cloud[i], u), b = dot(cloud[i], v);
    if (a < umin || a > umax || b < vmin || b > vmax) continue;
    out.push_back(i);
  }
  return out;
}

// Flatness statistics.
//  * rms / p95 are computed over all band points (cheap, but influenced by leftover outliers).
//  * maxPos / maxNeg / maxAbs use a ROBUST flatness map: bin the points into cellSize x cellSize
//    cells in the plane's (u,v) coordinates, take the MEDIAN signed deviation of each cell, and
//    report the extremes of those medians.  A lone outlier cannot move a median, but a real
//    bump (many points shifted together) does.
inline FlatnessStats flatnessStats(const Cloud& cloud, const Plane& pl, const std::vector<size_t>& idx,
                                   double cellSize = 0.05, size_t minCellPts = 3) {
  FlatnessStats s;
  s.n = idx.size();
  if (idx.empty()) return s;
  Vec3 u, v;
  planeBasis(pl.n, u, v);

  std::vector<double> absd;
  absd.reserve(idx.size());
  double sq = 0;
  std::unordered_map<int64_t, std::vector<double>> cells;
  const double inv = 1.0 / cellSize;
  for (size_t i : idx) {
    const double d = pl.dist(cloud[i]);
    absd.push_back(std::fabs(d));
    sq += d * d;
    const int64_t cu = static_cast<int64_t>(std::floor(dot(cloud[i], u) * inv));
    const int64_t cv = static_cast<int64_t>(std::floor(dot(cloud[i], v) * inv));
    cells[(cu << 32) ^ (cv & 0xffffffffLL)].push_back(d);
  }
  s.rms = std::sqrt(sq / idx.size());
  std::sort(absd.begin(), absd.end());
  s.p95Abs = absd[static_cast<size_t>(0.95 * (absd.size() - 1))];

  for (auto& kv : cells) {
    auto& dv = kv.second;
    if (dv.size() < minCellPts) continue;
    std::nth_element(dv.begin(), dv.begin() + dv.size() / 2, dv.end());
    const double med = dv[dv.size() / 2];
    s.maxPos = std::max(s.maxPos, med);
    s.maxNeg = std::min(s.maxNeg, med);
  }
  s.maxAbs = std::max(s.maxPos, -s.maxNeg);
  return s;
}

// Map a signed deviation to a blue(-) - white(0) - red(+) colour. `range` = deviation shown fully saturated.
inline void deviationColor(double dev, double range, uint8_t& r, uint8_t& g, uint8_t& b) {
  const double t = std::max(-1.0, std::min(1.0, dev / range));
  if (t >= 0) {
    r = 255;
    g = b = static_cast<uint8_t>(255 * (1.0 - t));
  } else {
    b = 255;
    r = g = static_cast<uint8_t>(255 * (1.0 + t));
  }
}

}  // namespace planeperc
