#pragma once
// End-to-end pipeline:  raw cloud -> voxel downsample -> outlier removal -> sequential RANSAC
//                       -> least-squares refinement -> flatness analysis per plane.
#include <chrono>
#include "planeperc/flatness.hpp"
#include "planeperc/outlier_removal.hpp"
#include "planeperc/ransac_plane.hpp"
#include "planeperc/voxel_grid.hpp"

namespace planeperc {

struct PipelineParams {
  double voxelLeaf = 0.03;   // m, 0 disables
  size_t sorK = 16;
  double sorStd = 1.5;
  RansacParams ransac;
  int maxPlanes = 6;
  double flatBand = 0.05;    // m, half-width of band used for flatness statistics
  double flatCell = 0.05;    // m, cell size of the robust flatness map (auto-raised to 2x voxel leaf)
  Vec3 viewpoint{0, 0, 0};   // normals are flipped to face this point (sensor origin)
};

struct PlaneReport {
  Plane plane;
  size_t inliers = 0;
  double inlierRatio = 0;     // inliers / cleaned cloud size
  int ransacIters = 0;
  FlatnessStats flat;
  std::vector<size_t> bandIdx;  // indices into `cleaned` used for flatness stats
};

struct Timing {
  double voxelMs = 0, sorMs = 0, ransacMs = 0, flatMs = 0, totalMs = 0;
};

struct PipelineOutput {
  Cloud cleaned;  // cloud after downsampling + outlier removal (all indices refer to this)
  std::vector<PlaneReport> planes;
  std::vector<int> label;       // per cleaned point: plane index, or -1
  std::vector<double> deviation;  // per cleaned point: signed distance to its plane (NaN if unlabeled)
  size_t nIn = 0, nVoxel = 0, nClean = 0;
  Timing t;
};

inline PipelineOutput runPipeline(const Cloud& raw, const PipelineParams& prm) {
  using clk = std::chrono::steady_clock;
  auto ms = [](clk::time_point a, clk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
  PipelineOutput out;
  out.nIn = raw.size();
  const auto t0 = clk::now();

  Cloud vox = voxelDownsample(raw, prm.voxelLeaf);
  out.nVoxel = vox.size();
  const auto t1 = clk::now();

  out.cleaned = statisticalOutlierRemoval(vox, prm.sorK, prm.sorStd);
  out.nClean = out.cleaned.size();
  const auto t2 = clk::now();

  auto found = extractPlanes(out.cleaned, prm.ransac, prm.maxPlanes);
  const auto t3 = clk::now();

  out.label.assign(out.cleaned.size(), -1);
  out.deviation.assign(out.cleaned.size(), std::nan(""));
  std::vector<int8_t> owned(out.cleaned.size(), 0);
  for (size_t k = 0; k < found.size(); ++k)
    for (size_t i : found[k].inliers) {
      out.label[i] = static_cast<int>(k);
      owned[i] = 1;
    }

  for (size_t k = 0; k < found.size(); ++k) {
    PlaneReport rep;
    rep.plane = found[k].plane;
    orientToward(rep.plane, prm.viewpoint);
    rep.inliers = found[k].inliers.size();
    rep.inlierRatio = static_cast<double>(rep.inliers) / std::max<size_t>(1, out.cleaned.size());
    rep.ransacIters = found[k].iterations;
    // "owned by OTHER planes" mask for this plane
    std::vector<int8_t> other(owned.size(), 0);
    for (size_t i = 0; i < owned.size(); ++i) other[i] = (out.label[i] >= 0 && out.label[i] != static_cast<int>(k));
    rep.bandIdx = surfaceBand(out.cleaned, rep.plane, found[k].inliers, other, prm.flatBand);
    rep.flat = flatnessStats(out.cleaned, rep.plane, rep.bandIdx, std::max(prm.flatCell, 2.0 * prm.voxelLeaf));
    for (size_t i : rep.bandIdx) {
      out.deviation[i] = rep.plane.dist(out.cleaned[i]);
      if (out.label[i] < 0) out.label[i] = static_cast<int>(k);
    }
    out.planes.push_back(std::move(rep));
  }
  const auto t4 = clk::now();
  out.t = {ms(t0, t1), ms(t1, t2), ms(t2, t3), ms(t3, t4), ms(t0, t4)};
  return out;
}

// Offset between nearly-parallel planes (e.g. a raised patch detected as its own plane in front of a wall).
// offsetM = mean signed distance of plane j's points to plane i; > 0 means plane j lies on the normal side
// of plane i, i.e. closer to the sensor.  Averaging over a whole plane's points makes this much less
// noisy than comparing the two plane equations directly.
struct ParallelOffset {
  size_t i = 0, j = 0;
  double angleDeg = 0;
  double offsetM = 0;
};

inline std::vector<ParallelOffset> parallelOffsets(const PipelineOutput& out, double maxAngleDeg = 4.0) {
  std::vector<ParallelOffset> res;
  const double cosMax = std::cos(maxAngleDeg * 3.14159265358979323846 / 180.0);
  for (size_t i = 0; i < out.planes.size(); ++i)
    for (size_t j = i + 1; j < out.planes.size(); ++j) {
      const double c = dot(out.planes[i].plane.n, out.planes[j].plane.n);
      if (c < cosMax) continue;
      double s = 0;
      size_t cnt = 0;
      for (size_t k = 0; k < out.cleaned.size(); ++k)
        if (out.label[k] == static_cast<int>(j)) {
          s += out.planes[i].plane.dist(out.cleaned[k]);
          ++cnt;
        }
      if (cnt == 0) continue;
      res.push_back({i, j, std::acos(std::min(1.0, c)) * 180.0 / 3.14159265358979323846, s / cnt});
    }
  return res;
}

}  // namespace planeperc
