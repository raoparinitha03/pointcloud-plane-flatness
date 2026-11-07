// Synthetic evaluation harness: builds a noisy room (floor + 2 walls + one 3 cm bump + outliers)
// with KNOWN ground truth, runs the pipeline, and reports accuracy / speed / parameter sweeps.
//   ./synthetic_eval [outdir]      -> writes sweep.csv and flatness_wall.ply into outdir
#include <cstdio>
#include <fstream>
#include <random>
#include <string>
#include "planeperc/pipeline.hpp"

using namespace planeperc;

namespace {

constexpr double kRoomX = 5.0, kRoomY = 4.0, kRoomZ = 2.6;
constexpr double kBumpAmp = 0.03, kBumpSigma = 0.15, kBumpY = 2.0, kBumpZ = 1.3;
constexpr double kNoiseSigma = 0.004;  // 4 mm, typical depth-camera noise at ~2 m

struct GtPlane {
  const char* name;
  Plane plane;
};
const GtPlane kGt[3] = {{"floor  z=0", {{0, 0, 1}, 0}}, {"wallA  x=0", {{1, 0, 0}, 0}}, {"wallB  y=0", {{0, 1, 0}, 0}}};

Cloud makeRoom(uint64_t seed, double density, double outlierFrac) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> U01(0, 1);
  std::normal_distribution<double> noise(0, kNoiseSigma);
  Cloud c;
  auto n = [&](double area) { return static_cast<size_t>(area * density); };
  for (size_t i = 0; i < n(kRoomX * kRoomY); ++i) c.push_back({U01(rng) * kRoomX, U01(rng) * kRoomY, noise(rng)});  // floor
  for (size_t i = 0; i < n(kRoomY * kRoomZ); ++i) {  // wall A with a bump
    const double y = U01(rng) * kRoomY, z = U01(rng) * kRoomZ;
    const double r2 = (y - kBumpY) * (y - kBumpY) + (z - kBumpZ) * (z - kBumpZ);
    c.push_back({noise(rng) + kBumpAmp * std::exp(-r2 / (2 * kBumpSigma * kBumpSigma)), y, z});
  }
  for (size_t i = 0; i < n(kRoomX * kRoomZ); ++i) c.push_back({U01(rng) * kRoomX, noise(rng), U01(rng) * kRoomZ});  // wall B
  const size_t nOut = static_cast<size_t>(c.size() * outlierFrac);
  for (size_t i = 0; i < nOut; ++i)  // uniform clutter / flying pixels
    c.push_back({-0.5 + U01(rng) * (kRoomX + 1), -0.5 + U01(rng) * (kRoomY + 1), -0.2 + U01(rng) * (kRoomZ + 0.4)});
  std::shuffle(c.begin(), c.end(), rng);
  return c;
}

struct Match {
  int gt = -1;
  double angleDeg = 0;
  double offsetErr = 0;
};
Match matchToGt(const Plane& p) {
  Match m;
  double bestDot = -1;
  for (int g = 0; g < 3; ++g) {
    const double d = std::fabs(dot(p.n, kGt[g].plane.n));
    if (d > bestDot) {
      bestDot = d;
      m.gt = g;
    }
  }
  m.angleDeg = std::acos(std::min(1.0, bestDot)) * 180.0 / M_PI;
  // compare offsets after flipping sign to the same hemisphere
  const double s = dot(p.n, kGt[m.gt].plane.n) >= 0 ? 1.0 : -1.0;
  m.offsetErr = std::fabs(s * p.d - kGt[m.gt].plane.d);
  return m;
}

void writePly(const std::string& path, const PipelineOutput& o, int planeIdx, double range) {
  std::ofstream f(path);
  std::vector<size_t> pts;
  for (size_t i = 0; i < o.cleaned.size(); ++i)
    if (o.label[i] == planeIdx && !std::isnan(o.deviation[i])) pts.push_back(i);
  f << "ply\nformat ascii 1.0\nelement vertex " << pts.size()
    << "\nproperty float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nend_header\n";
  for (size_t i : pts) {
    uint8_t r, g, b;
    deviationColor(o.deviation[i], range, r, g, b);
    f << o.cleaned[i].x << ' ' << o.cleaned[i].y << ' ' << o.cleaned[i].z << ' ' << int(r) << ' ' << int(g) << ' ' << int(b) << '\n';
  }
}

}  // namespace

int main(int argc, char** argv) {
  const std::string outdir = argc > 1 ? argv[1] : ".";
  const Vec3 roomCentre{kRoomX / 2, kRoomY / 2, kRoomZ / 2};

  // ---------- 1. single detailed run ----------
  {
    Cloud raw = makeRoom(0, 3000.0, 0.03);
    PipelineParams prm;
    prm.viewpoint = roomCentre;
    prm.voxelLeaf = 0.03;
    prm.ransac.threshold = 0.01;
    PipelineOutput out = runPipeline(raw, prm);

    std::printf("=== Single run: %zu raw pts -> %zu after voxel(%.0f mm) -> %zu after SOR ===\n", out.nIn, out.nVoxel,
                prm.voxelLeaf * 1000, out.nClean);
    std::printf("time: voxel %.1f ms | SOR %.1f ms | RANSAC %.1f ms | flatness %.1f ms | total %.1f ms\n", out.t.voxelMs,
                out.t.sorMs, out.t.ransacMs, out.t.flatMs, out.t.totalMs);
    int bumpPlane = -1;
    for (size_t k = 0; k < out.planes.size(); ++k) {
      const auto& r = out.planes[k];
      const Match m = matchToGt(r.plane);
      std::printf("plane %zu ~ %s | normal err %.3f deg | offset err %.2f mm | inlier ratio %.1f%% | RANSAC iters %d\n", k,
                  kGt[m.gt].name, m.angleDeg, m.offsetErr * 1000, r.inlierRatio * 100, r.ransacIters);
      std::printf("          flatness: rms %.2f mm | p95 %.2f mm | cell-median peak %+.2f mm | valley %+.2f mm (n=%zu)\n", r.flat.rms * 1000,
                  r.flat.p95Abs * 1000, r.flat.maxPos * 1000, r.flat.maxNeg * 1000, r.flat.n);
      if (m.gt == 1) bumpPlane = static_cast<int>(k);
    }
    std::printf("ground truth: bump on wallA = %.1f mm peak (toward room) | sensor noise sigma = %.1f mm\n", kBumpAmp * 1000,
                kNoiseSigma * 1000);
    if (bumpPlane >= 0) writePly(outdir + "/flatness_wall.ply", out, bumpPlane, 0.03);
  }

  // ---------- 2. parameter sweep ----------
  std::ofstream csv(outdir + "/sweep.csv");
  csv << "voxel_mm,ransac_thr_mm,mean_normal_err_deg,planes_found,mean_inlier_ratio,bump_est_mm,mean_total_ms\n";
  std::printf("\n=== Sweep (3 seeds each) ===\n");
  std::printf("voxel_mm thr_mm | normal_err_deg planes inlier%% bump_est_mm  time_ms\n");
  const double voxels[] = {0.0, 0.03, 0.05, 0.10};
  const double thrs[] = {0.005, 0.01, 0.02, 0.05};
  for (double v : voxels)
    for (double thr : thrs) {
      double errSum = 0, ratioSum = 0, bumpSum = 0, msSum = 0;
      int errCnt = 0, planesSum = 0, bumpCnt = 0;
      for (int seed = 0; seed < 3; ++seed) {
        Cloud raw = makeRoom(100 + seed, 3000.0, 0.03);
        PipelineParams prm;
        prm.viewpoint = roomCentre;
        prm.voxelLeaf = v;
        prm.ransac.threshold = thr;
        prm.ransac.seed = 1000 + seed;
        prm.flatBand = std::max(0.05, 3 * thr);
        PipelineOutput o = runPipeline(raw, prm);
        planesSum += static_cast<int>(o.planes.size());
        msSum += o.t.totalMs;
        // For each ground-truth plane keep only the LARGEST matching detection; any additional
        // detection (e.g. a patch fitted to the bump) counts as an "extra" plane, not as an error.
        size_t bestInl[3] = {0, 0, 0};
        const PlaneReport* bestRep[3] = {nullptr, nullptr, nullptr};
        double bestAng[3] = {0, 0, 0};
        for (const auto& r : o.planes) {
          const Match m = matchToGt(r.plane);
          if (r.inliers > bestInl[m.gt]) {
            bestInl[m.gt] = r.inliers;
            bestRep[m.gt] = &r;
            bestAng[m.gt] = m.angleDeg;
          }
        }
        for (int g = 0; g < 3; ++g) {
          if (!bestRep[g]) continue;
          errSum += bestAng[g];
          ++errCnt;
          ratioSum += bestRep[g]->inlierRatio;
          if (g == 1) {
            bumpSum += bestRep[g]->flat.maxPos;
            ++bumpCnt;
          }
        }
      }
      const double meanErr = errCnt ? errSum / errCnt : -1;
      const double meanRatio = errCnt ? ratioSum / errCnt : 0;
      const double meanBump = bumpCnt ? bumpSum / bumpCnt : -1;
      std::printf("%7.0f %6.0f | %14.3f %6.1f %7.1f %12.1f %8.1f\n", v * 1000, thr * 1000, meanErr, planesSum / 3.0,
                  meanRatio * 100, meanBump * 1000, msSum / 3);
      csv << v * 1000 << ',' << thr * 1000 << ',' << meanErr << ',' << planesSum / 3.0 << ',' << meanRatio << ','
          << meanBump * 1000 << ',' << msSum / 3 << '\n';
    }
  return 0;
}
