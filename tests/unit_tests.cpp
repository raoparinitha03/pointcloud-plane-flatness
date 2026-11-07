// Plain-assert unit tests (no framework needed) so CI can run them anywhere.
#include <cassert>
#include <cstdio>
#include <random>
#include "planeperc/pipeline.hpp"

using namespace planeperc;

static int g_fail = 0;
#define CHECK(cond)                                                      \
  do {                                                                   \
    if (!(cond)) {                                                       \
      std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);        \
      ++g_fail;                                                          \
    }                                                                    \
  } while (0)

static void testVoxel() {
  Cloud c = {{0.01, 0.01, 0.01}, {0.03, 0.03, 0.03}, {1.01, 0, 0}, {-0.01, -0.01, -0.01}};
  Cloud o = voxelDownsample(c, 0.1);
  CHECK(o.size() == 3);  // two points share a voxel, one is in a different x-voxel, one negative voxel
  CHECK(voxelDownsample(c, 0.0).size() == 4);
  // centroid correctness
  Cloud two = {{0.0, 0.0, 0.0}, {0.04, 0.0, 0.0}};
  Cloud m = voxelDownsample(two, 0.1);
  CHECK(m.size() == 1 && std::fabs(m[0].x - 0.02) < 1e-12);
}

static void testKdTree() {
  std::mt19937 rng(1);
  std::uniform_real_distribution<double> U(-5, 5);
  Cloud pts(2000);
  for (auto& p : pts) p = {U(rng), U(rng), U(rng)};
  KdTree tree(pts);
  for (int t = 0; t < 50; ++t) {
    Vec3 q{U(rng), U(rng), U(rng)};
    auto got = tree.knn(q, 10);
    std::vector<std::pair<double, size_t>> brute;
    for (size_t i = 0; i < pts.size(); ++i) {
      Vec3 d = pts[i] - q;
      brute.push_back({dot(d, d), i});
    }
    std::sort(brute.begin(), brute.end());
    CHECK(got.size() == 10);
    for (int j = 0; j < 10; ++j) CHECK(std::fabs(got[j].first - brute[j].first) < 1e-12);
  }
}

static void testJacobi() {
  double A[3][3] = {{4, 1, 0}, {1, 3, 0}, {0, 0, 2}};
  double w[3], V[3][3];
  jacobiEigen3(A, w, V);
  std::vector<double> ev = {w[0], w[1], w[2]};
  std::sort(ev.begin(), ev.end());
  const double r = std::sqrt(1.25);  // eigenvalues 3.5 -/+ sqrt(1.25)
  CHECK(std::fabs(ev[0] - 2.0) < 1e-9);
  CHECK(std::fabs(ev[1] - (3.5 - r)) < 1e-9);
  CHECK(std::fabs(ev[2] - (3.5 + r)) < 1e-9);
}

static void testPlaneFit() {
  // points on the plane z = 0.5x + 0.2y + 1
  std::mt19937 rng(2);
  std::uniform_real_distribution<double> U(-2, 2);
  Cloud c;
  std::vector<size_t> idx;
  for (int i = 0; i < 500; ++i) {
    double x = U(rng), y = U(rng);
    c.push_back({x, y, 0.5 * x + 0.2 * y + 1.0});
    idx.push_back(i);
  }
  Plane p;
  CHECK(fitPlaneLS(c, idx, p));
  Vec3 expect{-0.5, -0.2, 1.0};
  expect = expect * (1.0 / norm(expect));
  CHECK(std::fabs(std::fabs(dot(p.n, expect)) - 1.0) < 1e-9);
  for (auto& q : c) CHECK(std::fabs(p.dist(q)) < 1e-9);
}

static void testRansacWithOutliers() {
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> U(0, 4), Un(-0.003, 0.003), Uo(-1, 5);
  Cloud c;
  for (int i = 0; i < 3000; ++i) c.push_back({U(rng), U(rng), Un(rng)});  // plane z~0
  for (int i = 0; i < 1500; ++i) c.push_back({Uo(rng), Uo(rng), Uo(rng)});  // 33% junk
  RansacParams prm;
  prm.threshold = 0.01;
  auto res = ransacPlane(c, prm);
  CHECK(res.has_value());
  if (res) {
    CHECK(std::fabs(std::fabs(res->plane.n.z) - 1.0) < 1e-3);
    CHECK(res->inliers.size() > 2900);
  }
}

static void testSor() {
  std::mt19937 rng(4);
  std::uniform_real_distribution<double> U(0, 1), Uo(5, 20);
  Cloud c;
  for (int i = 0; i < 2000; ++i) c.push_back({U(rng), U(rng), U(rng)});
  for (int i = 0; i < 20; ++i) c.push_back({Uo(rng), Uo(rng), Uo(rng)});
  Cloud o = statisticalOutlierRemoval(c, 16, 1.5);
  size_t far = 0;
  for (auto& p : o)
    if (p.x > 3) ++far;
  CHECK(far == 0);
  CHECK(o.size() > 1500);
}


static void testParallelOffset() {
  // wall at x = 0 (normal +x) with a 30 mm raised patch at x = 0.03; sensor on the +x side
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> U(0, 1);
  std::normal_distribution<double> N(0, 0.005);
  Cloud c;
  for (int i = 0; i < 6000; ++i) {
    const double y = U(rng) * 1.6, z = U(rng) * 1.0;
    const bool inPatch = (y > 0.5 && y < 1.1 && z > 0.1 && z < 0.5);
    c.push_back({(inPatch ? 0.03 : 0.0) + N(rng), y, z});
  }
  PipelineParams prm;
  prm.voxelLeaf = 0.0;
  prm.ransac.threshold = 0.015;
  prm.ransac.minInliers = 60;
  prm.flatCell = 0.1;
  prm.viewpoint = {1.5, 0.8, 0.5};
  PipelineOutput out = runPipeline(c, prm);
  CHECK(out.planes.size() >= 2);
  auto offs = parallelOffsets(out);
  bool found = false;
  for (const auto& o : offs)
    if (std::fabs(std::fabs(o.offsetM) - 0.03) < 0.006) found = true;
  CHECK(found);
  for (const auto& o : offs) std::printf("  parallel pair %zu-%zu: angle %.2f deg, offset %+.1f mm\n", o.i, o.j, o.angleDeg, o.offsetM * 1000);
}

int main() {
  testVoxel();
  testKdTree();
  testJacobi();
  testPlaneFit();
  testRansacWithOutliers();
  testSor();
  testParallelOffset();
  if (g_fail == 0) std::printf("ALL TESTS PASSED\n");
  return g_fail == 0 ? 0 : 1;
}
