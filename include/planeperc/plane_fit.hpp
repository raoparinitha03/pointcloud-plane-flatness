#pragma once
// Least-squares plane fit from scratch.
// Idea: the best-fit plane passes through the centroid, and its normal is the direction of
// LEAST variance of the points = eigenvector of the 3x3 covariance matrix with the SMALLEST
// eigenvalue.  We find eigenvectors with the classic Jacobi rotation method (no library).
#include <cmath>
#include "planeperc/types.hpp"

namespace planeperc {

// Cyclic Jacobi eigen-decomposition of a symmetric 3x3 matrix.
// On return: w[i] = eigenvalues, column i of V = the matching unit eigenvector.
inline void jacobiEigen3(double A[3][3], double w[3], double V[3][3]) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) V[i][j] = (i == j) ? 1.0 : 0.0;

  for (int sweep = 0; sweep < 50; ++sweep) {
    const double off = std::fabs(A[0][1]) + std::fabs(A[0][2]) + std::fabs(A[1][2]);
    if (off < 1e-18) break;
    const int P[3] = {0, 0, 1}, Q[3] = {1, 2, 2};
    for (int r = 0; r < 3; ++r) {
      const int p = P[r], q = Q[r];
      if (std::fabs(A[p][q]) < 1e-300) continue;
      // choose rotation angle that zeroes A[p][q]
      const double theta = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
      const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
      const double c = 1.0 / std::sqrt(t * t + 1.0);
      const double s = t * c;
      // A <- A * J   (columns p,q)
      for (int k = 0; k < 3; ++k) {
        const double akp = A[k][p], akq = A[k][q];
        A[k][p] = c * akp - s * akq;
        A[k][q] = s * akp + c * akq;
      }
      // A <- J^T * A (rows p,q)
      for (int k = 0; k < 3; ++k) {
        const double apk = A[p][k], aqk = A[q][k];
        A[p][k] = c * apk - s * aqk;
        A[q][k] = s * apk + c * aqk;
      }
      // accumulate eigenvectors  V <- V * J
      for (int k = 0; k < 3; ++k) {
        const double vkp = V[k][p], vkq = V[k][q];
        V[k][p] = c * vkp - s * vkq;
        V[k][q] = s * vkp + c * vkq;
      }
    }
  }
  for (int i = 0; i < 3; ++i) w[i] = A[i][i];
}

// Fit a plane to the points `cloud[idx[i]]` by least squares (orthogonal distance).
// Returns false if there are fewer than 3 points.
inline bool fitPlaneLS(const Cloud& cloud, const std::vector<size_t>& idx, Plane& out) {
  if (idx.size() < 3) return false;
  Vec3 c{0, 0, 0};
  for (size_t i : idx) c = c + cloud[i];
  c = c * (1.0 / idx.size());

  double C[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};  // covariance (unnormalised)
  for (size_t i : idx) {
    const Vec3 d = cloud[i] - c;
    const double v[3] = {d.x, d.y, d.z};
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) C[a][b] += v[a] * v[b];
  }
  double w[3], V[3][3];
  jacobiEigen3(C, w, V);
  int m = 0;  // index of smallest eigenvalue
  for (int i = 1; i < 3; ++i)
    if (w[i] < w[m]) m = i;
  Vec3 n{V[0][m], V[1][m], V[2][m]};
  n = n * (1.0 / norm(n));
  out.n = n;
  out.d = -dot(n, c);
  return true;
}

// Flip the plane so its normal points TOWARD `viewpoint` (sign of a plane is ambiguous).
inline void orientToward(Plane& pl, const Vec3& viewpoint) {
  if (pl.dist(viewpoint) < 0) {
    pl.n = pl.n * -1.0;
    pl.d = -pl.d;
  }
}

}  // namespace planeperc
