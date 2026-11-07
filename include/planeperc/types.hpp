#pragma once
// Minimal 3D math types. Everything in this library is written from scratch:
// no PCL, no Open3D, no Eigen.
#include <cmath>
#include <cstddef>
#include <vector>

namespace planeperc {

struct Vec3 {
  double x = 0, y = 0, z = 0;
};

inline Vec3 operator+(Vec3 a, Vec3 b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(Vec3 a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline double dot(Vec3 a, Vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(Vec3 a, Vec3 b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(Vec3 a) { return std::sqrt(dot(a, a)); }
inline double coord(const Vec3& p, int axis) { return axis == 0 ? p.x : (axis == 1 ? p.y : p.z); }

using Cloud = std::vector<Vec3>;

// Plane in Hessian normal form:  n . p + d = 0,  |n| = 1.
// dist(p) is the SIGNED perpendicular distance from p to the plane.
struct Plane {
  Vec3 n{0, 0, 1};
  double d = 0;
  double dist(const Vec3& p) const { return dot(n, p) + d; }
};

}  // namespace planeperc
