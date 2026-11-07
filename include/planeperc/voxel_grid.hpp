#pragma once
// Voxel-grid downsampling: bucket points into cubes of side `leaf`, replace each
// occupied cube by the centroid of its points.  O(n) with a hash map.
// Effect: uniform density + averaging reduces sensor noise by ~sqrt(points per voxel).
#include <cstdint>
#include <unordered_map>
#include "planeperc/types.hpp"

namespace planeperc {

inline Cloud voxelDownsample(const Cloud& in, double leaf) {
  if (leaf <= 0.0) return in;  // leaf == 0 means "disabled"

  struct Key {
    int32_t i, j, k;
    bool operator==(const Key& o) const { return i == o.i && j == o.j && k == o.k; }
  };
  struct KeyHash {
    size_t operator()(const Key& key) const {
      // three large primes (Teschner et al. spatial hashing)
      uint64_t h = static_cast<uint64_t>(static_cast<uint32_t>(key.i)) * 73856093ull ^
                   static_cast<uint64_t>(static_cast<uint32_t>(key.j)) * 19349663ull ^
                   static_cast<uint64_t>(static_cast<uint32_t>(key.k)) * 83492791ull;
      return static_cast<size_t>(h);
    }
  };
  struct Acc {
    double sx = 0, sy = 0, sz = 0;
    uint32_t n = 0;
  };

  std::unordered_map<Key, Acc, KeyHash> cells;
  cells.reserve(in.size() / 2 + 1);
  const double inv = 1.0 / leaf;
  for (const Vec3& p : in) {
    // floor (not truncation!) so negative coordinates bucket correctly
    Key key{static_cast<int32_t>(std::floor(p.x * inv)), static_cast<int32_t>(std::floor(p.y * inv)),
            static_cast<int32_t>(std::floor(p.z * inv))};
    Acc& a = cells[key];
    a.sx += p.x;
    a.sy += p.y;
    a.sz += p.z;
    ++a.n;
  }
  Cloud out;
  out.reserve(cells.size());
  for (const auto& kv : cells) {
    const Acc& a = kv.second;
    out.push_back({a.sx / a.n, a.sy / a.n, a.sz / a.n});
  }
  return out;
}

}  // namespace planeperc
