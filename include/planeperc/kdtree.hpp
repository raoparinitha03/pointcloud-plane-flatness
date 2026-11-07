#pragma once
// 3D k-d tree for k-nearest-neighbour queries.
// Build: recursive median split on the axis of widest spread (nth_element) -> O(n log n).
// Query: depth-first, visit the far child only if the splitting plane is closer than the
//        current k-th best distance -> ~O(log n) per query for low dimensions.
// NOTE: the tree stores a REFERENCE to the cloud; the cloud must outlive the tree.
#include <algorithm>
#include <numeric>
#include <queue>
#include <utility>
#include "planeperc/types.hpp"

namespace planeperc {

class KdTree {
 public:
  explicit KdTree(const Cloud& pts) : pts_(pts) {
    idx_.resize(pts.size());
    std::iota(idx_.begin(), idx_.end(), size_t{0});
    nodes_.reserve(pts.size());
    root_ = build(0, idx_.size());
  }

  // Returns up to k (squared_distance, index) pairs sorted by ascending distance.
  std::vector<std::pair<double, size_t>> knn(const Vec3& q, size_t k) const {
    std::priority_queue<std::pair<double, size_t>> heap;  // max-heap: worst candidate on top
    if (k > 0) search(root_, q, k, heap);
    std::vector<std::pair<double, size_t>> out(heap.size());
    for (size_t i = out.size(); i > 0; --i) {
      out[i - 1] = heap.top();
      heap.pop();
    }
    return out;
  }

 private:
  struct Node {
    int axis;
    int left, right;  // indices into nodes_, -1 = none
    size_t point;     // index into pts_
  };

  int build(size_t lo, size_t hi) {
    if (lo >= hi) return -1;
    // choose the axis with the largest extent
    Vec3 mn = pts_[idx_[lo]], mx = mn;
    for (size_t i = lo; i < hi; ++i) {
      const Vec3& p = pts_[idx_[i]];
      mn = {std::min(mn.x, p.x), std::min(mn.y, p.y), std::min(mn.z, p.z)};
      mx = {std::max(mx.x, p.x), std::max(mx.y, p.y), std::max(mx.z, p.z)};
    }
    const double ex = mx.x - mn.x, ey = mx.y - mn.y, ez = mx.z - mn.z;
    const int axis = (ex >= ey && ex >= ez) ? 0 : (ey >= ez ? 1 : 2);
    const size_t mid = lo + (hi - lo) / 2;
    std::nth_element(idx_.begin() + lo, idx_.begin() + mid, idx_.begin() + hi,
                     [&](size_t a, size_t b) { return coord(pts_[a], axis) < coord(pts_[b], axis); });
    const int me = static_cast<int>(nodes_.size());
    nodes_.push_back({axis, -1, -1, idx_[mid]});
    const int l = build(lo, mid);
    const int r = build(mid + 1, hi);
    nodes_[me].left = l;
    nodes_[me].right = r;
    return me;
  }

  void search(int ni, const Vec3& q, size_t k, std::priority_queue<std::pair<double, size_t>>& heap) const {
    if (ni < 0) return;
    const Node& nd = nodes_[ni];
    const Vec3& p = pts_[nd.point];
    const Vec3 dv = q - p;
    const double d2 = dot(dv, dv);
    if (heap.size() < k) {
      heap.push({d2, nd.point});
    } else if (d2 < heap.top().first) {
      heap.pop();
      heap.push({d2, nd.point});
    }
    const double diff = coord(q, nd.axis) - coord(p, nd.axis);
    const int nearChild = diff < 0 ? nd.left : nd.right;
    const int farChild = diff < 0 ? nd.right : nd.left;
    search(nearChild, q, k, heap);
    if (heap.size() < k || diff * diff < heap.top().first) search(farChild, q, k, heap);
  }

  const Cloud& pts_;
  std::vector<size_t> idx_;
  std::vector<Node> nodes_;
  int root_ = -1;
};

}  // namespace planeperc
