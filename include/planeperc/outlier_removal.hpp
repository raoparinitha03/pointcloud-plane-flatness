#pragma once
// Statistical outlier removal (SOR).
// For each point: mean distance to its k nearest neighbours.  Over the whole cloud these
// mean distances have mean mu and std sigma.  Drop points with meanDist > mu + stdMult*sigma.
// Isolated "flying pixels" have large neighbour distances, so they get removed.
#include <cmath>
#include "planeperc/kdtree.hpp"

namespace planeperc {

inline Cloud statisticalOutlierRemoval(const Cloud& in, size_t k, double stdMult) {
  if (in.size() <= k + 1) return in;
  KdTree tree(in);
  std::vector<double> meanDist(in.size());
  for (size_t i = 0; i < in.size(); ++i) {
    auto nn = tree.knn(in[i], k + 1);  // +1 because the nearest neighbour is the point itself
    double s = 0;
    size_t cnt = 0;
    for (size_t j = 1; j < nn.size(); ++j) {  // skip j==0 (self)
      s += std::sqrt(nn[j].first);
      ++cnt;
    }
    meanDist[i] = cnt ? s / cnt : 0.0;
  }
  double mu = 0;
  for (double v : meanDist) mu += v;
  mu /= meanDist.size();
  double var = 0;
  for (double v : meanDist) var += (v - mu) * (v - mu);
  const double sigma = std::sqrt(var / meanDist.size());
  const double thr = mu + stdMult * sigma;

  Cloud out;
  out.reserve(in.size());
  for (size_t i = 0; i < in.size(); ++i)
    if (meanDist[i] <= thr) out.push_back(in[i]);
  return out;
}

}  // namespace planeperc
