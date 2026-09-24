#ifndef JETSORTER_H
#define JETSORTER_H

// JetSorter v1.0
// pT ordering of forest-style jet arrays with just this header
// Author: Nicholas Shawn Barnett

// USAGE
// JetSorter::Order(n, pt) returns an n-long index map, highest pT first:
// order[i] is the original index of the i-th hardest jet, so any jet array
// is read in pT order through it
//
//   std::vector<int> order = JetSorter::Order(nref, ptCorr);
//   float leadEta = jteta[order[0]];
//   for (int i = 0; i < nref; i++) {
//     int j = order[i];
//     // ptCorr[j], jteta[j], jtphi[j], jtPfCHF[j], ...
//   }
//
// Reorder copies an array into that order, when a sorted copy is handier:
//
//   std::vector<float> etaSorted = JetSorter::Reorder(order, jteta);

#include <algorithm>
#include <numeric>
#include <vector>

namespace JetSorter {

// equal pT keep their original order
template <typename T> std::vector<int> Order(int n, const T *pt) {
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(),
                   [&](int a, int b) { return pt[a] > pt[b]; });
  return order;
}

template <typename T>
std::vector<T> Reorder(const std::vector<int> &order, const T *x) {
  std::vector<T> out(order.size());
  for (size_t i = 0; i < order.size(); i++) {
    out[i] = x[order[i]];
  }
  return out;
}

} // namespace JetSorter

#endif
