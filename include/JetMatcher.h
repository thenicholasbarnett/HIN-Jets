#ifndef JETMATCHER_H
#define JETMATCHER_H

// JetMatcher v1.0
// reco to gen jet matching of forest-style jet arrays with just this header
// Author: Nicholas Shawn Barnett

// USAGE
// JetMatcher::Match(...) returns an nref-long index map, reco to gen:
// match[i] is the gen jet matched to reco jet i, or -1 for none, within
// dR < R/2 for cone radius R (JME convention)
//
//   std::vector<int> match = JetMatcher::Match(nref, jteta, jtphi,
//                                              ngen, geneta, genphi, 0.4);
//   int g = match[i];
//   double genPt = (g >= 0) ? genpt[g] : -1; // -1 = unmatched, for JetSmearer
//
// Unique (default): one to one, closest pairs first, so no jet is matched
// twice. Nearest: each reco jet takes its closest gen jet, which two reco
// jets can share. dRFraction sets the limit as a fraction of R; the forest's
// ref* matching is Unique with dR < R:
//
//   JetMatcher::Match(nref, jteta, jtphi, ngen, geneta, genphi, 0.4,
//                     JetMatcher::Mode::Unique, 1.0);
//
// Composes with JetSorter: match[order[0]] is the leading jet's gen jet

#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>

namespace JetMatcher {

enum class Mode { Nearest, Unique };

inline double DeltaR(double eta1, double phi1, double eta2, double phi2) {
  const double dphi = std::remainder(phi1 - phi2, 2 * M_PI);
  const double deta = eta1 - eta2;
  return std::sqrt(deta * deta + dphi * dphi);
}

template <typename T>
std::vector<int> Match(int nref, const T *jteta, const T *jtphi, int ngen,
                       const T *geneta, const T *genphi, double coneR,
                       Mode mode = Mode::Unique, double dRFraction = 0.5) {
  const double maxDR = dRFraction * coneR;
  std::vector<int> match(nref, -1);
  if (mode == Mode::Nearest) {
    for (int i = 0; i < nref; i++) {
      double best = maxDR;
      for (int g = 0; g < ngen; g++) {
        const double dr = DeltaR(jteta[i], jtphi[i], geneta[g], genphi[g]);
        if (dr < best) {
          best = dr;
          match[i] = g;
        }
      }
    }
    return match;
  }

  // unique: every pair within maxDR, smallest dR first
  std::vector<std::tuple<double, int, int>> pairs;
  for (int i = 0; i < nref; i++) {
    for (int g = 0; g < ngen; g++) {
      const double dr = DeltaR(jteta[i], jtphi[i], geneta[g], genphi[g]);
      if (dr < maxDR) {
        pairs.emplace_back(dr, i, g);
      }
    }
  }
  std::sort(pairs.begin(), pairs.end());
  std::vector<bool> genUsed(ngen, false);
  for (const auto &p : pairs) {
    const int i = std::get<1>(p);
    const int g = std::get<2>(p);
    if (match[i] < 0 && !genUsed[g]) {
      match[i] = g;
      genUsed[g] = true;
    }
  }
  return match;
}

} // namespace JetMatcher

#endif
