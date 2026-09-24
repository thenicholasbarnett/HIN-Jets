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
// JME: JER twiki / CMSSW SmearedJetProducerT matching -- each reco jet takes
// the closest gen jet with dR < R/2 and |pT - pT_gen| < 3 sigma_JER pT
// (JEC-corrected pT, sigma_JER per reco jet, e.g. JetSmearer::Resolution):
//
//   std::vector<double> sigma(nref);
//   for (int i = 0; i < nref; i++) {
//     sigma[i] = smearer.Resolution(ptCorr[i], jteta[i], rho);
//   }
//   JetMatcher::Match(nref, jteta, jtphi, ngen, geneta, genphi, 0.4,
//                     JetMatcher::Mode::JME, 0.5, ptCorr, genpt,
//                     sigma.data());
//
// Composes with JetSorter: match[order[0]] is the leading jet's gen jet

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace JetMatcher {

enum class Mode { Nearest, Unique, JME };

inline double DeltaR(double eta1, double phi1, double eta2, double phi2) {
  const double dphi = std::remainder(phi1 - phi2, 2 * M_PI);
  const double deta = eta1 - eta2;
  return std::sqrt(deta * deta + dphi * dphi);
}

template <typename T>
std::vector<int> Match(int nref, const T *jteta, const T *jtphi, int ngen,
                       const T *geneta, const T *genphi, double coneR,
                       Mode mode = Mode::Unique, double dRFraction = 0.5,
                       const T *jtpt = nullptr, const T *genpt = nullptr,
                       const double *sigmaJER = nullptr, double nSigma = 3.0) {
  const double maxDR = dRFraction * coneR;
  std::vector<int> match(nref, -1);
  if (mode == Mode::JME && (!jtpt || !genpt || !sigmaJER)) {
    throw std::invalid_argument(
        "JetMatcher: Mode::JME needs jtpt, genpt and sigmaJER");
  }
  if (mode == Mode::Nearest || mode == Mode::JME) {
    for (int i = 0; i < nref; i++) {
      double best = maxDR;
      for (int g = 0; g < ngen; g++) {
        if (mode == Mode::JME &&
            std::abs(jtpt[i] - genpt[g]) >= nSigma * sigmaJER[i] * jtpt[i]) {
          continue;
        }
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
