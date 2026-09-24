#ifndef JETMAPPER_H
#define JETMAPPER_H

// JetMapper v1.0
// index maps over forest-style jet arrays with just this header: pT order,
// and reco to gen matching
// Author: Nicholas Shawn Barnett

// USAGE
// Every map is an index array, read any jet array through it
//
// Order(n, pt): n-long, highest pT first -- order[i] is the original index
// of the i-th hardest jet
//
//   std::vector<int> order = JetMapper::Order(nref, ptCorr);
//   float leadEta = jteta[order[0]];
//
// Match(...): nref-long, reco to gen -- match[i] is the gen jet matched to
// reco jet i, or JetMapper::kUnmatched (-999, as the forest's refpt) for
// none. Check it before indexing, jtpt[-999] is out of bounds
//
//   std::vector<int> match = JetMapper::Match(nref, jteta, jtphi,
//                                             ngen, geneta, genphi, 0.4);
//
// Match modes, dR < R/2 for cone radius R by default (JME convention):
// OneToOne (default): one to one, closest pairs first, so no jet is matched
// twice. Nearest: each reco jet takes its closest gen jet, which two reco
// jets can share. dRFraction sets the limit as a fraction of R; the forest's
// ref* matching is OneToOne with dR < R:
//
//   JetMapper::Match(nref, jteta, jtphi, ngen, geneta, genphi, 0.4,
//                    JetMapper::Mode::OneToOne, 1.0);
//
// JME: JER twiki / CMSSW SmearedJetProducerT matching -- each reco jet takes
// the closest gen jet with dR < R/2 and |pT - pT_gen| < 3 sigma_JER pT
// (JEC-corrected pT, sigma_JER per reco jet, e.g. JetSmearer::Resolution):
//
//   std::vector<double> sigma(nref);
//   for (int i = 0; i < nref; i++) {
//     sigma[i] = smearer.Resolution(ptCorr[i], jteta[i], rho);
//   }
//   JetMapper::Match(nref, jteta, jtphi, ngen, geneta, genphi, 0.4,
//                    JetMapper::Mode::JME, 0.5, ptCorr, genpt, sigma.data());
//
// Maps compose: match[order[0]] is the leading jet's gen jet

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace JetMapper {

// unmatched entry of a map, same filler as the forest's refpt
constexpr int kUnmatched = -999;

// equal pT keep their original order
template <typename T> std::vector<int> Order(int n, const T *pt) {
  std::vector<int> order(n);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(),
                   [&](int a, int b) { return pt[a] > pt[b]; });
  return order;
}

enum class Mode { Nearest, OneToOne, JME };

inline double DeltaR(double eta1, double phi1, double eta2, double phi2) {
  const double dphi = std::remainder(phi1 - phi2, 2 * M_PI);
  const double deta = eta1 - eta2;
  return std::sqrt(deta * deta + dphi * dphi);
}

template <typename T>
std::vector<int> Match(int nref, const T *jteta, const T *jtphi, int ngen,
                       const T *geneta, const T *genphi, double coneR,
                       Mode mode = Mode::OneToOne, double dRFraction = 0.5,
                       const T *jtpt = nullptr, const T *genpt = nullptr,
                       const double *sigmaJER = nullptr, double nSigma = 3.0) {
  const double maxDR = dRFraction * coneR;
  std::vector<int> match(nref, kUnmatched);
  // arrays only needed with jets to compare (empty vectors give nullptr)
  if (mode == Mode::JME && nref > 0 && ngen > 0 &&
      (!jtpt || !genpt || !sigmaJER)) {
    throw std::invalid_argument(
        "JetMapper: Mode::JME needs jtpt, genpt and sigmaJER");
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

  // one to one: every pair within maxDR, smallest dR first
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

} // namespace JetMapper

#endif
