// TestJetMapper
// pT order maps on hand-made arrays, ties, reordering, and random arrays;
// reco to gen maps on hand-made jets: dR limit, phi wrap, nearest vs one to
// one, JME mode, unmatched filler, and composing order + match maps

#include "JetMapper.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

static int nFail = 0;
static int nCheck = 0;

static void Check(bool ok, const std::string &what) {
  nCheck++;
  if (!ok) {
    nFail++;
    std::printf("  FAIL: %s\n", what.c_str());
  }
}

using JetMapper::Match;
using JetMapper::Mode;
static const int kU = JetMapper::kUnmatched;

static void TestOrder() {
  std::printf("order map\n");
  const float pt[4] = {50, 80, 40, 30};
  const float eta[4] = {0.1f, -1.2f, 2.0f, 0.5f};
  std::vector<int> order = JetMapper::Order(4, pt);
  Check(order == std::vector<int>({1, 0, 2, 3}), "{50,80,40,30} -> {1,0,2,3}");
  Check(eta[order[0]] == -1.2f, "leading eta read through the map");

  std::vector<float> etaSorted;
  for (int j : order) {
    etaSorted.push_back(eta[j]);
  }
  Check(etaSorted == std::vector<float>({-1.2f, 0.1f, 2.0f, 0.5f}),
        "eta read in pT order through the map");

  const double ties[4] = {10, 20, 20, 5};
  Check(JetMapper::Order(4, ties) == std::vector<int>({1, 2, 0, 3}),
        "equal pT keep their original order (double input)");
  Check(JetMapper::Order(0, pt).empty(), "no jets, empty map");
}

static void TestOrderRandom() {
  std::printf("order, random arrays\n");
  std::mt19937 rng(7);
  std::uniform_real_distribution<float> u(0, 500);
  bool ok = true;
  for (int trial = 0; trial < 1000; trial++) {
    const int n = trial % 40;
    std::vector<float> pt(n);
    for (auto &p : pt) {
      p = u(rng);
    }
    std::vector<int> order = JetMapper::Order(n, pt.data());
    std::vector<bool> seen(n, false);
    for (int i = 0; i < n; i++) {
      ok = ok && order[i] >= 0 && order[i] < n && !seen[order[i]];
      seen[order[i]] = true;
      ok = ok && (i == 0 || pt[order[i - 1]] >= pt[order[i]]);
    }
  }
  Check(ok, "1000 random arrays: permutation, pT descending");
}

static void TestBasic() {
  std::printf("basic matching\n");
  const float jteta[3] = {0.0f, 1.0f, -2.0f};
  const float jtphi[3] = {0.0f, 1.0f, 2.0f};
  const float geneta[2] = {1.05f, 0.1f};
  const float genphi[2] = {1.0f, 0.0f};
  std::vector<int> m = Match(3, jteta, jtphi, 2, geneta, genphi, 0.4);
  Check(m == std::vector<int>({1, 0, kU}),
        "{1, 0, unmatched}: R = 0.4, within R/2");
  Check(std::fabs(JetMapper::DeltaR(0, 0, 0.1, 0) - 0.1) < 1e-9, "DeltaR");
  Check(Match(3, jteta, jtphi, 2, geneta, genphi, 0.08) ==
            std::vector<int>({kU, kU, kU}),
        "R = 0.08: nothing within 0.04");
  Check(Match(1, jteta, jtphi, 0, geneta, genphi, 0.4) ==
            std::vector<int>({kU}),
        "no gen jets, all unmatched");
}

static void TestPhiWrap() {
  std::printf("phi wrap\n");
  const double jteta[1] = {0.5};
  const double jtphi[1] = {3.1};
  const double geneta[1] = {0.5};
  const double genphi[1] = {-3.1};
  Check(std::fabs(JetMapper::DeltaR(0.5, 3.1, 0.5, -3.1) - (2 * M_PI - 6.2)) <
            1e-9,
        "dphi across +-pi");
  Check(Match(1, jteta, jtphi, 1, geneta, genphi, 0.4) == std::vector<int>({0}),
        "matched across +-pi");
}

static void TestModes() {
  std::printf("nearest vs one to one\n");
  // reco 0 at dR 0.05 from gen 0; reco 1 at dR 0.03 from gen 0 and 0.15
  // from gen 1
  const float jteta[2] = {0.05f, -0.03f};
  const float jtphi[2] = {0.0f, 0.0f};
  const float geneta[2] = {0.0f, -0.18f};
  const float genphi[2] = {0.0f, 0.0f};
  Check(Match(2, jteta, jtphi, 2, geneta, genphi, 0.4, Mode::Nearest) ==
            std::vector<int>({0, 0}),
        "nearest: both reco jets take gen 0");
  Check(Match(2, jteta, jtphi, 2, geneta, genphi, 0.4) ==
            std::vector<int>({kU, 0}),
        "one to one (default): gen 0 to the closer reco 1, reco 0 has nothing "
        "else "
        "in 0.2");
  Check(Match(2, jteta, jtphi, 2, geneta, genphi, 0.4, Mode::OneToOne, 0.75) ==
            std::vector<int>({1, 0}),
        "one to one, dR < 0.75 R = 0.3: reco 0 falls back to gen 1");
}

static void TestCompose() {
  std::printf("order + match maps\n");
  const float jtpt[3] = {40, 90, 60};
  const float jteta[3] = {0.0f, 1.0f, -1.0f};
  const float jtphi[3] = {0.0f, 1.0f, -1.0f};
  const float genpt[2] = {85, 38};
  const float geneta[2] = {1.02f, 0.01f};
  const float genphi[2] = {1.01f, 0.02f};
  std::vector<int> order = JetMapper::Order(3, jtpt);
  std::vector<int> match = Match(3, jteta, jtphi, 2, geneta, genphi, 0.4);
  Check(genpt[match[order[0]]] == 85, "leading reco jet's gen pT");
  Check(match[order[1]] == kU, "subleading reco jet unmatched");
  Check(genpt[match[order[2]]] == 38, "third reco jet's gen pT");
}

// JME: dR < R/2 and |pT - pT_gen| < 3 sigma pT, closest passing pair
static void TestJME() {
  std::printf("Mode::JME\n");
  // reco jet: pT 100, sigma 0.1 -> window |dpT| < 30
  const float jtpt[1] = {100};
  const float jteta[1] = {0.0f};
  const float jtphi[1] = {0.0f};
  const double sigma[1] = {0.1};
  // gen 0: dR 0.05 but pT 60 (fails the window); gen 1: dR 0.15, pT 95
  const float genpt[2] = {60, 95};
  const float geneta[2] = {0.05f, 0.15f};
  const float genphi[2] = {0.0f, 0.0f};
  Check(Match(1, jteta, jtphi, 2, geneta, genphi, 0.4) == std::vector<int>({0}),
        "OneToOne ignores pT: closest gen 0");
  Check(Match(1, jteta, jtphi, 2, geneta, genphi, 0.4, Mode::JME, 0.5, jtpt,
              genpt, sigma) == std::vector<int>({1}),
        "JME skips gen 0 (outside 3 sigma), takes gen 1");
  const float lowpt[2] = {60, 50};
  Check(Match(1, jteta, jtphi, 2, geneta, genphi, 0.4, Mode::JME, 0.5, jtpt,
              lowpt, sigma) == std::vector<int>({kU}),
        "JME: nothing inside the window, unmatched");
  bool threw = false;
  try {
    Match(1, jteta, jtphi, 2, geneta, genphi, 0.4, Mode::JME);
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  Check(threw, "JME without pT/sigma arrays throws");
  std::vector<double> none;
  Check(Match(0, jteta, jtphi, 2, geneta, genphi, 0.4, Mode::JME, 0.5,
              (const float *)nullptr, genpt, none.data())
            .empty(),
        "JME with no reco jets (empty sigma vector) doesn't throw");
}

// default mode on random events: no gen jet used twice, all within R/2
static void TestOneToOne() {
  std::printf("one to one on random events\n");
  std::mt19937 rng(11);
  std::uniform_real_distribution<float> ueta(-2.5, 2.5), uphi(-M_PI, M_PI),
      ushift(-0.25, 0.25);
  bool unique = true, inCone = true;
  int nShared = 0;
  for (int trial = 0; trial < 1000; trial++) {
    const int ngen = 1 + trial % 15;
    std::vector<float> geta(ngen), gphi(ngen);
    for (int g = 0; g < ngen; g++) {
      geta[g] = ueta(rng);
      gphi[g] = uphi(rng);
    }
    // reco jets scattered around the gen jets, some two per gen jet
    std::vector<float> reta, rphi;
    for (int g = 0; g < ngen; g++) {
      for (int k = 0; k < 1 + g % 2; k++) {
        reta.push_back(geta[g] + ushift(rng));
        rphi.push_back(gphi[g] + ushift(rng));
      }
    }
    const int nref = (int)reta.size();
    std::vector<int> m = Match(nref, reta.data(), rphi.data(), ngen,
                               geta.data(), gphi.data(), 0.4);
    std::vector<int> nearest =
        Match(nref, reta.data(), rphi.data(), ngen, geta.data(), gphi.data(),
              0.4, Mode::Nearest);
    std::vector<int> uses(ngen, 0), usesNearest(ngen, 0);
    for (int i = 0; i < nref; i++) {
      if (m[i] >= 0) {
        unique = unique && ++uses[m[i]] == 1;
        inCone = inCone && JetMapper::DeltaR(reta[i], rphi[i], geta[m[i]],
                                             gphi[m[i]]) < 0.2;
      }
      if (nearest[i] >= 0 && ++usesNearest[nearest[i]] == 2) {
        nShared++;
      }
    }
  }
  Check(unique, "OneToOne: no gen jet matched twice");
  Check(inCone, "OneToOne: every match within R/2");
  Check(nShared > 0, "the events do have shared gen jets under Nearest");
}

int main() {
  TestOrder();
  TestOrderRandom();
  TestBasic();
  TestPhiWrap();
  TestModes();
  TestCompose();
  TestOneToOne();
  TestJME();
  std::printf("%d/%d checks passed\n", nCheck - nFail, nCheck);
  return nFail == 0 ? 0 : 1;
}
