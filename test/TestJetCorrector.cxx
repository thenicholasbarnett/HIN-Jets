// TestJetCorrector
// hand-computable correction and uncertainty files written at runtime, the
// JES up/down variations, plus every txt/ correction and uncertainty file
// loads and gives a sane central-jet value

#include "JetCorrector.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
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

static void Near(double got, double want, const std::string &what) {
  bool ok = std::abs(got - want) <= 1e-9 * std::max(1.0, std::abs(want));
  Check(ok, what + ": got " + std::to_string(got) + ", want " +
                std::to_string(want));
}

static void Near(std::pair<double, double> got, double down, double up,
                 const std::string &what) {
  bool ok =
      std::abs(got.first - down) < 1e-9 && std::abs(got.second - up) < 1e-9;
  Check(ok, what + ": got (" + std::to_string(got.first) + ", " +
                std::to_string(got.second) + "), want (" +
                std::to_string(down) + ", " + std::to_string(up) + ")");
}

static std::string Write(const std::string &name, const std::string &text) {
  std::string path = std::string(std::getenv("TEST_TMPDIR")) + "/" + name;
  std::ofstream f(path);
  f << text;
  return path;
}

static void TestL2Relative() {
  std::printf("L2Relative-style, 1 bin var, 1 dependency\n");
  std::string file = Write("L2.txt", "{1 JetEta 1 JetPt [0]+[1]*log10(x) "
                                     "Correction L2Relative}\n"
                                     "-5.191 0 4 10 1000 1.1 0.01\n"
                                     "0 5.191 4 10 1000 1.2 -0.02\n");
  JetCorrector jec(file);
  jec.SetJetEta(-1.0);
  jec.SetJetPT(100.0);
  Near(jec.GetCorrectedPT(), 100.0 * (1.1 + 0.01 * 2.0), "negative eta bin");
  jec.SetJetEta(1.0);
  Near(jec.GetCorrectedPT(), 100.0 * (1.2 - 0.02 * 2.0), "positive eta bin");
  jec.SetJetPT(5.0);
  Near(jec.GetCorrectedPT(), 5.0 * (1.2 - 0.02 * 1.0), "pT clamped to low");
  jec.SetJetPT(2000.0);
  Near(jec.GetCorrectedPT(), 2000.0 * (1.2 - 0.02 * 3.0), "pT clamped to high");
  jec.SetJetPT(100.0);
  jec.SetJetEta(0.0);
  Near(jec.GetCorrectedPT(), 100.0 * (1.1 + 0.01 * 2.0),
       "shared edge goes to first listed bin");
  jec.SetJetEta(6.0);
  Check(jec.GetCorrectedPT() < 0, "eta outside every bin gives -1");
  Near(jec.GetCorrection(), -1.0, "GetCorrection -1 outside bins");
}

static void TestChain() {
  std::printf("L1FastJet (3 dependencies) -> L2Relative chain\n");
  std::string l1 = Write("L1.txt", "{1 JetEta 3 Rho JetPt JetA "
                                   "max(0.0001,1-(z/y)*([1]*(x-[0]))) "
                                   "Correction L1FastJet}\n"
                                   "-5.191 5.191 8 0 70 1 6500 0 10 1.5 0.5\n");
  std::string l2 = Write("L2chain.txt", "{1 JetEta 1 JetPt [0]+[1]*log10(x) "
                                        "Correction L2Relative}\n"
                                        "-5.191 5.191 4 10 1000 1.1 0.01\n");
  JetCorrector jec(std::vector<std::string>{l1, l2});
  jec.SetJetEta(0.5);
  jec.SetJetPhi(0.0);
  jec.SetJetPT(50.0);
  jec.SetRho(3.0);
  jec.SetJetArea(0.5);
  double afterL1 = 50.0 * (1.0 - (0.5 / 50.0) * (0.5 * (3.0 - 1.5)));
  double afterL2 = afterL1 * (1.1 + 0.01 * std::log10(afterL1));
  Near(jec.GetCorrectedPT(), afterL2, "L1 then L2");
  Near(jec.GetCorrection(), afterL2 / 50.0, "total correction factor");
}

// rows: eta range, count, then (pT, down, up) triplets
static void TestUncertainty() {
  std::printf("uncertainty: pT interpolation and edges\n");
  std::string file =
      Write("Unc.txt", "{1 JetEta 1 JetPt \"\" Correction Uncertainty}\n"
                       "-5 0 9 10 0.05 0.06 100 0.03 0.04 1000 0.01 0.02\n"
                       "0 5 6 10 0.10 0.10 100 0.02 0.02\n");
  JetUncertainty jeu(file);
  auto at = [&](double pt, double eta) {
    jeu.SetJetPT(pt);
    jeu.SetJetEta(eta);
    jeu.SetJetPhi(0.0);
    return jeu.GetUncertainty();
  };
  Near(at(10.0, -1.0), 0.05, 0.06, "first pT point");
  Near(at(55.0, -1.0), 0.04, 0.05, "linear between 10 and 100");
  Near(at(100.0, -1.0), 0.03, 0.04, "exact middle point");
  Near(at(550.0, -1.0), 0.02, 0.03, "linear between 100 and 1000");
  Near(at(5.0, -1.0), 0.05, 0.06, "below first point clamps");
  Near(at(2000.0, -1.0), 0.01, 0.02, "above last point clamps");
  Near(at(55.0, 1.0), 0.06, 0.06, "second eta bin");
  Near(at(55.0, 0.0), 0.04, 0.05, "shared edge goes to first listed bin");
  Near(at(55.0, 6.0), -1.0, -1.0, "eta outside every bin gives -1");
}

// Variation: corrected pT * (1 + up) / (1 - down), uncertainty at the
// corrected pT
static void TestVariations() {
  std::printf("JES variations\n");
  std::string l2 = Write("L2var.txt", "{1 JetEta 1 JetPt [0] "
                                      "Correction L2Relative}\n"
                                      "-5.191 5.191 4 1 5000 1.1\n");
  std::string unc =
      Write("UncVar.txt", "{1 JetEta 1 JetPt \"\" Correction Uncertainty}\n"
                          "-3 3 6 10 0.05 0.06 1000 0.05 0.06\n");
  JetCorrector jec(std::vector<std::string>{l2}, unc);
  jec.SetJetEta(0.5);
  jec.SetJetPhi(0.0);
  jec.SetJetPT(100.0);
  Near(jec.GetCorrectedPT(), 110.0, "nominal, unchanged method");
  Near(jec.GetCorrectedPT(Variation::NOMINAL), 110.0, "Variation::NOMINAL");
  Near(jec.GetCorrectedPT(Variation::UP), 110.0 * 1.06, "UP = pT (1 + up)");
  Near(jec.GetCorrectedPT(Variation::DOWN), 110.0 * 0.95,
       "DOWN = pT (1 - down)");
  std::pair<double, double> u = jec.GetUncertainty();
  Check(std::abs(u.first - 0.05) < 1e-9 && std::abs(u.second - 0.06) < 1e-9,
        "GetUncertainty {down, up} = {0.05, 0.06}");
  const double sNom = 104.5; // e.g. a smeared pT
  Near(sNom * (1 + u.second), 104.5 * 1.06, "JES up on a smeared pT");
  jec.SetJetEta(4.0);
  Near(jec.GetCorrectedPT(), 110.0, "nominal outside the uncertainty eta");
  Near(jec.GetCorrectedPT(Variation::UP), -1.0,
       "UP -1 outside the uncertainty eta");
  Check(jec.GetUncertainty().second == -1.0,
        "GetUncertainty -1 outside the uncertainty eta");

  JetCorrector noUnc(l2);
  noUnc.SetJetEta(0.5);
  noUnc.SetJetPT(100.0);
  bool threw = false;
  try {
    noUnc.GetCorrectedPT(Variation::UP);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  Check(threw, "UP without an uncertainty file throws");
  threw = false;
  try {
    noUnc.GetUncertainty();
  } catch (const std::runtime_error &) {
    threw = true;
  }
  Check(threw, "GetUncertainty without an uncertainty file throws");
  Near(noUnc.GetCorrectedPT(Variation::NOMINAL), 110.0,
       "NOMINAL without an uncertainty file");
}

static std::vector<std::string> TxtFiles(const std::string &dir) {
  std::vector<std::string> out;
  DIR *d = opendir(dir.c_str());
  if (!d) {
    return out;
  }
  while (dirent *e = readdir(d)) {
    std::string name = e->d_name;
    if (name == "." || name == "..") {
      continue;
    }
    std::string path = dir + "/" + name;
    if (e->d_type == DT_DIR) {
      for (const auto &p : TxtFiles(path)) {
        out.push_back(p);
      }
    } else if (name.size() > 4 && name.substr(name.size() - 4) == ".txt") {
      out.push_back(path);
    }
  }
  closedir(d);
  return out;
}

static void TestRealFiles() {
  std::printf("every txt/ correction file\n");
  std::vector<std::string> files = TxtFiles("txt");
  Check(!files.empty(), "found txt/ files");
  for (const auto &file : files) {
    std::ifstream f(file);
    std::string first;
    std::getline(f, first);
    // uncertainty files are tested below, JER files by TestJetSmearer
    if (first.find("Correction") == std::string::npos ||
        first.find("Uncertainty") != std::string::npos) {
      continue;
    }
    JetCorrector jec(file);
    jec.SetJetEta(0.1);
    jec.SetJetPhi(0.0);
    jec.SetJetPT(50.0);
    jec.SetRho(2.0);
    jec.SetJetArea(0.5);
    double c = jec.GetCorrection();
    Check(std::isfinite(c) && c > 0.5 && c < 2.0,
          file + " central correction " + std::to_string(c));
  }
}

static void TestRealUncertaintyFiles() {
  std::printf("every txt/ uncertainty file\n");
  int n = 0;
  for (const auto &file : TxtFiles("txt")) {
    std::ifstream f(file);
    std::string first;
    std::getline(f, first);
    if (first.find("Uncertainty") == std::string::npos) {
      continue;
    }
    JetUncertainty jeu(file);
    jeu.SetJetPT(100.0);
    jeu.SetJetEta(0.1);
    jeu.SetJetPhi(0.0);
    std::pair<double, double> u = jeu.GetUncertainty();
    Check(std::isfinite(u.first) && std::isfinite(u.second) && u.first > 0 &&
              u.second > 0 && u.first < 0.2 && u.second < 0.2,
          file + " central (" + std::to_string(u.first) + ", " +
              std::to_string(u.second) + ")");
    n++;
  }
  Check(n > 0, "found uncertainty files in txt/");
}

int main() {
  if (!std::getenv("TEST_TMPDIR")) {
    std::printf("TEST_TMPDIR not set, run through test/run_tests.sh\n");
    return 1;
  }
  TestL2Relative();
  TestChain();
  TestUncertainty();
  TestVariations();
  TestRealFiles();
  TestRealUncertaintyFiles();
  std::printf("%d/%d checks passed\n", nCheck - nFail, nCheck);
  return nFail == 0 ? 0 : 1;
}
