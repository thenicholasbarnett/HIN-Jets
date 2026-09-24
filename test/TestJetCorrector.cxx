// TestJetCorrector
// hand-computable correction files written at runtime, plus every txt/ file
// loads and gives a sane central-jet correction

#include "JetCorrector.h"

#include <cmath>
#include <cstdio>
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
    // uncertainty files are JetUncertainty.h's job, JER files JetSmearer.h's
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

int main() {
  if (!std::getenv("TEST_TMPDIR")) {
    std::printf("TEST_TMPDIR not set, run through test/run_tests.sh\n");
    return 1;
  }
  TestL2Relative();
  TestChain();
  TestRealFiles();
  std::printf("%d/%d checks passed\n", nCheck - nFail, nCheck);
  return nFail == 0 ? 0 : 1;
}
