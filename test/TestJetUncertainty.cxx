// TestJetUncertainty
// hand-computable uncertainty file written at runtime, plus every txt/
// uncertainty file loads and gives a sane central-jet value

#include "JetUncertainty.h"

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

static void Near(std::pair<double, double> got, double down, double up,
                 const std::string &what) {
  bool ok = std::abs(got.first - down) < 1e-9 && std::abs(got.second - up) < 1e-9;
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

// rows: eta range, count, then (pT, down, up) triplets
static void TestInterpolation() {
  std::printf("pT interpolation and edges\n");
  std::string file = Write("Unc.txt", "{1 JetEta 1 JetPt \"\" Correction Uncertainty}\n"
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
  TestInterpolation();
  TestRealFiles();
  std::printf("%d/%d checks passed\n", nCheck - nFail, nCheck);
  return nFail == 0 ? 0 : 1;
}
