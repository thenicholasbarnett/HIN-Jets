// TestJetSmearer
// hand-computable resolution + scale factor files written at runtime,
// hybrid smearing checked against the SmearedJetProducerT formulas by hand,
// plus every real txt/ PtResolution + SF pair loads

#include "JetSmearer.h"

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

static void Near(double got, double want, double tol, const std::string &what) {
  Check(std::abs(got - want) <= tol, what + ": got " + std::to_string(got) +
                                         ", want " + std::to_string(want));
}

static std::string Write(const std::string &name, const std::string &text) {
  std::string path = std::string(std::getenv("TEST_TMPDIR")) + "/" + name;
  std::ofstream f(path);
  f << text;
  return path;
}

// sigma_JER 0.1 for eta < 0, 0.2 for eta >= 0
static std::string ResolutionFile() {
  return Write("PtResolution.txt", "{1 JetEta 1 JetPt [0] Resolution}\n"
                                   "-5.191 0 3 0 10000 0.1\n"
                                   "0 5.191 3 0 10000 0.2\n");
}

// official SF format, nominal down up
static std::string ScaleFactorFile() {
  return Write("SF.txt", "{1 JetEta 0 None ScaleFactor}\n"
                         "-5.191 0 3 0.95 0.9 1.0\n"
                         "0 5.191 3 1.2 1.1 1.3\n");
}

static void TestMatched() {
  std::printf("gen-matched scaling branch\n");
  JetSmearer s(ResolutionFile(), ScaleFactorFile());
  JetSmearing::Result r = s.Smear(100.0, 0.5, 1.5, 98.0);
  Check(r.matched, "close gen jet is matched");
  Near(r.resolution, 0.2, 1e-6, "resolution read from file");
  Near(r.scaleFactor, 1.2, 1e-6, "nominal scale factor read from file");
  Near(r.smearFactor, 1.0 + 0.2 * 2.0 / 100.0, 1e-6, "scaling formula");
  Near(s.SmearedPt(100.0, 0.5, 1.5, 98.0), 100.4, 1e-4, "SmearedPt");

  Near(s.Smear(100.0, 0.5, 1.5, 98.0, Variation::UP).smearFactor,
       1.0 + 0.3 * 2.0 / 100.0, 1e-6, "UP variation");
  Near(s.Smear(100.0, 0.5, 1.5, 98.0, Variation::DOWN).smearFactor,
       1.0 + 0.1 * 2.0 / 100.0, 1e-6, "DOWN variation");
  Near(s.Smear(100.0, -0.5, 1.5, 98.0).smearFactor,
       1.0 + (0.95 - 1.0) * 2.0 / 100.0, 1e-6, "SF < 1 narrows when matched");

  // |reco - gen| < 3 sigma pT is the match window
  Check(s.Smear(100.0, 0.5, 1.5, 41.0).matched, "59 GeV off, inside 3 sigma");
  Check(!s.Smear(100.0, 0.5, 1.5, 39.0).matched, "61 GeV off, outside 3 sigma");
}

static void TestStochastic() {
  std::printf("stochastic branch\n");
  JetSmearer s(ResolutionFile(), ScaleFactorFile());
  const int n = 200000;
  double sum = 0, sum2 = 0;
  for (int i = 0; i < n; i++) {
    JetSmearing::Result r = s.Smear(100.0, 0.5, 1.5, -1.0);
    if (i == 0) {
      Check(!r.matched, "genPt < 0 is unmatched");
    }
    sum += r.smearFactor;
    sum2 += r.smearFactor * r.smearFactor;
  }
  double mean = sum / n;
  double sigma = std::sqrt(sum2 / n - mean * mean);
  double want = 0.2 * std::sqrt(1.2 * 1.2 - 1.0);
  Near(mean, 1.0, 0.002, "stochastic mean");
  Near(sigma, want, 0.01 * want, "stochastic width sigma*sqrt(SF^2-1)");

  Near(s.Smear(100.0, -0.5, 1.5, -1.0).smearFactor, 1.0, 1e-12,
       "SF < 1 unmatched is left alone");

  JetSmearer a(ResolutionFile(), ScaleFactorFile(), 7);
  JetSmearer b(ResolutionFile(), ScaleFactorFile(), 7);
  JetSmearer c(ResolutionFile(), ScaleFactorFile(), 8);
  bool same = true, diff = false;
  for (int i = 0; i < 10; i++) {
    double fa = a.Smear(100.0, 0.5, 1.5, -1.0).smearFactor;
    double fb = b.Smear(100.0, 0.5, 1.5, -1.0).smearFactor;
    double fc = c.Smear(100.0, 0.5, 1.5, -1.0).smearFactor;
    same = same && fa == fb;
    diff = diff || fa != fc;
  }
  Check(same, "same seed, same sequence");
  Check(diff, "different seed, different sequence");
}

// Scaling / Stochastic only, vs the Hybrid default
static void TestMethods() {
  std::printf("smearing methods\n");
  using JetSmearing::Method;
  JetSmearer hybrid(ResolutionFile(), ScaleFactorFile());
  Check(hybrid.GetMethod() == Method::Hybrid, "default method is hybrid");

  JetSmearer scaling(ResolutionFile(), ScaleFactorFile(),
                     JetSmearer::kDefaultSeed, Method::Scaling);
  Near(scaling.Smear(100.0, 0.5, 1.5, 98.0).smearFactor, 1.004, 1e-6,
       "scaling: matched jet scaled");
  Near(scaling.Smear(100.0, 0.5, 1.5, -1.0).smearFactor, 1.0, 1e-12,
       "scaling: unmatched jet left alone");
  Near(scaling.Smear(100.0, 0.5, 1.5, 39.0).smearFactor, 1.0, 1e-12,
       "scaling: outside the match window left alone");

  JetSmearer stochastic(ResolutionFile(), ScaleFactorFile(),
                        JetSmearer::kDefaultSeed, Method::Stochastic);
  const int n = 200000;
  double sum = 0, sum2 = 0;
  bool anyMatched = false;
  for (int i = 0; i < n; i++) {
    JetSmearing::Result r = stochastic.Smear(100.0, 0.5, 1.5, 98.0);
    anyMatched = anyMatched || r.matched;
    sum += r.smearFactor;
    sum2 += r.smearFactor * r.smearFactor;
  }
  const double mean = sum / n;
  const double want = 0.2 * std::sqrt(1.2 * 1.2 - 1.0);
  Check(!anyMatched, "stochastic: gen match ignored");
  Near(mean, 1.0, 0.002, "stochastic: mean");
  Near(std::sqrt(sum2 / n - mean * mean), want, 0.01 * want,
       "stochastic: width even for a matched jet");

  stochastic.SetMethod(Method::Hybrid);
  Check(stochastic.Smear(100.0, 0.5, 1.5, 98.0).matched, "SetMethod");

  Check(JetSmearing::MethodFromString("stochastic") == Method::Stochastic,
        "MethodFromString");
  bool threw = false;
  try {
    JetSmearing::MethodFromString("gaussian");
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  Check(threw, "MethodFromString rejects unknown names");
}

static void TestEdges() {
  std::printf("edges\n");
  JetSmearer s(ResolutionFile(), ScaleFactorFile());
  Near(s.Smear(100.0, 6.0, 1.5, 98.0).smearFactor, 1.0, 1e-12,
       "eta outside the files is not smeared");
  Check(JetSmearing::SmearedPt(10.0, -5.0) > 0.0, "pT floor");
  Near(JetSmearing::SmearedPt(10.0, 1.1), 11.0, 1e-12, "SmearedPt helper");
}

// every txt/ PtResolution file has a matching SF file that loads and gives a
// sane central-jet resolution and scale factor
static void TestRealFiles() {
  std::printf("every txt/ JER pair\n");
  std::vector<std::string> dirs;
  if (DIR *top = opendir("txt")) {
    while (dirent *e = readdir(top)) {
      if (e->d_type == DT_DIR && e->d_name[0] != '.') {
        dirs.push_back(std::string("txt/") + e->d_name);
      }
    }
    closedir(top);
  }
  int nPairs = 0;
  for (const std::string &dir : dirs) {
    DIR *d = opendir(dir.c_str());
    if (!d) {
      continue;
    }
    while (dirent *e = readdir(d)) {
      std::string name = e->d_name;
      size_t at = name.find("_PtResolution_");
      if (at == std::string::npos) {
        continue;
      }
      std::string res = dir + "/" + name;
      std::string sf = dir + "/" + name.substr(0, at) + "_SF_" +
                       name.substr(at + 14);
      JetSmearer s(res, sf);
      JetSmearing::Result r = s.Smear(100.0, 0.5, 5.0, 98.0);
      Check(r.resolution > 0.03 && r.resolution < 0.3,
            res + " sigma " + std::to_string(r.resolution));
      Check(r.scaleFactor > 0.9 && r.scaleFactor < 1.5,
            sf + " SF " + std::to_string(r.scaleFactor));
      nPairs++;
    }
    closedir(d);
  }
  Check(nPairs > 0, "found JER pairs in txt/");
}

int main() {
  if (!std::getenv("TEST_TMPDIR")) {
    std::printf("TEST_TMPDIR not set, run through test/run_tests.sh\n");
    return 1;
  }
  TestMatched();
  TestStochastic();
  TestMethods();
  TestEdges();
  TestRealFiles();
  std::printf("%d/%d checks passed\n", nCheck - nFail, nCheck);
  return nFail == 0 ? 0 : 1;
}
