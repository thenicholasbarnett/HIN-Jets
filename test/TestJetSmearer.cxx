// TestJetSmearer
// hand-computable resolution + scale factor files written at runtime,
// smearing checked against the hybrid formulas by hand, JERSmear against
// correctionlib's own output on JME's jer_smear.json, plus every real txt/
// PtResolution + SF pair loads

#include "JetSmearer.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
  JetSmearing::Result r = s.Smear(100.0, 0.5, 1.5, 98.0, 1);
  Check(r.matched, "close gen jet is matched");
  Near(r.resolution, 0.2, 1e-6, "resolution read from file");
  Near(r.scaleFactor, 1.2, 1e-6, "nominal scale factor read from file");
  Near(r.smearFactor, 1.0 + 0.2 * 2.0 / 100.0, 1e-6, "scaling formula");
  Near(s.SmearedPt(100.0, 0.5, 1.5, 98.0, 1), 100.4, 1e-4, "SmearedPt");

  Near(s.Smear(100.0, 0.5, 1.5, 98.0, 1, Variation::UP).smearFactor,
       1.0 + 0.3 * 2.0 / 100.0, 1e-6, "UP variation");
  Near(s.Smear(100.0, 0.5, 1.5, 98.0, 1, Variation::DOWN).smearFactor,
       1.0 + 0.1 * 2.0 / 100.0, 1e-6, "DOWN variation");
  Near(s.Smear(100.0, -0.5, 1.5, 98.0, 1).smearFactor,
       1.0 + (0.95 - 1.0) * 2.0 / 100.0, 1e-6, "SF < 1 narrows when matched");

  // |reco - gen| < 3 sigma pT is the match window
  Check(s.Smear(100.0, 0.5, 1.5, 41.0, 1).matched,
        "59 GeV off, inside 3 sigma");
  Check(!s.Smear(100.0, 0.5, 1.5, 39.0, 1).matched,
        "61 GeV off, outside 3 sigma");
}

static void TestStochastic() {
  std::printf("stochastic branch\n");
  JetSmearer s(ResolutionFile(), ScaleFactorFile());
  const int n = 200000;
  double sum = 0, sum2 = 0;
  for (int i = 0; i < n; i++) {
    JetSmearing::Result r = s.Smear(100.0, 0.5, 1.5, -1.0, i);
    if (i == 0) {
      Check(!r.matched, "genPt < 0 is unmatched");
    }
    sum += r.smearFactor;
    sum2 += r.smearFactor * r.smearFactor;
  }
  double mean = sum / n;
  double sigma = std::sqrt(sum2 / n - mean * mean);
  double want = 0.2 * std::sqrt(1.2 * 1.2 - 1.0);
  Near(mean, 1.0, 0.002, "stochastic mean over event IDs");
  Near(sigma, want, 0.01 * want, "stochastic width sigma*sqrt(SF^2-1)");

  Near(s.Smear(100.0, -0.5, 1.5, -1.0, 1).smearFactor, 1.0, 1e-12,
       "SF < 1 unmatched is left alone");

  // the draw depends only on (pT, eta, rho, eventID): no hidden state
  JetSmearer other(ResolutionFile(), ScaleFactorFile());
  const double f1 = s.Smear(100.0, 0.5, 1.5, -1.0, 42).smearFactor;
  for (int i = 0; i < 5; i++) {
    other.Smear(80.0, 0.3, 1.5, -1.0, i);
  }
  Check(other.Smear(100.0, 0.5, 1.5, -1.0, 42).smearFactor == f1,
        "same jet, same event: same smear, whatever came before");
  Check(s.Smear(100.0, 0.5, 1.5, -1.0, 43).smearFactor != f1,
        "different event ID, different smear");
  Check(s.Smear(100.0, 0.5, 1.6, -1.0, 42).smearFactor != f1,
        "different rho, different smear");
}

// JME's correctionlib JERSmear on jer_smear.json (correctionlib 2.9.0):
// {JetPt, JetEta, GenPt, Rho, EventID, JER, JERSF, output}
static void TestJERSmearReference() {
  std::printf("JERSmear vs correctionlib\n");
  struct Row {
    double pt, eta, gen, rho;
    long long eid;
    double jer, sf, out;
  };
  const Row rows[] = {
      {1901.175074170241, 1.2003216515760828, -1.7605236605147883,
       49.996041183624136, 498466514399LL, 0.27191330920905893,
       1.2697298083646846, 0.9158300351062498},
      {292.5984273756693, -2.9957116995516833, -1.3369702368047562,
       5.3796858018933875, 251184097034LL, 0.0813456794597646,
       1.5732070721663405, 0.9781077384817323},
      {1897.5556470388015, 2.2178430356851724, -0.8897957898451898,
       39.72878710804939, 73745862620LL, 0.31865445043870977,
       1.5394817815316748, 1.283631533141122},
      {627.1037467609185, -4.06746032951771, -1.968805762920349,
       24.383448444178285, 787931576372LL, 0.1514453582191491,
       1.1991746744491527, 0.9643934955832941},
      {1026.0841412770121, 0.5590372795786862, -2.4167598158391037,
       42.949664325946955, 67542619370LL, 0.1272446636891928,
       0.9017280067451869, 1.0},
      {219.7850830459586, -4.060994193998628, 150.43631497551985,
       57.052220298905816, 12914606801LL, 0.20807274131368877,
       1.4548616135213392, 1.1435224452136308},
      {329.18086494644473, 3.0548883742955906, 173.54084198906205,
       29.278166134811695, 839826462277LL, 0.17091251445441266,
       1.0084658294555655, 1.0040027292930038},
      {1393.7750427427804, 1.52252903549085, 1310.947289856692,
       29.3045044750378, 786162560040LL, 0.19510196312948058,
       1.1800754725466016, 1.0107013300450438},
      {1086.9413704608742, -1.9940269312924461, 0.0, 22.428722902239823,
       246936153312LL, 0.05889133843847026, 1.0070368556329086,
       1.0070368556329086},
  };
  int exact = 0;
  for (const Row &r : rows) {
    exact += JetSmearing::JERSmear(r.pt, r.eta, r.gen, r.rho, r.eid, r.jer,
                                   r.sf) == r.out;
  }
  Check(exact == (int)(sizeof(rows) / sizeof(rows[0])),
        "bit-identical to correctionlib: " + std::to_string(exact) + "/" +
            std::to_string(sizeof(rows) / sizeof(rows[0])));

  // XXH64 reference vectors (xxHash)
  Check(JetSmearing::detail::XXH64("", 0, 0) == 0xEF46DB3751D8E999ULL,
        "XXH64 of empty input");
  Check(JetSmearing::detail::XXH64("abc", 3, 0) == 0x44BC2CF5AD770999ULL,
        "XXH64 of \"abc\"");
}

// JME / Scaling / Stochastic, vs the Hybrid default
static void TestMethods() {
  std::printf("smearing methods\n");
  using JetSmearing::Method;
  JetSmearer hybrid(ResolutionFile(), ScaleFactorFile());
  Check(hybrid.GetMethod() == Method::Hybrid, "default method is hybrid");

  // JME: genPt as given, no pT window
  JetSmearer jme(ResolutionFile(), ScaleFactorFile(), Method::JME);
  Check(jme.Smear(100.0, 0.5, 1.5, 39.0, 1).matched,
        "JME: 61 GeV off still scaled (matching is the caller's)");
  Near(jme.Smear(100.0, 0.5, 1.5, 39.0, 1).smearFactor,
       1.0 + 0.2 * 61.0 / 100.0, 1e-6, "JME: scaling formula");
  Check(jme.Smear(100.0, 0.5, 1.5, -1.0, 7).smearFactor ==
            hybrid.Smear(100.0, 0.5, 1.5, -1.0, 7).smearFactor,
        "JME and Hybrid share the stochastic draw");

  JetSmearer scaling(ResolutionFile(), ScaleFactorFile(), Method::Scaling);
  Near(scaling.Smear(100.0, 0.5, 1.5, 98.0, 1).smearFactor, 1.004, 1e-6,
       "scaling: matched jet scaled");
  Near(scaling.Smear(100.0, 0.5, 1.5, -1.0, 1).smearFactor, 1.0, 1e-12,
       "scaling: unmatched jet left alone");
  Near(scaling.Smear(100.0, 0.5, 1.5, 39.0, 1).smearFactor, 1.0, 1e-12,
       "scaling: outside the match window left alone");

  JetSmearer stochastic(ResolutionFile(), ScaleFactorFile(),
                        Method::Stochastic);
  const int n = 200000;
  double sum = 0, sum2 = 0;
  bool anyMatched = false;
  for (int i = 0; i < n; i++) {
    JetSmearing::Result r = stochastic.Smear(100.0, 0.5, 1.5, 98.0, i);
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
  Check(stochastic.Smear(100.0, 0.5, 1.5, 98.0, 1).matched, "SetMethod");

  Check(JetSmearing::MethodFromString("jme") == Method::JME,
        "MethodFromString(\"jme\")");
  Check(JetSmearing::MethodFromString("stochastic") == Method::Stochastic,
        "MethodFromString(\"stochastic\")");
  bool threw = false;
  try {
    JetSmearing::MethodFromString("gaussian");
  } catch (const std::invalid_argument &) {
    threw = true;
  }
  Check(threw, "MethodFromString rejects unknown names");
}

static void TestResolution() {
  std::printf("Resolution accessor\n");
  JetSmearer s(ResolutionFile(), ScaleFactorFile());
  Near(s.Resolution(100.0, 0.5, 1.5), 0.2, 1e-6, "sigma_JER eta >= 0");
  Near(s.Resolution(100.0, -0.5, 1.5), 0.1, 1e-6, "sigma_JER eta < 0");
}

static void TestEdges() {
  std::printf("edges\n");
  JetSmearer s(ResolutionFile(), ScaleFactorFile());
  Near(s.Smear(100.0, 6.0, 1.5, 98.0, 1).smearFactor, 1.0, 1e-12,
       "eta outside the files is not smeared");
  Near(JetSmearing::SmearedPt(10.0, -5.0), -50.0, 1e-12,
       "no pT floor (as JERSmear)");
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
      std::string sf =
          dir + "/" + name.substr(0, at) + "_SF_" + name.substr(at + 14);
      JetSmearer s(res, sf);
      JetSmearing::Result r = s.Smear(100.0, 0.5, 5.0, 98.0, 1);
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
  TestJERSmearReference();
  TestMethods();
  TestResolution();
  TestEdges();
  TestRealFiles();
  std::printf("%d/%d checks passed\n", nCheck - nFail, nCheck);
  return nFail == 0 ? 0 : 1;
}
