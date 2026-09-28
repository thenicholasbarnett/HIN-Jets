// TestJetSelector
// jet ID boundaries from the Run 3 AK4 CHS table (cms-jme-jmar docs),
// Run 2 pre-UL ID JSONs vs the twiki selection code on random jets
// veto spot checks taken from the JECDatabase ROOT maps of the same tags

#include "JetSelector.h"

#include <cmath>
#include <cstdio>
#include <random>
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

static const std::string kJetID = "json/jetid.json";
static std::string VetoFile(const std::string &tag) {
  return "json/jetvetomaps_" + tag + ".json";
}

// passes TightLeptonVeto in every eta region
struct Jet {
  double CHF = 0.5, NHF = 0.1, CEF = 0.1, NEF = 0.1, MUF = 0.0;
  int CHM = 5, NHM = 3, CEM = 0, NEM = 10, MUM = 0;
};

static bool ID(const JetSelector &js, double eta, const Jet &j) {
  return js.PassesID(eta, j.CHF, j.NHF, j.CEF, j.NEF, j.MUF, j.CHM, j.NHM,
                     j.CEM, j.NEM, j.MUM);
}

struct Cut {
  const char *name;
  double eta;
  Jet jet;
  bool ppPass;
  bool ionPass;
};

static std::vector<Cut> IDCases() {
  std::vector<Cut> v;
  auto add = [&](const char *name, double eta, bool pp, bool ion,
                 auto modify) {
    Jet j;
    modify(j);
    v.push_back({name, eta, j, pp, ion});
  };

  for (double eta : {0.0, 1.5, -2.0, 2.599, 4.0, -4.9}) {
    add("baseline", eta, true, true, [](Jet &) {});
  }

  // 0 <= |eta| < 2.6
  add("central NHF 0.989", 1.0, true, true, [](Jet &j) { j.NHF = 0.989; });
  add("central NHF 0.99", 1.0, false, true, [](Jet &j) { j.NHF = 0.99; });
  add("central NEF 0.899", 1.0, true, true, [](Jet &j) { j.NEF = 0.899; });
  add("central NEF 0.9", 1.0, false, false, [](Jet &j) { j.NEF = 0.9; });
  add("central CHF 0.01", -1.0, true, true, [](Jet &j) { j.CHF = 0.01; });
  add("central CHF 0.009", -1.0, false, true, [](Jet &j) { j.CHF = 0.009; });
  add("central chMult 0", 1.0, false, true, [](Jet &j) {
    j.CHM = 0;
    j.CEM = 0;
    j.MUM = 0;
  });
  add("central mult 1", 1.0, false, true, [](Jet &j) {
    j.CHM = 1;
    j.NHM = 0;
    j.NEM = 0;
  });
  add("central chMult from electron", 1.0, true, true, [](Jet &j) {
    j.CHM = 0;
    j.CEM = 1;
  });
  add("central MUF 0.8", 1.0, false, false, [](Jet &j) { j.MUF = 0.8; });
  add("central MUF 0.799", 1.0, true, true, [](Jet &j) { j.MUF = 0.799; });
  add("central CEF 0.8", 1.0, false, false, [](Jet &j) { j.CEF = 0.8; });

  // 2.6 <= |eta| < 2.7, 2.6 itself is in here
  add("eta 2.599 NHF 0.95", 2.599, true, true, [](Jet &j) { j.NHF = 0.95; });
  add("eta 2.6 NHF 0.95", 2.6, false, true, [](Jet &j) { j.NHF = 0.95; });
  add("2.6-2.7 NHF 0.9", -2.65, false, true, [](Jet &j) { j.NHF = 0.9; });
  add("2.6-2.7 NEF 0.95", 2.65, true, true, [](Jet &j) { j.NEF = 0.95; });
  add("2.6-2.7 NEF 0.99", 2.65, false, false, [](Jet &j) { j.NEF = 0.99; });
  add("2.6-2.7 CHF 0", 2.65, true, true, [](Jet &j) { j.CHF = 0.0; });
  add("2.6-2.7 chMult 0", 2.65, false, true, [](Jet &j) { j.CHM = 0; });
  add("2.6-2.7 MUF 0.8", 2.65, false, false, [](Jet &j) { j.MUF = 0.8; });
  add("2.6-2.7 CEF 0.8", 2.65, false, false, [](Jet &j) { j.CEF = 0.8; });

  // 2.7 <= |eta| < 3.0
  add("2.7-3.0 NHF 0.99", 2.85, false, true, [](Jet &j) { j.NHF = 0.99; });
  add("2.7-3.0 NEF 0.99", -2.85, false, false, [](Jet &j) { j.NEF = 0.99; });
  add("2.7-3.0 neMult 1", 2.85, false, true, [](Jet &j) {
    j.NHM = 0;
    j.NEM = 1;
  });
  add("2.7-3.0 no charged", 2.85, true, true, [](Jet &j) {
    j.CHF = 0.0;
    j.CHM = 0;
  });
  add("2.7-3.0 MUF 0.9", 2.85, true, true, [](Jet &j) { j.MUF = 0.9; });
  add("2.7-3.0 CEF 0.9", 2.85, true, true, [](Jet &j) { j.CEF = 0.9; });

  // 3.0 <= |eta|
  add("forward NEF 0.399", 4.0, true, true, [](Jet &j) { j.NEF = 0.399; });
  add("forward NEF 0.4", 4.0, false, false, [](Jet &j) { j.NEF = 0.4; });
  add("forward neMult 10", -4.0, false, true, [](Jet &j) {
    j.NHM = 0;
    j.NEM = 10;
  });
  add("forward neMult 11", -4.0, true, true, [](Jet &j) {
    j.NHM = 1;
    j.NEM = 10;
  });
  add("forward NHF 0.999", 4.0, true, true, [](Jet &j) { j.NHF = 0.999; });
  add("forward MUF 0.9", 4.0, true, true, [](Jet &j) { j.MUF = 0.9; });
  add("eta 3.0 NEF 0.5", 3.0, false, false, [](Jet &j) { j.NEF = 0.5; });
  add("eta 2.999 NEF 0.5", 2.999, true, true, [](Jet &j) { j.NEF = 0.5; });
  return v;
}

struct VetoPoint {
  const char *tag;
  double eta, phi;
  bool analysis, calibration;
};

// bin centers, from the JECDatabase ROOT files, HIN-used maps only
static const std::vector<VetoPoint> kVetoPoints = {
    {"Summer23BPixPrompt23_RunD_V1", -2.4110, -1.3526, false, false},
    {"Summer23BPixPrompt23_RunD_V1", -2.4110, 1.8762, false, false},
    {"Summer23BPixPrompt23_RunD_V1", -2.2470, -1.3526, true, true},
    {"Summer23BPixPrompt23_RunD_V1", -2.1075, -2.3998, true, true},
    {"Summer23BPixPrompt23_RunD_V1", -1.4355, -1.1781, false, true},
    {"Summer24Prompt24_RunBCDEFGHI_V1", -3.9260, 2.5744, true, true},
    {"Summer24Prompt24_RunBCDEFGHI_V1", -3.4015, 2.9234, true, true},
    {"Summer24Prompt24_RunBCDEFGHI_V1", -2.4110, 1.8762, false, false},
    {"Summer24Prompt24_RunBCDEFGHI_V1", -2.2470, 1.7890, false, false},
    {"Summer24Prompt24_RunBCDEFGHI_V1", -1.9865, 2.4871, false, true},
    {"Summer24Prompt25_RunCDEFG_V1", -2.4110, 1.8762, false, false},
    {"Summer24Prompt25_RunCDEFG_V1", -2.2470, -1.3526, true, true},
    {"Summer24Prompt25_RunCDEFG_V1", -2.2470, 1.7890, false, false},
    {"Summer24Prompt25_RunCDEFG_V1", -2.1075, -2.3998, true, true},
    {"Summer24Prompt25_RunCDEFG_V1", -1.9865, 2.5744, false, true},
    {"Summer24Prompt26_RunBCD_V1", -2.4110, 1.8762, false, false},
    {"Summer24Prompt26_RunBCD_V1", -2.2470, -1.3526, true, true},
    {"Summer24Prompt26_RunBCD_V1", -2.2470, 1.7890, false, false},
    {"Summer24Prompt26_RunBCD_V1", -2.1075, -2.3998, true, true},
    {"Summer24Prompt26_RunBCD_V1", -1.9865, 2.5744, false, true},
    // Run 2 UL: calibration = recommended map (+ h2hot_mc for UL16)
    {"Summer19UL16_V1", -3.5765, -3.0107, true, true},
    {"Summer19UL16_V1", -2.4110, -1.3526, false, false},
    {"Summer19UL16_V1", -2.4110, 1.8762, false, false},
    {"Summer19UL16_V1", -1.0875, 0.4800, true, true},
    {"Summer19UL16_V1", 0.0435, 2.0508, false, true},
    {"Summer19UL17_V1", -2.9085, 2.3998, true, true},
    {"Summer19UL17_V1", -2.7515, 2.2253, true, true},
    {"Summer19UL17_V1", -2.4110, -1.3526, false, false},
    {"Summer19UL17_V1", -2.4110, 1.8762, false, false},
    {"Summer19UL18_V1", -2.9085, -1.5272, true, true},
    {"Summer19UL18_V1", -2.9085, -1.0908, true, true},
    {"Summer19UL18_V1", -2.4110, 1.8762, false, false},
    {"Summer19UL18_V1", -2.2470, 1.7890, false, false},
};

// pre-UL Run 2 TightLepVeto, transcribed from the twiki selection code
// (CMS/JetID13TeVRun2016, 2018) and table (2017); pat-style multiplicities
struct PatJet {
  double eta, CHF, NHF, CEF, NEF, MUF;
  int CHM, NM;
};

static bool TwikiID2016(const PatJet &j) {
  double a = std::abs(j.eta);
  int nConst = j.CHM + j.NM;
  if (a <= 2.7) {
    return j.NHF < 0.90 && j.NEF < 0.90 && nConst > 1 && j.MUF < 0.8 &&
           ((a <= 2.4 && j.CHF > 0 && j.CHM > 0 && j.CEF < 0.90) || a > 2.4);
  }
  if (a <= 3.0) {
    return j.NHF < 0.98 && j.NEF > 0.01 && j.NM > 2;
  }
  return j.NEF < 0.90 && j.NM > 10;
}

static bool TwikiID2017(const PatJet &j) {
  double a = std::abs(j.eta);
  int nConst = j.CHM + j.NM;
  if (a <= 2.7) {
    return j.NHF < 0.90 && j.NEF < 0.90 && nConst > 1 && j.MUF < 0.80 &&
           ((a <= 2.4 && j.CHF > 0 && j.CHM > 0 && j.CEF < 0.80) || a > 2.4);
  }
  if (a <= 3.0) {
    return j.NEF > 0.02 && j.NEF < 0.99 && j.NM > 2;
  }
  return j.NEF < 0.90 && j.NHF > 0.02 && j.NM > 10;
}

static bool TwikiID2018(const PatJet &j, bool ion) {
  double a = std::abs(j.eta);
  int nConst = j.CHM + j.NM;
  if (a <= 2.6) {
    return j.CEF < 0.8 && j.NEF < 0.9 && j.MUF < 0.8 &&
           (ion || (j.CHM > 0 && j.CHF > 0 && nConst > 1 && j.NHF < 0.9));
  }
  if (a <= 2.7) {
    return j.CEF < 0.8 && j.NEF < 0.99 && j.MUF < 0.8 &&
           (ion || (j.CHM > 0 && j.NHF < 0.9));
  }
  if (a <= 3.0) {
    return j.NEF > 0.02 && j.NEF < 0.99 && (ion || j.NM > 2);
  }
  if (a <= 5.0) {
    return j.NEF < 0.90 && (ion || (j.NHF > 0.2 && j.NM > 10));
  }
  return false;
}

static void TestID() {
  std::printf("jet ID\n");
  const std::string veto = VetoFile("Summer24Prompt24_RunBCDEFGHI_V1");
  JetSelector pp(JetSelector::System::pp, kJetID, veto);
  JetSelector ion(JetSelector::System::Ion, kJetID, veto);
  for (const Cut &c : IDCases()) {
    Check(ID(pp, c.eta, c.jet) == c.ppPass,
          std::string("pp ") + c.name + " eta " + std::to_string(c.eta));
    Check(ID(ion, c.eta, c.jet) == c.ionPass,
          std::string("ion ") + c.name + " eta " + std::to_string(c.eta));
  }
}

// random jets, JetSelector + json/jetid_Run*_preUL.json vs the twiki code
static void TestRun2ID() {
  std::printf("run 2 pre-UL jet ID vs twiki\n");
  const std::string veto = VetoFile("Summer19UL18_V1");
  struct Year {
    const char *name, *file;
    JetSelector::System sys;
    bool (*ref)(const PatJet &);
  };
  static auto ref2018pp = [](const PatJet &j) { return TwikiID2018(j, false); };
  static auto ref2018ion = [](const PatJet &j) { return TwikiID2018(j, true); };
  const Year years[] = {
      {"2016 pp", "json/jetid_Run2016_preUL.json", JetSelector::System::pp,
       TwikiID2016},
      {"2017 pp", "json/jetid_Run2017_preUL.json", JetSelector::System::pp,
       TwikiID2017},
      {"2018 pp", "json/jetid_Run2018_preUL.json", JetSelector::System::pp,
       +ref2018pp},
      {"2018 ion", "json/jetid_Run2018_preUL.json", JetSelector::System::Ion,
       +ref2018ion},
  };

  std::mt19937 rng(20260923);
  std::uniform_real_distribution<double> uEta(-5.5, 5.5), u01(0.0, 1.0);
  // exact zeros matter for the strict "> 0" cuts
  auto frac = [&]() { return u01(rng) < 0.15 ? 0.0 : u01(rng); };
  for (const Year &y : years) {
    JetSelector js(y.sys, y.file, veto);
    const int n = 200000;
    int mism = 0, pass = 0;
    for (int i = 0; i < n; i++) {
      double eta = uEta(rng);
      double CHF = frac(), NHF = frac(), CEF = frac(), NEF = frac(),
             MUF = frac();
      int CHM = rng() % 4, NHM = rng() % 7, CEM = rng() % 2, NEM = rng() % 15,
          MUM = rng() % 2;
      PatJet pj{eta, CHF, NHF, CEF, NEF, MUF, CHM + CEM + MUM, NHM + NEM};
      bool want = y.ref(pj);
      bool got = js.PassesID(eta, CHF, NHF, CEF, NEF, MUF, CHM, NHM, CEM, NEM,
                             MUM);
      mism += (want != got);
      pass += want;
    }
    Check(mism == 0, std::string(y.name) + ": " + std::to_string(mism) +
                         " of " + std::to_string(n) + " random jets differ");
    Check(pass > n / 100 && pass < n * 99 / 100,
          std::string(y.name) + ": both outcomes exercised");
  }
}

static void TestVetoMaps() {
  std::printf("veto maps\n");
  std::string current;
  JetSelector *ana = nullptr;
  JetSelector *cal = nullptr;
  for (const VetoPoint &p : kVetoPoints) {
    if (p.tag != current) {
      delete ana;
      delete cal;
      current = p.tag;
      ana = new JetSelector(JetSelector::System::pp, kJetID, VetoFile(p.tag));
      cal = new JetSelector(JetSelector::System::pp, kJetID, VetoFile(p.tag),
                            JetSelector::Purpose::Calibration);
      Check(ana->VetoMapTag() == p.tag, "tag read from " + VetoFile(p.tag));
      Check(ana->VetoMapType() == "jetvetomap", "analysis map type");
      Check(cal->VetoMapType() == "jetvetomap_all", "calibration map type");
    }
    std::string where = std::string(p.tag) + " eta " + std::to_string(p.eta) +
                        " phi " + std::to_string(p.phi);
    Check(ana->InVetoRegion(p.eta, p.phi) == p.analysis, "analysis " + where);
    Check(cal->InVetoRegion(p.eta, p.phi) == p.calibration,
          "calibration " + where);
  }
  delete ana;
  delete cal;

  // outside the map returns the file's flow value, 0, not vetoed
  JetSelector js(JetSelector::System::pp, kJetID,
                 VetoFile("Summer24Prompt24_RunBCDEFGHI_V1"));
  Check(!js.InVetoRegion(5.3, 0.0), "eta beyond map not vetoed");
  Check(!js.InVetoRegion(-5.3, 0.0), "eta below map not vetoed");
}

static void TestSelectionAndEventVeto() {
  std::printf("selection + event veto\n");
  JetSelector js(JetSelector::System::pp, kJetID,
                 VetoFile("Summer24Prompt24_RunBCDEFGHI_V1"));
  Jet j;
  const double vEta = -3.9260, vPhi = 2.5744; // vetoed
  const double cEta = -2.4110, cPhi = 1.8762; // clean
  auto sel = [&](double eta, double phi, const Jet &x) {
    return !js.VetoJet(eta, phi, x.CHF, x.NHF, x.CEF, x.NEF, x.MUF, x.CHM,
                       x.NHM, x.CEM, x.NEM, x.MUM);
  };
  auto ev = [&](double pt, double eta, double phi, const Jet &x) {
    return js.VetoEvent(pt, eta, phi, x.CHF, x.NHF, x.CEF, x.NEF, x.MUF,
                        x.CHM, x.NHM, x.CEM, x.NEM, x.MUM);
  };

  Check(sel(cEta, cPhi, j), "good jet, clean region selected");
  Check(!sel(vEta, vPhi, j), "good jet, vetoed region rejected");
  Jet bad = j;
  bad.MUF = 0.9;
  Check(!sel(cEta, cPhi, bad), "failing ID rejected");

  Check(ev(20.0, vEta, vPhi, j), "event vetoed by jet in vetoed region");
  Check(!ev(20.0, cEta, cPhi, j), "clean region never vetoes event");
  Check(!ev(15.0, vEta, vPhi, j), "pT 15 does not veto");
  Check(ev(15.01, vEta, vPhi, j), "pT 15.01 vetoes");
  // lepton cuts stop at |eta| 2.7, fail the forward jet on NEF instead
  Jet badFwd = j;
  badFwd.NEF = 0.5;
  Check(!ev(20.0, vEta, vPhi, badFwd), "failing ID does not veto");
  Jet em = j;
  em.CEF = 0.0;
  em.NEF = 0.3;
  Check(ev(20.0, vEta, vPhi, em), "EM fraction 0.3 vetoes");
  // forward ID needs NEF < 0.4, so test the EM cut in the central region
  const double bEta = -2.2470, bPhi = -1.3526; // vetoed, central
  em.CEF = 0.05;
  em.NEF = 0.8;
  Check(ev(20.0, bEta, bPhi, em), "EM fraction 0.85 vetoes");
  em.CEF = 0.2;
  em.NEF = 0.75;
  Check(sel(cEta, cPhi, em), "EM fraction 0.95 jet still passes ID");
  Check(!ev(20.0, bEta, bPhi, em), "EM fraction 0.95 does not veto");

  // event-level overload on forest-style arrays
  float eta[3] = {(float)cEta, (float)vEta, (float)cEta};
  float phi[3] = {(float)cPhi, (float)vPhi, (float)cPhi};
  double pt[3] = {50.0, 20.0, 30.0};
  float CHF[3], NHF[3], CEF[3], NEF[3], MUF[3];
  int CHM[3], NHM[3], CEM[3], NEM[3], MUM[3];
  for (int k = 0; k < 3; k++) {
    CHF[k] = j.CHF;
    NHF[k] = j.NHF;
    CEF[k] = j.CEF;
    NEF[k] = j.NEF;
    MUF[k] = j.MUF;
    CHM[k] = j.CHM;
    NHM[k] = j.NHM;
    CEM[k] = j.CEM;
    NEM[k] = j.NEM;
    MUM[k] = j.MUM;
  }
  auto evArr = [&](int n) {
    return js.VetoEvent(n, pt, eta, phi, CHF, NHF, CEF, NEF, MUF, CHM, NHM,
                        CEM, NEM, MUM);
  };
  Check(evArr(3), "event with one vetoing jet among three is vetoed");
  Check(!evArr(1), "event with only the clean first jet is kept");
  pt[1] = 10.0;
  Check(!evArr(3), "vetoed-region jet below 15 GeV doesn't veto the event");
  Check(!evArr(0), "empty event is kept");
}

static void TestErrors() {
  std::printf("errors\n");
  auto throws = [](auto f) {
    try {
      f();
    } catch (const std::runtime_error &) {
      return true;
    }
    return false;
  };
  const std::string veto = VetoFile("Summer24Prompt24_RunBCDEFGHI_V1");
  Check(throws([&] {
          JetSelector js(JetSelector::System::pp, "json/nope.json", veto);
        }),
        "missing jet ID file throws");
  Check(throws([&] {
          JetSelector js(JetSelector::System::pp, kJetID, "json/nope.json");
        }),
        "missing veto file throws");
  Check(throws([&] {
          JetSelector js(JetSelector::System::pp, veto, veto);
        }),
        "veto file as jet ID file throws");
  Check(throws([&] {
          JetSelector js(JetSelector::System::pp, kJetID, kJetID);
        }),
        "jet ID file as veto file throws");
}

int main() {
  TestID();
  TestRun2ID();
  TestVetoMaps();
  TestSelectionAndEventVeto();
  TestErrors();
  std::printf("%d/%d checks passed\n", nCheck - nFail, nCheck);
  return nFail == 0 ? 0 : 1;
}
