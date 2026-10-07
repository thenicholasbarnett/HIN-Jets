
***Mission Statement:*** *Help analyzers in the HIN PAG select, correct, and smear jets with JME reccomendations using a standalone codebase.*

| File | Use |
| - | - |
| `include/JetCorrector.h` | Get corrected, and up & down variations, of jet p<sub>T</sub> |
| `include/JetSmearer.h` | Get smeared, and up & down variations, of jet p<sub>T</sub>; Match two jet collections |
| `include/JetSelector.h` | Apply jet ID, jet veto map, and JME recommended event veto |

Classes written in the header files matching the glob pattern `include/Jet*.h` in this repository correct JES, smear JER, provide JES and JER uncertainties, as well as veto, match, and sort jets.

Jet analysis workflows should follow this order of operations:
1. Correct - correct jet p<sub>T</sub>
2. Uncertainty - get up and down variation in jet p<sub>T</sub> due to JES uncerainty
3. Match - match reconstructed and generated jets
4. Smear - smear simulated JER to resemble data
2. Uncertainty - get up and down variation in jet p<sub>T</sub> due to JER uncerainty
5. Sort - order jets by descending p<sub>T</sub>
6. Select - select jets and events with jet ID and veto map JSON

# JetCorrector

JetCorrector class in `include/JetCorrector.h` corrects the energy scale of jets as well as provides the uncertainty. Get the corrected p<sub>T</sub>, with the up and down variations, for a jet using this class.

## Usage

Examples below is for ak4 jets from 2018 PbPb in data.

### 

Include the header and instantiate JetCorrector object.

```cpp
#include "JetCorrector.h"

JetCorrector JES({"txt/2018PbPb/Autumn18_HI_V8_MC_L2Relative_AK4PF.txt",
                  "txt/2018PbPb/Autumn18_HI_V8_DATA_L2L3Residual_AK4PF.txt"},
                  "txt/2018PbPb/Autumn18_HI_V8_DATA_Uncertainty_AK4PF.txt");
```

### Correct Jet Energy Scale

Scale a jet's transverse momentum.

> NOTE: returns jet p<sub>T</sub> of -1 outside of range in correction plain text file(s) applied.

#### Inline Arguments

```cpp
double ptCorr = JES.CorrectedPt(rawpt[j], jteta[j], jtphi[j]);
double ptUp = JES.CorrectedPt(rawpt[j], jteta[j], jtphi[j], Variation::UP); // ptCorr * (1 + up)
double ptDown = JES.CorrectedPt(rawpt[j], jteta[j], jtphi[j], Variation::DOWN); // ptCorr * (1 - down)
std::pair<double, double> unc = JES.Uncertainty(rawpt[j], jteta[j], jtphi[j]); // {down, up}
```

#### Setters and Getters
```cpp
// input for correction factor
JES.SetJetPT(rawpt[j]);
JES.SetJetEta(jteta[j]);
JES.SetJetPhi(jtphi[j]);

// output corrected jet pT and variations
double ptCorr = JES.CorrectedPt();
double ptUp = JES.CorrectedPt(Variation::UP); // ptCorr * (1 + up)
double ptDown = JES.CorrectedPt(Variation::DOWN); // ptCorr * (1 - down)
std::pair<double, double> unc = JES.Uncertainty(); // {down, up}
```

### Quantify Jet Energy Scale Uncertainty

Get uncertainty, as well as up and down variations, in corrected jet p<sub>T</sub>.
```cpp
std::pair<double, double> unc = JES.GetUncertainty(); // {down, up}
double ptUp = JES.GetCorrectedPT(Variation::UP); // ptNom * (1 + up)
double ptDown = JES.GetCorrectedPT(Variation::DOWN); // ptNom * (1 - down)
```

JetUncertainty class can also give uncertainty with corrected jet p<sub>T</sub>.
```cpp
JetUncertainty JEU("txt/2023PbPb/Spring23PbPb_TotalUncertainties.txt");
JEU.SetJetPT(ptCorr);
JEU.SetJetEta(jteta[j]);
JEU.SetJetPhi(jtphi[j]);
std::pair<double, double> unc = JEU.GetUncertainty();
double ptDown = ptCorr * (1 - unc.first);
double ptUp = ptCorr * (1 + unc.second);
```
> NOTE: GetUncertainty() returns uncertainty of 1 outside of η range in uncertainty plain text file used.

# JetSmearer
Smear energy resolution of simulated jets with the JetSmearer class in `include/JetSmearer.h`. Get smeared p<sub>T</sub>, and its variations, for a jet with using this class.

## Usage

Example below is for ak4 jets from 2024 ppRef MC.

Include header.
```cpp
#include "JetSmearer.h"
```
Instantiate class.
```cpp
// input: relative pT resolution in MC, JER scale factors
JetSmearer smearer("RelativePtResolutionMC_AK4PF.txt", "JERSF_AK4PF.txt");
```
### Match Reconstructed to Generated Jets

### Smear Simulated Jet Energy Resolution
Get a nominal smeared jet p<sub>T</sub>, and p<sub>T</sub> variations, with JetSmearer instance.
```cpp
double ptSmeared = smearer.SmearedPt(ptCorr, jteta[j], rho, genPt, evt);
double ptUp = smearer.SmearedPt(ptCorr, jteta[j], rho, genPt, evt, Variation::UP);
double ptDown = smearer.SmearedPt(ptCorr, jteta[j], rho, genPt, evt, Variation::DOWN);
```
### Smearing Information
Details for any jet smearing can be obtained. Only possible after smearing with default, or other deterministic, RNG seed.
```
JetSmearing::Result r = smearer.Smear(ptCorr, jteta[j], rho, genPt, evt);
r.smearFactor // smeared p<sub>T</sub> = ptCorr * r.smearFactor
r.resolution // sigmaJER
r.scaleFactor // JER SF
r.matched // true = scaling, false = stochastic (for Hybrid/default)
```
Users can optionally provide the random number generator seed and/or specify to use scaling, stochastic, or hybrid procedure. JetSmearer uses hybrid smearing by default, which is the official recommendation of JME. Please see [[JERC - Application Methods for JER](https://cms-jme-jerc.docs.cern.ch/application/jer/)] for more details on these JER smearing implementations. 
JetSmearer v2.1, in the default setting, and correctionlib 2.9.0 have given bit-identical output for the same input smearing factor JSON and 200k jet sample.

## Select Jets
`include/JetSelector.h`

### Apply Jet Veto Maps

### Apply Jet ID Criteria

### Veto Events 

## Mapping Jets

JetMapper namespace from `include/JetMapper.h` has functions to match reconstructed to generated jets and sort jets by descending p<sub>T</sub>. 
```cpp
#include "JetMapper.h"
```

### Matching Jets
Given jet kinematics and multiplicity for an event JetMapper::Match returns an output vector of integers with nref entries. Each index in the output is the index of the reconstructed jet, and each entry in the output is the index to the matched generated jet. Each reconstructed jet without a matched generated jet has -999 in its output entry. 

Example below using one-to-one matching by ΔR closest proximity within 0.2 in η-ϕ space.
```cpp
std::vector<int> match = JetMapper::Match(nref, jteta, jtphi, ngen, geneta, genphi, 0.2)
```

Example below uses official JME recommendation that the following two matching conditions are satisfied.
1. ΔR < Rcone/2
2. |pTJEC - pTgen| < 3 σJER pTJEC is satisfied.
```cpp
```

### Ordering Jets

# HIN Datasets

## Run 2

| Dataset | √s_NN [TeV] | Era | Run(s) | JEC | JER | JES Uncertainty | Jet Veto Map |
| - | - | - | - | - | - | - | - |
| 2016 pPb | 8.16 | PARun2016C | 285479-286496 | `Autumn16_HI_pPb_{pgoing,Pbgoing}_Embedded_MC_L2Relative_AK4PF`, `Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PF` (pp 2016 Run H) | `Summer16_25nsV1_MC_*_AK4PF` (pp 2016) | `Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PF` | `jetvetomaps_Summer19UL16_V1.json` |
| 2017 ppRef | 5.02 | 2017G | 306546-306826 | `Spring18_ppRef5TeV_V6_MC_L2Relative_AK{2,3,4,5,6}PF`, `Spring18_ppRef5TeV_V6_DATA_L2L3Residual_AK{2,3,4,5,6}PF` | `Fall17_V3_MC_*_AK4PF` (pp 2017, 94X) | `Spring18_ppRef5TeV_V6_DATA_Uncertainty_AK{2,3,4,5,6}PF` | `jetvetomaps_Summer19UL17_V1.json` |
| 2018 PbPb | 5.02 | HIRun2018A | 326381-327564 | `Autumn18_HI_V8_MC_L2Relative_AK{2,3,4,5,6}PF`, `Autumn18_HI_V8_DATA_L2L3Residual_AK{2,3,4,5,6}PF` | `Autumn18_RunD_V7b_MC_*_AK4PF` (pp 2018) | `Autumn18_HI_V8_DATA_Uncertainty_AK{2,3,4,5,6}PF` | `jetvetomaps_Summer19UL18_V1.json` |

## Run 3

| Dataset | √s_NN [TeV] | Era | Run(s) | JEC | JER | JES Uncertainty | Jet Veto Map |
| - | - | - | - | - | - | - | - |
| 2023 ppRef | 5.36 | Run2023F | 373710 | `L2Relative_AK{2,3,4,5,6}PF_ppReco_v1.txt`, `L2Residuals_2023ppRef_ppReco_v0.txt` | n/a | n/a | `Summer23BPixPrompt23_RunD_V1` |
| 2023 PbPb | 5.36 | HIRun2023A | 374288–375823 | `Prompt23HIPbPb_V1_DATA_L2Residual_AK{2,3,4,5,6}PF.txt`, `Prompt23HIPbPb_V1_MC_L2Relative_AK{2,3,4,5,6}PF.txt` | n/a | `Spring23PbPb_TotalUncertainties.txt` | `Summer23BPixPrompt23_RunD_V1` |
| 2024 lowPUpp | 13.6 | Run2024I | 386{642,749,753} | ? | ? | n/a | `Summer24Prompt24_RunBCDEFGHI_V1` |
| 2024 ppRef | 5.36 | Run2024J | 387474-387721 | `Prompt24HIpp_V1_MC_L1FastJet_AK{2,3,4,5,6}PF.txt`, `Prompt24HIpp_V2_MC_L2Relative_AK{2,3,4,5,6}PF.txt`, `Prompt24HIpp_V1_DATA_L2Residual_AK{2,3,4,5,6}PF.txt` | ? | n/a | `Summer24Prompt24_RunBCDEFGHI_V1` |
| 2024 PbPb | 5.36 | HIRun2024{A,B} | 387853-388784 | ? | n/a | n/a | `Summer24Prompt24_RunBCDEFGHI_V1` |
| 2025 OO | 5.36 | OORun2025 | 394153-394217 | `Prompt25HIOO_V1_MC_L2Relative_AK{2,3,4,5,6}PF.txt`, `2026_05_25_OO_cent_80_90_run3_AK{2,4}_ptFitnokFSR.txt` | n/a | n/a | `Summer24Prompt25_RunCDEFG_V1` |
| 2025 lowPUpp | 13.6 | Run2025G | 398{682,683,803} | n/a | n/a | n/a | `Summer24Prompt25_RunCDEFG_V1` |
| 2025 PbPb | 5.36 | HIRun2025A | 399465-400426 | n/a | n/a | n/a | `Summer24Prompt25_RunCDEFG_V1` |
| 2026 lowPUpp | 13.6 | Run2026D | 40386{3,6} | n/a | n/a | n/a | `Summer24Prompt26_RunBCD_V1` |
| 2026 PbPb | 5.36 | HIRun2026A | 404469-404926 | n/a | n/a | n/a | `Summer24Prompt26_RunBCD_V1` |