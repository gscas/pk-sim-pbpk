---
name: pk-sim-pbpk
description: |
  PK-Sim PBPK modeling and simulation skill for Open Systems Pharmacology.
  Use when: (1) Building PBPK models from compound physicochemical/ADME data,
  (2) Running PK simulations (IV/PO, single/multiple dose, DDI, population),
  (3) Performing parameter identification/optimization against observed data,
  (4) Setting up drug-drug interaction simulations (CYP inhibition/induction),
  (5) Creating population PK simulations (special populations, covariates),
  (6) Sensitivity analysis of PBPK model parameters,
  (7) Converting in vitro data to in vivo inputs (IVIVE),
  (8) Generating PK parameters (AUC, Cmax, Tmax, half-life, clearance, Vd),
  (9) Exporting simulation results and creating PK reports.
  Covers both PK-Sim GUI workflow and CLI/Python API automation.
---

# PK-Sim PBPK Modeling and Simulation

PK-Sim is a whole-body physiologically based pharmacokinetic (PBPK) modeling tool
from Open Systems Pharmacology (OSP). This skill covers model building, simulation,
parameter optimization, DDI, population PK, and result analysis.

## Prerequisites

- PK-Sim installed (https://github.com/Open-Systems-Pharmacology/PK-Sim/releases)
- Python 3.8+ with ospsuite package: pip install ospsuite
- Or R with ospsuite package (see references/r-api.md)

## Workflow Overview

1. Define Compound
2. Define Individual
3. Define Protocol
4. Create Simulation
5. Run and Compare
6. Optimize
7. Export

## 1. Define Compound

Required parameters (minimum viable model):

| Parameter | Unit | Typical Range | Notes |
|---|---|---|---|
| Molweight | g/mol | 100-1000 | Molecular weight |
| Lipophilicity (logP/logMA) | - | -2 to 6 | Prefer membrane affinity (logMA) |
| Fraction unbound (fu) | - | 0.001-1 | Species-specific |
| Solubility at ref pH | mg/L | varies | For oral absorption |
| pKa | - | 1-14 | Acid/base; affects dissolution |

ADME properties:

| Process | Input | Method |
|---|---|---|
| Hepatic clearance | in vitro Clint (uL/min/mg) | IVIVE with well-stirred model |
| Renal clearance | GFR fraction or measured CLrenal | Scaling |
| CYP inhibition | Ki, IC50 | Competitive/mechanism-based |
| CYP induction | EC50, Emax | Link to expression profile |
| Transport | Km, Vmax or CLint | Active uptake/efflux |
| Permeability | Caco-2 Papp or specific permeability | Intestinal absorption |

Model type selection:
- Small molecule: 4-compartment per organ (plasma/interstitial/intracellular/blood cells)
- Large molecule (mAb/protein): 2-pore model with FcRn recycling

## 2. Define Individual / Population

Species: Human, Mouse, Rat, Minipig, Dog, Monkey, Beagle, Rabbit

Individual parameters:
- Age, Weight, Height, Gender
- Population (European, American, Japanese, etc.)
- Disease state (renal impairment, hepatic impairment)

Population simulation:
- Define covariate distributions (weight, age, enzyme activity)
- Number of individuals (typically 100-1000)
- Ontogeny functions auto-applied for CYPs/UGTs in children

Special populations:
- Renal impairment: reduce GFR, adjust organ volumes
- Hepatic impairment: reduce CYP activities (Child-Pugh scaling)
- Pediatric: ontogeny-based maturation of enzymes/transporters
- Elderly: adjust organ blood flow, GFR decline
- Pregnancy: specific PBPK model with placental transfer

## 3. Define Administration Protocol

| Route | Key Parameters |
|---|---|
| Intravenous (IV) | Dose [mg/kg], infusion time |
| Oral (PO) | Dose, formulation, absorption model |
| Subcutaneous | Dose, depot model |
| Inhalation | Dose, deposited fraction |

Formulation types for oral:
- Dissolved: immediate release
- Weibull: controlled release (dissolution curve f = 1 - exp(-alpha*t^beta))
- Particle dissolution: particle size dependent
- Table: user-defined dissolution profile

Multiple dosing: set dosing interval and number of doses; steady-state assessment requires >=5 half-lives

## 4. Create and Run Simulation

### Via PK-Sim CLI (Windows)

PKSim.CLI run --project project.pksim5 --simulation Sim1 --output results.csv
PKSim.CLI export --project project.pksim5 --simulation Sim1 --output sim1.pkml
PKSim.CLI snap --project project.pksim5 --simulation Sim1 --output sim1_snap.json

### Via Python (ospsuite)

import ospsuite

sim = ospsuite.Simulation("path/to/simulation.pkml")
sim.setParameterValue("Compound1|Lipophilicity", 2.5)
sim.setParameterValue("Compound1|Fraction unbound (plasma, reference value)", 0.1)
sim.run()

results = sim.results
time = results.time
conc_plasma = results["Organism|PeripheralVenousBlood|Compound1|Plasma (Peripheral Venous Blood)"]

### Via R (ospsuite)

See references/r-api.md for R workflow.

## 5. Parameter Identification (Optimization)

Fit model parameters to observed data.

Common parameters to optimize:
1. Lipophilicity (distribution)
2. Specific clearance / Clint (elimination)
3. Specific intestinal permeability (absorption)
4. Fraction unbound (if uncertain)
5. kcat for enzymatic metabolism

Optimization methods: Nelder-Mead, Levenberg-Marquardt

## 6. Drug-Drug Interaction (DDI) Simulation

Steps:
1. Build well-validated perpetrator model (with inhibition/induction parameters)
2. Build victim (substrate) model with metabolic pathway
3. Define expression profiles for involved enzymes
4. Create DDI simulation with both compounds
5. Run with and without perpetrator; calculate AUCR and CmaxR

CYP inhibition input: Ki value (uM), inhibition type (competitive, noncompetitive, uncompetitive, mechanism-based). For MBI: Kinact, KI.

CYP induction input: Emax (fold induction), EC50 (uM).

DDI assessment criteria:
- AUCR >= 1.25: weak DDI
- AUCR >= 2.0: moderate DDI
- AUCR >= 5.0: strong DDI

## 7. Sensitivity Analysis

Identify parameters with greatest impact on PK outputs. Vary each parameter by +-10% and compute sensitivity coefficient.

## 8. PK Parameter Calculation

From concentration-time data (simulated or observed):

| Parameter | Formula | Unit |
|---|---|---|
| Cmax | max(C) | mg/L or uM |
| Tmax | t at Cmax | h |
| AUC0-t | Trapezoidal rule | mg*h/L |
| AUC0-inf | AUC0-t + Ct/lambdaz | mg*h/L |
| t1/2 | ln(2)/lambdaz | h |
| CL | Dose/AUC | L/h |
| Vd | CL/lambdaz | L |
| Vss | From moment analysis | L |

Use scripts/calculate_pk.py for automated PK parameter calculation from simulation CSV output.

## 9. IVIVE: In Vitro to In Vivo Extrapolation

Convert in vitro assay data to PBPK model inputs.

Hepatic clearance (well-stirred model):
  CLhep = Qh * fu * Clint / (Qh + fu * Clint)
  Clint = in_vitro_clint * MPPGL * liver_weight
  MPPGL = microsomal protein per gram liver (~40 mg/g human)

Intestinal permeability: Peff = f(fu_gut, Caco-2_Papp) correlation-based scaling.

Renal clearance: CLrenal = fu * GFR + CLR_secretion - CLR_reabsorption

See references/ivive_reference.md for detailed scaling factors by species.

## 10. Export and Reporting

Export formats:
- PKML: simulation interchange format
- CSV: results data
- JSON snapshot: version-control friendly

Report template:
1. Compound physicochemical properties table
2. ADME characterization
3. Model structure and assumptions
4. Observed vs predicted overlay plots
5. PK parameter comparison table
6. DDI prediction results (if applicable)
7. Sensitivity analysis tornado plot
8. Population simulation statistics

## Project File Conventions

- PK-Sim project: .pksim5 (binary, all building blocks)
- Simulation export: .pkml (XML, single simulation)
- Snapshot: .json (human-readable, version-control friendly)
- Observed data: .csv or .xls (time, concentration, error)
- Results: .csv (time, all outputs)

## Key References

- PK-Sim docs: https://docs.open-systems-pharmacology.org/working-with-pk-sim/pk-sim-documentation
- OSP Python API: https://github.com/Open-Systems-Pharmacology/OSPSuite.Python
- OSP R API: https://github.com/Open-Systems-Pharmacology/OSPSuite.R
- PK-Sim releases: https://github.com/Open-Systems-Pharmacology/PK-Sim/releases
- Community forum: https://github.com/Open-Systems-Pharmacology/Forum
