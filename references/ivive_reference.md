# IVIVE Reference: In Vitro to In Vivo Extrapolation Scaling Factors

## Hepatic Clearance Scaling

### Microsomal scaling (human)

| Parameter | Value | Unit |
|---|---|---|
| MPPGL (microsomal protein per gram liver) | 40 | mg/g |
| Liver weight | 25.7 (male), 24.1 (female) | g/kg body weight |
| Hepatic blood flow (Qh) | 87 (male), 91 (female) | L/h |
| Hepatocyte number | 110 | million cells/g liver |

**Calculation:**
```
CLint_hep = in_vitro_clint [uL/min/mg microsomal protein] * MPPGL * liver_weight * body_weight
CLint_hep [L/h] = Clint_uL_min_mg * 40 * 25.7 * BW / 1000 * 60 / 1000

CLhep (well-stirred) = Qh * fu * CLint / (Qh + fu * CLint)
CLhep (parallel tube) = Qh * (1 - exp(-fu * CLint / Qh))
CLhep (dispersion) = Qh * (1 - 4a / (1+a)^2)  where a = sqrt(1 + 4*fu*CLint/Qh/Dn)
```

### Hepatocyte scaling

```
CLint_hep = in_vitro_clint [uL/min/10^6 cells] * hepatocyte_number * liver_weight * BW
```

## Species-Specific Scaling Factors

| Species | MPPGL (mg/g) | Liver wt (g/kg) | Qh (L/h/kg) | Body weight (kg) |
|---|---|---|---|---|
| Human | 40 | 25.7 | 21 | 70 |
| Rat | 45 | 40 | 55 | 0.25 |
| Mouse | 55 | 43 | 90 | 0.025 |
| Dog | 35 | 30 | 41 | 10 |
| Monkey | 30 | 27 | 44 | 5 |
| Minipig | 35 | 25 | 36 | 20 |
| Rabbit | 40 | 27 | 46 | 3 |

## Intestinal Absorption Scaling

### Caco-2 to Peff correlation

```
Peff (cm/s * 10^-4) = 0.54 * log10(Caco2_Papp) + 1.56  (Caco2 in cm/s)
```

Reference: Sugano et al. (2010)

### Fa (fraction absorbed) from Peff

```
Fa = 1 - exp(-2 * Peff * tau / R)  (plug flow model)
tau = residence time in small intestine (~3.5 h)
R = small intestine radius (~1.75 cm)
```

Or simplified: Fa = 1 - exp(-Peff_human * An)
where An = dose number, depends on solubility and permeability.

## Renal Clearance Scaling

```
CLrenal = fu * (GFR * ffiltration + CLR_secretion) - CLR_reabsorption
```

| Species | GFR (mL/min/kg) |
|---|---|
| Human | 1.8 |
| Rat | 5.2 |
| Mouse | 9.6 |
| Dog | 3.2 |
| Monkey | 2.0 |

## CYP Abundance Scaling

| Enzyme | pmol/mg microsomal protein (human) | Inter-individual CV |
|---|---|---|
| CYP1A2 | 52 | 70% |
| CYP2C8 | 24 | 60% |
| CYP2C9 | 96 | 60% |
| CYP2C19 | 19 | 90% |
| CYP2D6 | 8 | 80% |
| CYP3A4 | 137 | 60% |
| CYP3A5 | 21 (expressors) | 100% |

**IVIVE for metabolic clearance:**
```
CLint_CYP = Vmax/Km  (or Clint per enzyme)
CLint_total = sum(CLint_CYPi * abundance_i) for all contributing CYPs
```

## Blood-to-Plasma Ratio

```
CLblood = CLplasma / B:P ratio
fu_blood = fu_plasma / B:P ratio
```

Typical B:P ratios: acidic drugs ~0.55, basic drugs ~1.2-1.5, neutral ~1.0

## Key References

1. Barter et al. (2007) Scaling factors for drug metabolism
2. Sugano et al. (2010) Biopharmaceutics modeling
3. Ito & Houston (2005) Prediction of human drug clearance
4. OSP Documentation: https://docs.open-systems-pharmacology.org
