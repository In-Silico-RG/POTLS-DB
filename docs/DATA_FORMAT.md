# DATA_FORMAT.md — POTLS-DB Structure Files

## Overview

All generated structures are stored in JSON format, one file per LGS
configuration. Files are organized as JSON arrays, where each element
represents one lignin oligomer.

---

## File naming convention

```
POTLS_{BIOMASS}_{METHOD}_{SET}.json
```

| Token    | Values |
|----------|--------|
| BIOMASS  | EFB, MF, KS |
| METHOD   | KR (kraft), SD (soda), OR (organosolv) |
| SET      | C (central), L (sensitivity_low), H (sensitivity_high) |

Examples:
- `POTLS_EFB_KR_C.json`  — EFB Kraft, central parameters, QC-passed
- `POTLS_KS_OR_H_rejected.json` — KS Organosolv, high sensitivity, QC-failed

---

## JSON record schema

Each record contains the following fields:

```json
{
  "lg_id":           "POTLS_EFB_KR_C_00001",
  "config_id":       "POTLS_EFB_KR_C",
  "smiles":          "COc1cc(...) ccc1O",
  "DP":              13,
  "MW_lgs":          2641.5,
  "MW_rdkit":        2639.8,
  "HC_rdkit":        1.108,
  "OC_rdkit":        0.351,
  "SG_rdkit":        2.27,
  "n_phenolic_OH":   3,
  "n_methoxy":       24,
  "bond_frequencies": {
    "beta-O-4":  0.731,
    "beta-beta": 0.118,
    "beta-5":    0.083,
    "5-5":       0.012,
    "4-O-5":     0.048,
    "DBDO":      0.008
  },
  "monomer_count": {
    "G": 4,
    "S": 9,
    "H": 0
  },
  "pass_F1":       true,
  "pass_F2":       true,
  "pass_F3":       true,
  "F2_deviation":  0.013,
  "F3_deviation":  0.0,
  "qc_passed":     true
}
```

### Field definitions

| Field | Type | Description |
|-------|------|-------------|
| `lg_id` | string | Unique structure identifier |
| `config_id` | string | Source LGS configuration |
| `smiles` | string | Canonical SMILES (RDKit-validated) |
| `DP` | int | Degree of polymerization (monomer count) |
| `MW_lgs` | float | Molecular weight reported by LGS (g/mol) |
| `MW_rdkit` | float | Exact MW calculated by RDKit (g/mol) |
| `HC_rdkit` | float | Atomic H/C ratio (RDKit) |
| `OC_rdkit` | float | Atomic O/C ratio (RDKit) |
| `SG_rdkit` | float | S/G ratio from SMARTS matching (RDKit) |
| `n_phenolic_OH` | int | Free phenolic OH groups (SMARTS: `[OH]c`) |
| `n_methoxy` | int | Methoxy groups |
| `bond_frequencies` | dict | Normalized interunit linkage fractions (sum = 1.0) |
| `monomer_count` | dict | Absolute count of G, S, H monomers |
| `pass_F1` | bool | F1 filter: valid SMILES + DP in range |
| `pass_F2` | bool | F2 filter: S/G within ±5% of target |
| `pass_F3` | bool | F3 filter: MW within calibrated range ±20% |
| `F2_deviation` | float | Absolute fractional deviation of S/G from target |
| `F3_deviation` | float | Fractional deviation of MW from calibration window |
| `qc_passed` | bool | True if F1 ∧ F2 ∧ F3 all pass |

---

## Convergence files

One CSV per configuration in `data/convergence/`:

```
POTLS_{ID}_convergence.csv
```

| Column | Description |
|--------|-------------|
| `batch` | Batch number |
| `n_batch` | Structures generated in this batch |
| `mean_MW` | Cumulative mean MW of passed structures |
| `std_MW` | Cumulative standard deviation of MW |
| `mean_CX` | Mean complexity index (condensed linkage fraction) |
| `std_CX` | Std deviation of complexity index |
| `delta_mean_MW` | Relative change in mean MW from previous batch |
| `delta_std_MW` | Relative change in std MW from previous batch |
| `delta_mean_CX` | Relative change in mean CX from previous batch |
| `delta_std_CX` | Relative change in std CX from previous batch |
| `converged` | True when all four deltas < 0.01 |

---

## Metadata table

`POTLS_DB_metadata.csv` — one row per configuration:

| Column | Description |
|--------|-------------|
| `dataset_id` | Configuration identifier |
| `biomass_source` | EFB / MF / KS |
| `pulping_method` | Kraft / Soda / Organosolv |
| `parameter_set` | central / sensitivity_low / sensitivity_high |
| `archetype` | SG-rico / SG-intermedio / G-dominante |
| `SG_target` | Target S/G ratio from HSQC |
| `DP_min`, `DP_max` | Calibrated DP window |
| `MW_min`, `MW_max` | Corresponding MW window (g/mol) |
| `f_betaO4` | Central β-O-4 linkage fraction |
| `n_generated` | Total structures generated |
| `n_passed` | QC-passed structures |
| `pass_rate` | n_passed / n_generated |
| `n_fail_F1/F2/F3` | Failures per filter |
| `converged` | Whether convergence criterion was met |
| `mean_MW_final` | Final converged mean MW |
| `std_MW_final` | Final converged std MW |
| `run_timestamp` | ISO timestamp of generation run |
