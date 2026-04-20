# PARAMETER_BASIS.md — Rationale for LGS Input Parameters

## 1. S/G molar ratios

Source: quantitative HSQC-NMR (Oñate-Gutiérrez et al. 2025)
Integration region: aromatic C2-H2 correlations of S (δC/δH 104/6.67)
and G (δC/δH 115/6.78) units.

| Configuration | S/G exp. | f_G   | f_S   |
|---------------|----------|-------|-------|
| EFB Kraft     | 2.30     | 0.303 | 0.697 |
| EFB Soda      | 1.50     | 0.400 | 0.600 |
| EFB Org       | 0.40     | 0.714 | 0.286 |
| MF Kraft      | 0.70     | 0.588 | 0.412 |
| MF Soda       | 0.00     | 1.000 | 0.000 |
| MF Org        | 0.20     | 0.833 | 0.167 |
| KS Kraft      | 0.10     | 0.909 | 0.091 |
| KS Soda       | 0.00     | 1.000 | 0.000 |
| KS Org        | 0.10     | 0.909 | 0.091 |

## 2. Interunit linkage proportions

Source: Eswaran et al. 2022 (Sci. Data), Table 1. Values represent
the biosynthetic architecture of the native lignin polymer, not the
post-pulping technical lignin structure. DP calibration (see §3)
captures the fragmentation induced by pulping.

Three archetypes used, based on S/G and biological classification
of Elaeis guineensis biomass fractions:

SG-rico (EFB Kraft, EFB Soda): hardwood-like (high f_S)
  β-O-4 central: 0.74 | range: 0.68–0.79
  β-β central:   0.12 | range: 0.04–0.19
  β-5 central:   0.08 | range: 0.04–0.12
  5-5 central:   0.01 | range: 0.00–0.04
  4-O-5 central: 0.05 | range: 0.02–0.08
  DBDO central:  0.00 | range: 0.00–0.01

SG-intermedio (EFB Org, MF Kraft): hardwood/grass interpolated
  β-O-4 central: 0.71 | range: 0.65–0.77
  β-β central:   0.08 | range: 0.03–0.15
  β-5 central:   0.11 | range: 0.07–0.14
  5-5 central:   0.03 | range: 0.01–0.06
  4-O-5 central: 0.05 | range: 0.03–0.07
  DBDO central:  0.01 | range: 0.00–0.02

G-dominante (MF Soda/Org, KS all): softwood-like (low/zero f_S)
  β-O-4 central: 0.64 | range: 0.59–0.69
  β-β central:   0.05 | range: 0.03–0.08
  β-5 central:   0.14 | range: 0.10–0.17
  5-5 central:   0.07 | range: 0.03–0.13
  4-O-5 central: 0.07 | range: 0.03–0.10
  DBDO central:  0.04 | range: 0.01–0.06

## 3. Degree of polymerization (DP) calibration

DP range determined by triangulation of three independent MW sources:

Source A — MALDI-TOF (this work):
  Provides DP_min (lower bound). MALDI underdetects high-MW oligomers
  due to ionization efficiency. Used as conservative floor estimate.
  DP_min = Mw_MALDI / M0_weighted

Source B — GPC literature (palm oil-specific):
  Mohamad Ibrahim et al. 2004: EFB soda Mw = 2444–3279 g/mol
  Hussin et al. 2013/2014: OPF kraft/soda/organosolv Mw trend confirmed
  (kraft > soda > organosolv regardless of biomass source)

Source C — GPC literature (equivalent non-wood biomass):
  Latif et al. 2021 (coconut husk): kraft 769, soda 959, org 606 g/mol
  Constant et al. 2016 (wheat straw soda P1000): 4000–7000 g/mol
  Constant et al. 2016 (Alcell hardwood org): 5000–8000 g/mol

M0 calculation per configuration:
  M0 = f_G × 180.2 + f_S × 210.2   (g/mol)

Final DP windows:
| Configuration | M0 (g/mol) | DP_min | DP_max |
|---------------|-----------|--------|--------|
| EFB Kraft     | 201.3     | 10     | 16     |
| EFB Soda      | 198.2     | 10     | 16     |
| EFB Org       | 188.8     |  3     |  8     |
| MF Kraft      | 192.4     | 10     | 17     |
| MF Soda       | 180.2     | 11     | 17     |
| MF Org        | 185.2     |  3     |  8     |
| KS Kraft      | 182.9     | 11     | 18     |
| KS Soda       | 180.2     | 11     | 17     |
| KS Org        | 182.9     |  3     |  8     |

## 4. van Krevelen analytical verification

Predicted H/C and O/C centroids from monomer composition (C9H10O3 per G
unit, C11H12O4 per S unit in backbone):

  H/C = (10·f_G + 12·f_S) / (9·f_G + 11·f_S)
  O/C = ( 3·f_G +  4·f_S) / (9·f_G + 11·f_S)

All predictions fall within the computational van Krevelen envelope
reported in Figure 7a of Oñate-Gutiérrez et al. 2025 (H/C: 1.096–1.111;
O/C: 0.333–0.356). Directional trends (EFB > KS in O/C; kraft/soda >
organosolv) are reproduced correctly.
