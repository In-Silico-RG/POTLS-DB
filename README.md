# POTLS-DB: Palm Oil Technical Lignin Structural Dataset

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19685101.svg)](https://doi.org/10.5281/zenodo.19685101)
[![License: CC BY 4.0](https://img.shields.io/badge/Data%20License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![License: MIT](https://img.shields.io/badge/Code%20License-MIT-blue.svg)](LICENSE_CODE)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![RDKit](https://img.shields.io/badge/RDKit-2024%2B-green.svg)](https://www.rdkit.org/)

**POTLS-DB** is a computational dataset of lignin molecular models derived from palm oil processing by-products — Empty Fruit Bunch (EFB), Mesocarp Fiber (MF), and Kernel Shell (KS) — extracted by three industrial pulping methods (kraft, soda, and organosolv).

The dataset was generated using the [Lignin Structure Generator (LGS)](https://github.com/sudhacheran/lignin-structure-generator) tool (Eswaran et al., *Sci. Data* 2022) with input parameters calibrated from experimental HSQC-NMR, MALDI-TOF, and GPC literature data. All generated structures are validated by a three-filter QC pipeline and characterized using RDKit molecular descriptors.

---

## Key numbers

| Item | Value |
|------|-------|
| Biomass sources | 3 (EFB, MF, KS) |
| Pulping methods | 3 (kraft, soda, organosolv) |
| LGS configurations | 27 (9 central + 9 low + 9 high sensitivity) |
| QC-passed structures | ~50,000–130,000 (version-dependent) |
| DP range | 3–18 (method- and source-specific) |
| Descriptors per structure | 20+ (MW, H/C, O/C, S/G, linkages, functional groups) |

---

## Repository structure

```
POTLS-DB/
├── configs/                    # LGS input parameter files (27 YAML)
│   ├── central/                # Central parameter estimates (9 configs)
│   ├── sensitivity_low/        # Lower-bound linkage variants (9 configs)
│   └── sensitivity_high/       # Upper-bound linkage variants (9 configs)
│
├── data/                       # Generated dataset (see Data Access below)
│   ├── structures/             # QC-passed + rejected structures (JSON)
│   └── convergence/            # Convergence statistics per config (CSV)
│
├── scripts/
│   ├── potls_batch_generator.py   # Main generation pipeline
│   └── potls_validation.py        # Multi-descriptor validation (Etapa 4)
│
├── notebooks/
│   ├── 01_parameter_calibration.ipynb   # MW triangulation + linkage matrix
│   ├── 02_qc_analysis.ipynb             # QC filter performance analysis
│   ├── 03_van_krevelen.ipynb            # van Krevelen diagrams
│   ├── 04_pearson_correlation.ipynb     # Descriptor correlation matrix
│   └── 05_chemical_space_tmap.ipynb     # TMAP visualization
│
├── docs/
│   ├── DATA_FORMAT.md          # JSON schema documentation
│   ├── PARAMETER_BASIS.md      # Rationale for each input parameter
│   └── QC_PROTOCOL.md          # QC filter specifications
│
├── .github/
│   └── workflows/
│       └── zenodo_release.yml  # Auto-deposit to Zenodo on GitHub release
│
├── README.md                   # This file
├── CITATION.cff                # Machine-readable citation
├── LICENSE                     # CC BY 4.0 (data)
├── LICENSE_CODE                # MIT (code)
├── requirements.txt            # Python dependencies
├── environment.yml             # Conda environment
└── .gitignore
```

---

## Data access

Large structure files (JSON) are **not stored in this repository** due to size. They are deposited on Zenodo:

> **Zenodo DOI:** [10.5281/zenodo.19685101](https://doi.org/10.5281/zenodo.19685101)

To download and set up locally:

```bash
git clone https://github.com/In-Silico-RG/POTLS-DB.git
cd POTLS-DB
pip install -r requirements.txt

# Download data from Zenodo (requires zenodo_get)
pip install zenodo_get
zenodo_get 10.5281/zenodo.19685101 -o data/
```

---

## Reproducing the dataset

### Prerequisites

- Python ≥ 3.10 with RDKit, pandas, numpy, tqdm, PyYAML
- Java ≥ 11 (for LGS)
- LGS jar: download from [Figshare DOI: 10.6084/m9.figshare.16915672.v2](https://doi.org/10.6084/m9.figshare.16915672.v2)

### Run central configurations

```bash
python scripts/potls_batch_generator.py \
    --configs  configs/central \
    --lgs-jar  /path/to/lgs.jar \
    --output   data/ \
    --set      central \
    --batch    500 \
    --max-str  10000 \
    --delta    0.01
```

### Run full sensitivity analysis (27 configurations)

```bash
python scripts/potls_batch_generator.py \
    --configs  configs \
    --lgs-jar  /path/to/lgs.jar \
    --output   data/ \
    --set      all \
    --batch    500 \
    --max-str  10000
```

### Dry-run (validate configs only)

```bash
python scripts/potls_batch_generator.py --configs configs/central --dry-run
```

---

## QC pipeline summary

Each generated structure passes three sequential filters:

| Filter | Criterion | Tolerance |
|--------|-----------|-----------|
| F1 — LGS-native | Valid SMILES + DP within configured range | — |
| F2 — Compositional | RDKit S/G within ±5% of target | 5% |
| F3 — MW calibration | MW within [DP_min × M₀ × 0.8, DP_max × M₀ × 1.2] | 20% |

Rejected structures are preserved in `*_rejected.json` files with per-structure QC annotations.

---

## Parameter basis

Input parameters for LGS were calibrated from three independent sources:

| Parameter | Source |
|-----------|--------|
| S/G ratios | HSQC-NMR experimental (Oñate-Gutiérrez et al. 2025) |
| Bond proportions | Eswaran et al. 2022 (hardwood/softwood ranges by archetype) |
| DP range (min) | MALDI-TOF Mw / M₀ (this work) |
| DP range (max) | GPC literature Mw P₇₅ / M₀ (Mohamad Ibrahim 2004; Hussin 2014; Constant 2016) |

See `docs/PARAMETER_BASIS.md` for full rationale and source table.

---

## Citation

If you use POTLS-DB in your research, please cite both the dataset and the parent paper:

**Dataset:**
```bibtex
@dataset{potls_db_2025,
  author    = {Combariza, Aldo F. and Oñate-Gutiérrez, Jesús A. and
               Vargas-Vergara, Sebastián D. and Combariza, Marianny Y. and
               Blanco-Tirado, Cristian and Pinzón, Julio R.},
  title     = {{POTLS-DB}: Palm Oil Technical Lignin Structural Dataset},
  year      = {2025},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.19685101},
  url       = {https://doi.org/10.5281/zenodo.19685101}
}
```

**Parent paper:**
```bibtex
@article{onate2025palm,
  author  = {Oñate-Gutiérrez, Jesús A. and Vargas-Vergara, Sebastián D. and
             Combariza, Aldo F. and Combariza, Marianny Y. and
             Blanco-Tirado, Cristian and Pinzón, Julio R.},
  title   = {Comparative study of palm oil lignins: Impact of biomass source
             and extraction method},
  journal = {[Journal TBD]},
  year    = {2025},
  doi     = {[DOI TBD]}
}
```

**LGS tool (required dependency):**
```bibtex
@article{eswaran2022lgs,
  author  = {Eswaran, Sudha cheranma devi and Subramaniam, Senthil and
             Sanyal, Udishnu and Rallo, Robert and Zhang, Xiao},
  title   = {Molecular structural dataset of lignin macromolecule
             elucidating experimental structural compositions},
  journal = {Scientific Data},
  volume  = {9},
  pages   = {647},
  year    = {2022},
  doi     = {10.1038/s41597-022-01709-4}
}
```

---

## License

- **Data** (all files in `data/`, `configs/`): [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE)
- **Code** (all `.py` files, notebooks, workflows): [MIT License](LICENSE_CODE)

---

## Acknowledgments

This work was supported by SGR-Gobernación de Santander Grant BPIN 2020000100251 and SGR-Minciencias BPIN 2021000100022. Computational resources provided by IN SILICO Research Group (Minciencias Category C), Universidad de Sucre.

---

## Contact

**IN SILICO Research Group**  
Facultad de Educación y Ciencias, Universidad de Sucre  
Sincelejo, Sucre, Colombia  
✉ insilico@unisucre.edu.co  
🔗 https://github.com/In-Silico-RG
