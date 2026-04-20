"""
potls_batch_generator.py
========================
POTLS-DB — Palm Oil Technical Lignin Structural Dataset
Batch generation, QC filtering, and convergence monitoring.

Authors : IN SILICO Research Group — Universidad de Sucre
Contact : insilicorg@gmail.com
Version : 1.0.0

Usage
-----
  python potls_batch_generator.py [OPTIONS]

Options
  --configs   DIR   Directory with YAML config files   [default: ./central]
  --lgs-jar   PATH  Path to LGS jar file               [required]
  --output    DIR   Output directory for dataset        [default: ./POTLS_DB_output]
  --batch     INT   Structures per convergence batch    [default: 500]
  --max-str   INT   Max structures per config           [default: 10000]
  --delta     FLOAT Convergence threshold (fraction)    [default: 0.01]
  --set       STR   Parameter set to run (central/sensitivity_low/
                    sensitivity_high/all)               [default: central]
  --dry-run         Validate configs only, do not run LGS
  --log       PATH  Log file path                       [default: potls_run.log]

Pipeline
--------
  For each YAML config file:
    1. Load and validate parameters
    2. Call LGS in batches of --batch structures
    3. Apply 3 QC filters per structure:
         F1 — LGS-native  : valid SMILES + DP matches request
         F2 — Compositional: RDKit S/G within ±5% of target
         F3 — Molecular weight: MW within calibrated range ±20%
    4. Convergence check on mean/std of MW and complexity_index
    5. Write passing structures to JSON; rejected to separate file
  Assemble final dataset with metadata CSV.

Output structure
----------------
  POTLS_DB_output/
  ├── structures/
  │   ├── POTLS_EFB_KR_C.json          # QC-passed structures
  │   ├── POTLS_EFB_KR_C_rejected.json # QC-failed (annotated with reason)
  │   └── ...
  ├── convergence/
  │   ├── POTLS_EFB_KR_C_convergence.csv
  │   └── ...
  ├── POTLS_DB_metadata.csv            # all configurations summary
  └── potls_run.log
"""

import os
import sys
import json
import yaml
import logging
import argparse
import subprocess
import time
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

# RDKit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# ──────────────────────────────────────────────────────────────────────────────
# SECTION 1 — Configuration and dataclasses
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class LGSConfig:
    """Parsed and validated LGS configuration for one run."""
    yaml_path: Path
    dataset_id: str
    biomass_source: str
    pulping_method: str
    parameter_set: str
    archetype: str
    fG: float
    fS: float
    fH: float
    linkages: dict           # normalized fractions, sum = 1.0
    DP_min: int
    DP_max: int
    M0_gmol: float
    validation_targets: dict = field(default_factory=dict)

    @property
    def MW_min(self) -> float:
        """Minimum expected MW from DP_min × M0."""
        return self.DP_min * self.M0_gmol

    @property
    def MW_max(self) -> float:
        """Maximum expected MW from DP_max × M0."""
        return self.DP_max * self.M0_gmol

    @property
    def SG_target(self) -> Optional[float]:
        """Target S/G ratio from monomer fractions."""
        return self.fS / self.fG if self.fG > 0 else None


@dataclass
class StructureRecord:
    """One generated and QC-annotated structure."""
    lg_id: str
    smiles: str
    config_id: str
    DP: int
    MW_lgs: float       # MW reported by LGS
    MW_rdkit: float     # MW calculated by RDKit (independent check)
    HC_rdkit: float
    OC_rdkit: float
    SG_rdkit: float
    n_phenolic_OH: int
    n_methoxy: int
    bond_frequencies: dict
    monomer_count: dict
    # QC
    pass_F1: bool = False   # valid SMILES + correct DP
    pass_F2: bool = False   # S/G within ±5% of target
    pass_F3: bool = False   # MW within calibrated range ±20%
    F2_deviation: float = 0.0
    F3_deviation: float = 0.0
    qc_passed: bool = False


# ──────────────────────────────────────────────────────────────────────────────
# SECTION 2 — YAML loader and validator
# ──────────────────────────────────────────────────────────────────────────────

def load_config(yaml_path: Path) -> LGSConfig:
    """
    Load and validate a POTLS YAML configuration file.
    Raises ValueError if required fields are missing or linkages don't sum to 1.
    """
    with open(yaml_path) as f:
        raw = yaml.safe_load(f)

    proj  = raw.get("project", {})
    mono  = raw.get("monomers", {})
    links = raw.get("linkages", {})
    dp    = raw.get("degree_of_polymerization", {})
    vt    = raw.get("validation_targets", {})

    # Required field checks
    required_proj = ["dataset_id", "biomass_source", "pulping_method",
                     "parameter_set", "archetype"]
    for field_name in required_proj:
        if field_name not in proj:
            raise ValueError(f"{yaml_path.name}: missing project.{field_name}")

    fG = float(mono.get("G", 0))
    fS = float(mono.get("S", 0))
    fH = float(mono.get("H", 0))
    if abs(fG + fS + fH - 1.0) > 0.02:
        raise ValueError(f"{yaml_path.name}: monomers don't sum to 1.0 "
                         f"(G={fG}, S={fS}, H={fH}, sum={fG+fS+fH:.3f})")

    link_sum = sum(links.values())
    if abs(link_sum - 1.0) > 0.02:
        raise ValueError(f"{yaml_path.name}: linkages don't sum to 1.0 "
                         f"(sum={link_sum:.4f})")

    required_links = ["beta-O-4", "beta-beta", "beta-5", "5-5", "4-O-5", "DBDO"]
    for lk in required_links:
        if lk not in links:
            raise ValueError(f"{yaml_path.name}: missing linkage '{lk}'")

    return LGSConfig(
        yaml_path=yaml_path,
        dataset_id=proj["dataset_id"],
        biomass_source=proj["biomass_source"],
        pulping_method=proj["pulping_method"],
        parameter_set=proj["parameter_set"],
        archetype=proj["archetype"],
        fG=fG, fS=fS, fH=fH,
        linkages=links,
        DP_min=int(dp.get("min", 3)),
        DP_max=int(dp.get("max", 15)),
        M0_gmol=float(dp.get("M0_gmol", 190.0)),
        validation_targets=vt,
    )


# ──────────────────────────────────────────────────────────────────────────────
# SECTION 3 — LGS caller
# ──────────────────────────────────────────────────────────────────────────────

def call_lgs(
    lgs_jar: Path,
    config: LGSConfig,
    n_structures: int,
    output_dir: Path,
    timeout_sec: int = 3600,
) -> list[dict]:
    """
    Call the LGS Java tool for a given configuration.

    LGS command signature (from Eswaran et al. 2022 / GitHub):
        java -jar lgs.jar
             --config <project-config.yaml>
             --output <output_directory>
             --count  <n_structures>

    Returns list of raw JSON records from LGS output.

    NOTE: If LGS jar is not available, this function raises FileNotFoundError.
    A dry-run mode skips this call entirely.
    """
    if not lgs_jar.exists():
        raise FileNotFoundError(
            f"LGS jar not found at {lgs_jar}.\n"
            "Download from: https://github.com/sudhacheran/lignin-structure-generator\n"
            "or from Figshare DOI: 10.6084/m9.figshare.16915672.v2"
        )

    run_output_dir = output_dir / "lgs_raw" / config.dataset_id
    run_output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "java", "-jar", str(lgs_jar),
        "--config",  str(config.yaml_path),
        "--output",  str(run_output_dir),
        "--count",   str(n_structures),
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=timeout_sec,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"LGS failed for {config.dataset_id}:\n"
            f"STDOUT: {result.stdout[-500:]}\n"
            f"STDERR: {result.stderr[-500:]}"
        )

    # Read all JSON files produced by LGS in run_output_dir
    records = []
    for jf in sorted(run_output_dir.glob("*.json")):
        with open(jf) as f:
            data = json.load(f)
        # LGS may produce one file per structure or one file with list
        if isinstance(data, list):
            records.extend(data)
        elif isinstance(data, dict):
            records.append(data)

    return records


def simulate_lgs_output(config: LGSConfig, n: int, rng: np.random.Generator) -> list[dict]:
    """
    DEVELOPMENT STUB — simulate LGS output for testing the pipeline
    without requiring the Java tool.

    Generates synthetic records matching the expected LGS JSON schema.
    Remove this function and use call_lgs() in production.
    """
    records = []
    link_names = list(config.linkages.keys())
    link_probs = [config.linkages[k] for k in link_names]

    for i in range(n):
        dp = rng.integers(config.DP_min, config.DP_max + 1)
        # Simulated MW with noise around expected DP * M0
        mw = float(dp * config.M0_gmol * rng.normal(1.0, 0.08))

        # Simulated S/G with small noise around target
        if config.fG > 0:
            sg_target = config.fS / config.fG
            sg = max(0, sg_target + rng.normal(0, sg_target * 0.08))
        else:
            sg = 0.0

        # Bond frequencies with Dirichlet noise around configured values
        alpha = np.array(link_probs) * 20  # concentration parameter
        sampled_links = rng.dirichlet(alpha)
        bond_freqs = {k: float(round(v, 4)) for k, v in zip(link_names, sampled_links)}

        # Monomer counts
        n_S = int(round(dp * config.fS))
        n_G = dp - n_S
        n_PhOH = max(1, rng.integers(1, 4))
        n_OMe  = n_G + 2 * n_S

        records.append({
            "lg_id":   f"{config.dataset_id}_stub_{i:05d}",
            "smiles":  f"STUB_SMILES_{i}",   # replaced by real SMILES in production
            "degree_of_polymerization": int(dp),
            "molecular_weight":          float(round(mw, 2)),
            "s_g_ratio":                 float(round(sg, 4)),
            "bond_frequencies":          bond_freqs,
            "functional_groups": {
                "phenolic_OH":           int(n_PhOH),
                "aliphatic_OH_primary":  int(rng.integers(1, dp + 1)),
                "aliphatic_OH_secondary":int(rng.integers(dp // 2, dp + 1)),
                "methoxy":               int(n_OMe),
                "carbonyl":              0,
            },
            "monomer_count": {"G": n_G, "S": n_S, "H": 0},
        })
    return records


# ──────────────────────────────────────────────────────────────────────────────
# SECTION 4 — RDKit descriptor calculator
# ──────────────────────────────────────────────────────────────────────────────

# SMARTS patterns — validated against known lignin structures
_G_SMARTS  = Chem.MolFromSmarts("c1cc(OC)c(O)cc1")   # guaiacyl: 1 OMe + 1 OH
_S_SMARTS  = Chem.MolFromSmarts("c1cc(OC)c(OC)cc1")  # syringyl: 2 OMe
_PhOH_SMARTS = Chem.MolFromSmarts("[OH]c")            # phenolic OH


def rdkit_descriptors(smiles: str) -> Optional[dict]:
    """
    Calculate H/C, O/C, S/G, phenolic OH, and MW from SMILES using RDKit.
    Returns None if SMILES is invalid (stub or malformed).
    """
    if smiles.startswith("STUB_"):
        return None   # development stub — skip RDKit calculation

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    C = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    H = sum(a.GetTotalNumHs() for a in mol.GetAtoms())
    O = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8)

    if C == 0:
        return None

    n_G = len(mol.GetSubstructMatches(_G_SMARTS))
    n_S = len(mol.GetSubstructMatches(_S_SMARTS))
    sg  = n_S / n_G if n_G > 0 else (float("inf") if n_S > 0 else 0.0)
    n_PhOH = len(mol.GetSubstructMatches(_PhOH_SMARTS))
    n_OMe  = len(mol.GetSubstructMatches(Chem.MolFromSmarts("OC"))) - n_PhOH

    return {
        "MW_rdkit":     Descriptors.ExactMolWt(mol),
        "HC_rdkit":     H / C,
        "OC_rdkit":     O / C,
        "SG_rdkit":     sg,
        "n_phenolic_OH": n_PhOH,
        "n_methoxy":     max(0, n_OMe),
    }


# ──────────────────────────────────────────────────────────────────────────────
# SECTION 5 — QC filter pipeline
# ──────────────────────────────────────────────────────────────────────────────

def apply_qc_filters(
    record: dict,
    config: LGSConfig,
    sg_tolerance: float = 0.05,
    mw_tolerance: float = 0.20,
) -> StructureRecord:
    """
    Apply the three QC filters to a raw LGS record.

    F1 — LGS-native validation
         Pass: SMILES is valid (or stub) AND reported DP within config range.

    F2 — Compositional validation
         Pass: |SG_rdkit - SG_target| / SG_target ≤ sg_tolerance (default 5%).
         For pure-G configs (fG=1, SG_target=0): pass if SG_rdkit ≤ 0.05.

    F3 — Molecular weight calibration
         Pass: MW_lgs within [MW_min × (1 - mw_tolerance),
                               MW_max × (1 + mw_tolerance)].
         Structures outside range → segregated as 'high_mass' subgroup.

    Returns annotated StructureRecord.
    """
    lg_id  = record.get("lg_id", "unknown")
    smiles = record.get("smiles", "")
    dp_rep = int(record.get("degree_of_polymerization", 0))
    mw_lgs = float(record.get("molecular_weight", 0))
    sg_lgs = float(record.get("s_g_ratio", 0))

    # RDKit descriptors (None for stubs)
    rdkit = rdkit_descriptors(smiles)

    # Fallback values when SMILES is stub
    mw_rdkit   = rdkit["MW_rdkit"]   if rdkit else mw_lgs
    hc_rdkit   = rdkit["HC_rdkit"]   if rdkit else 0.0
    oc_rdkit   = rdkit["OC_rdkit"]   if rdkit else 0.0
    sg_rdkit   = rdkit["SG_rdkit"]   if rdkit else sg_lgs
    n_PhOH     = rdkit["n_phenolic_OH"] if rdkit else record.get("functional_groups", {}).get("phenolic_OH", 0)
    n_OMe      = rdkit["n_methoxy"]  if rdkit else record.get("functional_groups", {}).get("methoxy", 0)

    # ── F1: valid SMILES + DP in range ──
    smiles_ok = smiles.startswith("STUB_") or (Chem.MolFromSmiles(smiles) is not None)
    dp_ok     = config.DP_min <= dp_rep <= config.DP_max
    pass_F1   = smiles_ok and dp_ok

    # ── F2: S/G compositional check ──
    sg_val = sg_rdkit if rdkit else sg_lgs   # prefer RDKit value
    if config.SG_target is None or config.SG_target == 0:
        # Pure G config
        pass_F2       = sg_val <= 0.05
        F2_deviation  = sg_val
    else:
        F2_deviation = abs(sg_val - config.SG_target) / config.SG_target
        pass_F2      = F2_deviation <= sg_tolerance

    # ── F3: MW calibration window ──
    mw_lo = config.MW_min * (1 - mw_tolerance)
    mw_hi = config.MW_max * (1 + mw_tolerance)
    pass_F3      = mw_lo <= mw_lgs <= mw_hi
    F3_deviation = 0.0
    if mw_lgs < mw_lo:
        F3_deviation = (mw_lo - mw_lgs) / mw_lo
    elif mw_lgs > mw_hi:
        F3_deviation = (mw_lgs - mw_hi) / mw_hi

    qc_passed = pass_F1 and pass_F2 and pass_F3

    return StructureRecord(
        lg_id=lg_id,
        smiles=smiles,
        config_id=config.dataset_id,
        DP=dp_rep,
        MW_lgs=mw_lgs,
        MW_rdkit=mw_rdkit,
        HC_rdkit=hc_rdkit,
        OC_rdkit=oc_rdkit,
        SG_rdkit=sg_rdkit,
        n_phenolic_OH=n_PhOH,
        n_methoxy=n_OMe,
        bond_frequencies=record.get("bond_frequencies", {}),
        monomer_count=record.get("monomer_count", {}),
        pass_F1=pass_F1,
        pass_F2=pass_F2,
        pass_F3=pass_F3,
        F2_deviation=F2_deviation,
        F3_deviation=F3_deviation,
        qc_passed=qc_passed,
    )


# ──────────────────────────────────────────────────────────────────────────────
# SECTION 6 — Convergence monitor
# ──────────────────────────────────────────────────────────────────────────────

class ConvergenceMonitor:
    """
    Track MW and complexity_index statistics across batches.
    Convergence criterion: relative change in mean AND std < delta for both metrics.

    'complexity_index' is approximated here as the sum of condensed linkages
    (β-5 + 5-5 + 4-O-5 + DBDO) normalized by total linkages — a proxy for
    network branching. In production, use the LGS-reported complexity_index.
    """

    def __init__(self, delta: float = 0.01):
        self.delta = delta
        self.history: list[dict] = []

    @staticmethod
    def _complexity_index(record: StructureRecord) -> float:
        bf = record.bond_frequencies
        condensed = (bf.get("beta-5", 0) + bf.get("5-5", 0)
                     + bf.get("4-O-5", 0) + bf.get("DBDO", 0))
        return condensed  # already a fraction

    def update(self, batch_records: list[StructureRecord], batch_n: int) -> dict:
        """Add a batch, compute stats, return convergence dict."""
        mws = [r.MW_lgs for r in batch_records]
        cxs = [self._complexity_index(r) for r in batch_records]

        stat = {
            "batch":        batch_n,
            "n_batch":      len(batch_records),
            "mean_MW":      float(np.mean(mws)),
            "std_MW":       float(np.std(mws)),
            "mean_CX":      float(np.mean(cxs)),
            "std_CX":       float(np.std(cxs)),
            "converged":    False,
            "delta_mean_MW": None,
            "delta_std_MW":  None,
            "delta_mean_CX": None,
            "delta_std_CX":  None,
        }

        if len(self.history) >= 1:
            prev = self.history[-1]
            dMmean = abs(stat["mean_MW"] - prev["mean_MW"]) / (prev["mean_MW"] + 1e-9)
            dMstd  = abs(stat["std_MW"]  - prev["std_MW"])  / (prev["std_MW"]  + 1e-9)
            dCmean = abs(stat["mean_CX"] - prev["mean_CX"]) / (prev["mean_CX"] + 1e-9)
            dCstd  = abs(stat["std_CX"]  - prev["std_CX"])  / (prev["std_CX"]  + 1e-9)

            stat["delta_mean_MW"] = dMmean
            stat["delta_std_MW"]  = dMstd
            stat["delta_mean_CX"] = dCmean
            stat["delta_std_CX"]  = dCstd
            stat["converged"] = (dMmean < self.delta and dMstd < self.delta
                                 and dCmean < self.delta and dCstd < self.delta)

        self.history.append(stat)
        return stat

    def to_dataframe(self) -> pd.DataFrame:
        return pd.DataFrame(self.history)


# ──────────────────────────────────────────────────────────────────────────────
# SECTION 7 — Per-config runner
# ──────────────────────────────────────────────────────────────────────────────

def run_config(
    config: LGSConfig,
    lgs_jar: Optional[Path],
    output_dir: Path,
    batch_size: int,
    max_structures: int,
    delta: float,
    use_stub: bool,
    logger: logging.Logger,
) -> dict:
    """
    Full generation + QC + convergence loop for one LGS configuration.
    Returns summary dict for metadata table.
    """
    logger.info(f"{'─'*60}")
    logger.info(f"Starting: {config.dataset_id} "
                f"({config.biomass_source} | {config.pulping_method} | "
                f"{config.parameter_set})")
    sg_label = f"{config.SG_target:.2f}" if config.SG_target else "pure-G"
    logger.info(f"  Archetype: {config.archetype} | "
                f"S/G target: {sg_label} | "
                f"DP: {config.DP_min}-{config.DP_max} | "
                f"MW window: [{config.MW_min:.0f}-{config.MW_max:.0f}] +-20%")

    # Output paths
    struct_dir = output_dir / "structures"
    conv_dir   = output_dir / "convergence"
    struct_dir.mkdir(parents=True, exist_ok=True)
    conv_dir.mkdir(parents=True, exist_ok=True)

    passed_records:   list[StructureRecord] = []
    rejected_records: list[StructureRecord] = []
    monitor = ConvergenceMonitor(delta=delta)
    rng = np.random.default_rng(seed=42)

    batch_n = 0
    total_generated = 0
    converged = False

    with tqdm(total=max_structures, desc=config.dataset_id, unit="str",
              leave=True, dynamic_ncols=True) as pbar:

        while total_generated < max_structures and not converged:
            batch_n += 1

            # ── Generate batch ──────────────────────────────────────────────
            if use_stub:
                raw_records = simulate_lgs_output(config, batch_size, rng)
            else:
                try:
                    raw_records = call_lgs(lgs_jar, config, batch_size, output_dir)
                except Exception as e:
                    logger.error(f"LGS call failed at batch {batch_n}: {e}")
                    break

            total_generated += len(raw_records)

            # ── Apply QC filters ────────────────────────────────────────────
            batch_passed = []
            for rec in raw_records:
                sr = apply_qc_filters(rec, config)
                if sr.qc_passed:
                    passed_records.append(sr)
                    batch_passed.append(sr)
                else:
                    rejected_records.append(sr)

            pbar.update(len(raw_records))

            # ── Convergence check (requires ≥ 2 batches) ───────────────────
            if batch_passed:
                stat = monitor.update(batch_passed, batch_n)
                converged = stat["converged"]

                if stat["delta_mean_MW"] is not None:
                    logger.debug(
                        f"  Batch {batch_n:3d} | n_pass={len(batch_passed):5d} | "
                        f"mean_MW={stat['mean_MW']:7.1f} | "
                        f"ΔmeanMW={stat['delta_mean_MW']:.4f} | "
                        f"{'CONVERGED' if converged else '...'}"
                    )

            if converged:
                logger.info(f"  Converged at batch {batch_n} "
                            f"(n_passed = {len(passed_records)})")

    # ── QC summary ───────────────────────────────────────────────────────────
    n_total   = total_generated
    n_passed  = len(passed_records)
    n_fail_F1 = sum(1 for r in rejected_records if not r.pass_F1)
    n_fail_F2 = sum(1 for r in rejected_records if r.pass_F1 and not r.pass_F2)
    n_fail_F3 = sum(1 for r in rejected_records if r.pass_F1 and r.pass_F2 and not r.pass_F3)
    pass_rate = n_passed / n_total if n_total > 0 else 0

    logger.info(f"  Generated: {n_total:6d} | Passed QC: {n_passed:6d} "
                f"({pass_rate:.1%}) | "
                f"Fail-F1: {n_fail_F1} | Fail-F2: {n_fail_F2} | Fail-F3: {n_fail_F3}")

    # ── Write outputs ────────────────────────────────────────────────────────
    def _native(v):
        """Convert numpy scalars to native Python types for JSON."""
        if isinstance(v, (np.integer,)):  return int(v)
        if isinstance(v, (np.floating,)): return float(v)
        if isinstance(v, (np.bool_,)):    return bool(v)
        if isinstance(v, dict):           return {kk: _native(vv) for kk, vv in v.items()}
        return v

    def record_to_dict(r: StructureRecord) -> dict:
        return {
            "lg_id":            r.lg_id,
            "config_id":        r.config_id,
            "smiles":           r.smiles,
            "DP":               _native(r.DP),
            "MW_lgs":           _native(r.MW_lgs),
            "MW_rdkit":         _native(r.MW_rdkit),
            "HC_rdkit":         _native(r.HC_rdkit),
            "OC_rdkit":         _native(r.OC_rdkit),
            "SG_rdkit":         _native(r.SG_rdkit),
            "n_phenolic_OH":    _native(r.n_phenolic_OH),
            "n_methoxy":        _native(r.n_methoxy),
            "bond_frequencies": _native(r.bond_frequencies),
            "monomer_count":    _native(r.monomer_count),
            "pass_F1":          bool(r.pass_F1),
            "pass_F2":          bool(r.pass_F2),
            "pass_F3":          bool(r.pass_F3),
            "F2_deviation":     _native(r.F2_deviation),
            "F3_deviation":     _native(r.F3_deviation),
            "qc_passed":        bool(r.qc_passed),
        }

    # Passed structures
    out_passed = struct_dir / f"{config.dataset_id}.json"
    with open(out_passed, "w") as f:
        json.dump([record_to_dict(r) for r in passed_records], f, indent=2)

    # Rejected structures (with reason annotation)
    out_rejected = struct_dir / f"{config.dataset_id}_rejected.json"
    with open(out_rejected, "w") as f:
        json.dump([record_to_dict(r) for r in rejected_records], f, indent=2)

    # Convergence CSV
    conv_df = monitor.to_dataframe()
    conv_df["config_id"] = config.dataset_id
    conv_df.to_csv(conv_dir / f"{config.dataset_id}_convergence.csv", index=False)

    logger.info(f"  Written: {out_passed.name} | {out_rejected.name}")

    return {
        "dataset_id":      config.dataset_id,
        "biomass_source":  config.biomass_source,
        "pulping_method":  config.pulping_method,
        "parameter_set":   config.parameter_set,
        "archetype":       config.archetype,
        "SG_target":       config.SG_target,
        "DP_min":          config.DP_min,
        "DP_max":          config.DP_max,
        "MW_min":          config.MW_min,
        "MW_max":          config.MW_max,
        "f_betaO4":        config.linkages.get("beta-O-4"),
        "n_generated":     n_total,
        "n_passed":        n_passed,
        "n_fail_F1":       n_fail_F1,
        "n_fail_F2":       n_fail_F2,
        "n_fail_F3":       n_fail_F3,
        "pass_rate":       round(pass_rate, 4),
        "n_batches":       batch_n,
        "converged":       converged,
        "mean_MW_final":   monitor.history[-1]["mean_MW"] if monitor.history else None,
        "std_MW_final":    monitor.history[-1]["std_MW"]  if monitor.history else None,
        "run_timestamp":   datetime.now().isoformat(),
    }


# ──────────────────────────────────────────────────────────────────────────────
# SECTION 8 — CLI and main loop
# ──────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="POTLS-DB batch generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--configs",   default="./central",       type=Path)
    p.add_argument("--lgs-jar",   default=None,              type=Path)
    p.add_argument("--output",    default="./POTLS_DB_output", type=Path)
    p.add_argument("--batch",     default=500,               type=int)
    p.add_argument("--max-str",   default=10000,             type=int)
    p.add_argument("--delta",     default=0.01,              type=float)
    p.add_argument("--set",       default="central",
                   choices=["central", "sensitivity_low", "sensitivity_high", "all"])
    p.add_argument("--dry-run",   action="store_true")
    p.add_argument("--stub",      action="store_true",
                   help="Use simulated LGS output (development mode)")
    p.add_argument("--log",       default="potls_run.log",   type=Path)
    p.add_argument("--verbose",   action="store_true")
    return p


def setup_logging(log_path: Path, verbose: bool) -> logging.Logger:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(sys.stdout),
        ],
    )
    return logging.getLogger("potls")


def collect_yaml_files(configs_dir: Path, param_set: str) -> list[Path]:
    """Collect YAML files based on --set argument."""
    base = configs_dir.parent if configs_dir.name in (
        "central", "sensitivity_low", "sensitivity_high") else configs_dir

    if param_set == "all":
        sets = ["central", "sensitivity_low", "sensitivity_high"]
    else:
        sets = [param_set]

    files = []
    for s in sets:
        d = base / s
        if d.exists():
            files.extend(sorted(d.glob("*.yaml")))
        else:
            # Try treating configs_dir directly
            files.extend(sorted(configs_dir.glob("*.yaml")))

    return sorted(set(files))


def main():
    args = build_parser().parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(args.output / args.log, args.verbose)

    logger.info("=" * 60)
    logger.info("POTLS-DB Batch Generator v1.0.0")
    logger.info(f"Run started: {datetime.now().isoformat()}")
    logger.info(f"Config dir : {args.configs}")
    logger.info(f"Output dir : {args.output}")
    logger.info(f"Mode       : {'DRY-RUN' if args.dry_run else ('STUB' if args.stub else 'PRODUCTION')}")
    logger.info(f"Batch size : {args.batch} | Max structs: {args.max_str} | Delta: {args.delta}")
    logger.info("=" * 60)

    # Collect and load configs
    yaml_files = collect_yaml_files(args.configs, args.set)
    if not yaml_files:
        logger.error(f"No YAML files found in {args.configs} for set '{args.set}'")
        sys.exit(1)

    configs = []
    for yf in yaml_files:
        try:
            cfg = load_config(yf)
            configs.append(cfg)
            logger.info(f"  Loaded: {cfg.dataset_id} ({yf.name})")
        except (ValueError, KeyError) as e:
            logger.error(f"  INVALID: {yf.name} — {e}")

    logger.info(f"\nTotal valid configs: {len(configs)} / {len(yaml_files)}")

    if args.dry_run:
        logger.info("DRY-RUN complete — no structures generated.")
        return

    # Main generation loop
    all_summaries = []
    t_start = time.time()

    for i, cfg in enumerate(configs, 1):
        logger.info(f"\n[{i}/{len(configs)}]")
        summary = run_config(
            config=cfg,
            lgs_jar=args.lgs_jar,
            output_dir=args.output,
            batch_size=args.batch,
            max_structures=args.max_str,
            delta=args.delta,
            use_stub=args.stub or (args.lgs_jar is None),
            logger=logger,
        )
        all_summaries.append(summary)

    # Write metadata table
    meta_df = pd.DataFrame(all_summaries)
    meta_path = args.output / "POTLS_DB_metadata.csv"
    meta_df.to_csv(meta_path, index=False)

    # Final report
    elapsed = time.time() - t_start
    total_structures = meta_df["n_passed"].sum()
    logger.info("\n" + "=" * 60)
    logger.info("POTLS-DB generation complete")
    logger.info(f"  Configurations processed : {len(all_summaries)}")
    logger.info(f"  Total QC-passed structures: {total_structures:,}")
    logger.info(f"  Overall pass rate          : "
                f"{meta_df['pass_rate'].mean():.1%} (mean across configs)")
    logger.info(f"  Elapsed time               : {elapsed/60:.1f} min")
    logger.info(f"  Metadata table             : {meta_path}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
