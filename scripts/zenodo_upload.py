"""
zenodo_upload.py
================
Helper script for depositing a new version of POTLS-DB to Zenodo.
Called by the GitHub Actions workflow on each release.

Can also be run manually for the initial deposit or large data uploads.

Usage (GitHub Actions - new version):
    python zenodo_upload.py \
        --token   $ZENODO_TOKEN \
        --concept $ZENODO_CONCEPT_DOI \
        --version v1.0.1 \
        --notes   "Patch notes..."

Usage (manual - initial deposit with data):
    python zenodo_upload.py \
        --token  $ZENODO_TOKEN \
        --init \
        --data-dir ./data \
        --version v1.0.0

Zenodo API reference:
    https://developers.zenodo.org/
"""

import argparse
import json
import os
import sys
import glob
import requests
from pathlib import Path

# ──────────────────────────────────────────────────────────────────────────────
# Zenodo API helpers
# ──────────────────────────────────────────────────────────────────────────────

BASE_URL = "https://zenodo.org/api"
SANDBOX_URL = "https://sandbox.zenodo.org/api"   # for testing


def get_session(token: str, sandbox: bool = False) -> requests.Session:
    s = requests.Session()
    s.headers.update({"Authorization": f"Bearer {token}"})
    return s, SANDBOX_URL if sandbox else BASE_URL


def create_new_version(session, base_url: str, concept_doi: str) -> dict:
    """Create a new version of an existing Zenodo record."""
    # Get concept record ID from DOI
    concept_id = concept_doi.split(".")[-1]
    r = session.post(f"{base_url}/deposit/depositions/{concept_id}/actions/newversion")
    r.raise_for_status()
    data = r.json()
    # New version draft URL
    draft_url = data["links"]["latest_draft"]
    r2 = session.get(draft_url)
    r2.raise_for_status()
    return r2.json()


def update_metadata(session, base_url: str, deposition_id: int,
                    metadata_file: str, version: str, notes: str) -> None:
    """Update deposition metadata from .zenodo.json."""
    with open(metadata_file) as f:
        meta = json.load(f)

    meta["version"] = version
    if notes:
        meta["notes"] = (meta.get("notes", "") + f"\n\nRelease notes: {notes}").strip()

    payload = {"metadata": meta}
    r = session.put(f"{base_url}/deposit/depositions/{deposition_id}",
                    json=payload)
    r.raise_for_status()
    print(f"  Metadata updated for deposition {deposition_id}")


def upload_file(session, bucket_url: str, filepath: Path) -> None:
    """Upload a single file to a Zenodo deposition bucket."""
    fname = filepath.name
    print(f"  Uploading {fname} ({filepath.stat().st_size / 1024:.1f} KB)...")
    with open(filepath, "rb") as f:
        r = session.put(f"{bucket_url}/{fname}", data=f)
    r.raise_for_status()
    print(f"    -> {r.status_code} OK")


def publish(session, base_url: str, deposition_id: int) -> str:
    """Publish (finalize) a deposition. Returns the assigned DOI."""
    r = session.post(f"{base_url}/deposit/depositions/{deposition_id}/actions/publish")
    r.raise_for_status()
    data = r.json()
    doi = data.get("doi", "")
    print(f"  Published! DOI: {doi}")
    return doi


# ──────────────────────────────────────────────────────────────────────────────
# File collectors
# ──────────────────────────────────────────────────────────────────────────────

def collect_code_files(version: str) -> list[Path]:
    """Collect code + config archive (created by GitHub Actions)."""
    pattern = f"POTLS_DB_code_{version}.zip"
    files = list(Path(".").glob(pattern))
    if not files:
        raise FileNotFoundError(f"Release archive not found: {pattern}")
    return files


def collect_data_files(data_dir: Path) -> list[Path]:
    """
    Collect large structure JSON files for initial/manual upload.
    Files are chunked into 500 MB archives to respect Zenodo file size limits.
    """
    jsons = sorted(data_dir.glob("structures/*.json"))
    csvs  = sorted(data_dir.glob("convergence/*.csv"))
    meta  = list(data_dir.glob("POTLS_DB_metadata.csv"))
    return jsons + csvs + meta


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description="POTLS-DB Zenodo uploader")
    p.add_argument("--token",    required=True,  help="Zenodo API token")
    p.add_argument("--concept",  default=None,   help="Zenodo concept DOI (for new version)")
    p.add_argument("--init",     action="store_true", help="Create initial deposit (no concept)")
    p.add_argument("--version",  required=True,  help="Version string, e.g. v1.0.0")
    p.add_argument("--notes",    default="",     help="Release notes string")
    p.add_argument("--data-dir", default=None,   type=Path, help="data/ dir for large files")
    p.add_argument("--sandbox",  action="store_true", help="Use Zenodo sandbox for testing")
    p.add_argument("--no-publish", action="store_true", help="Skip publish step (leave as draft)")
    return p.parse_args()


def main():
    args = parse_args()
    session, base_url = get_session(args.token, sandbox=args.sandbox)

    print(f"POTLS-DB Zenodo Upload — version {args.version}")
    print(f"Target: {'sandbox' if args.sandbox else 'PRODUCTION'} Zenodo")

    # ── Create or retrieve deposition ────────────────────────────────────────
    if args.init:
        # Initial deposit: create new blank deposition
        r = session.post(f"{base_url}/deposit/depositions", json={})
        r.raise_for_status()
        dep = r.json()
        print(f"  New deposition created: ID={dep['id']}")
    else:
        if not args.concept:
            print("ERROR: --concept required for new version upload")
            sys.exit(1)
        dep = create_new_version(session, base_url, args.concept)
        print(f"  New version draft created: ID={dep['id']}")

    dep_id    = dep["id"]
    bucket_url = dep["links"]["bucket"]

    # ── Update metadata ───────────────────────────────────────────────────────
    if Path(".zenodo.json").exists():
        update_metadata(session, base_url, dep_id, ".zenodo.json",
                        args.version, args.notes)
    else:
        print("  WARNING: .zenodo.json not found — metadata not updated")

    # ── Upload code archive ───────────────────────────────────────────────────
    print("  Uploading code + config archive...")
    try:
        for f in collect_code_files(args.version):
            upload_file(session, bucket_url, f)
    except FileNotFoundError as e:
        print(f"  WARNING: {e}")

    # ── Upload data files (manual/initial only) ───────────────────────────────
    if args.data_dir and args.data_dir.exists():
        print(f"  Uploading data files from {args.data_dir}...")
        for f in collect_data_files(args.data_dir):
            upload_file(session, bucket_url, f)

    # ── Publish ───────────────────────────────────────────────────────────────
    if not args.no_publish:
        doi = publish(session, base_url, dep_id)
        # Write DOI to file for GitHub Actions to read
        Path("zenodo_doi.txt").write_text(doi)
        print(f"\nSuccess! DOI: {doi}")
        print(f"URL: https://doi.org/{doi}")
    else:
        print(f"\nDraft saved (not published). ID: {dep_id}")
        print(f"Review at: https://zenodo.org/deposit/{dep_id}")


if __name__ == "__main__":
    main()
