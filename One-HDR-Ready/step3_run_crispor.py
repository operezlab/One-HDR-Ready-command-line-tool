#!/usr/bin/env python3
"""
Step 3 — Run CRISPOR and Select Best sgRNA (automatic)
------------------------------------------------------

This version is patched to:

✔ Always use the local venv python
✔ Use genome directory from Makefile → crisporGenomes/hg38/
✔ Fail gracefully with clear error messages
✔ Fully portable to any system

Outputs:
    results/step3_crispor_output.tsv
    results/step3_best_sgRNA.json
"""

from __future__ import annotations
import json
import sys
import subprocess
from pathlib import Path
from typing import Dict, Any, Optional

import pandas as pd
from Bio.Seq import Seq


# ──────────────────────────────────────────────────────────────
# PATHS
# ──────────────────────────────────────────────────────────────

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

SCRIPT_DIR = Path(__file__).resolve().parent

# Project root = OneHDRSingleInstall/
PROJECT_ROOT = SCRIPT_DIR.parents[1]

# CRISPOR location
CRISPOR_SCRIPT = PROJECT_ROOT / "crisporWebsite" / "crispor.py"

# Genomes installed by Makefile:
GENOME_DIR = PROJECT_ROOT / "crisporGenomes"

# Step 1 input files
INPUT_FASTA = RESULTS_DIR / "step1_crispor_input.fa"
SEQS_JSON   = RESULTS_DIR / "step1_sequences.json"
CRISPOR_TSV = RESULTS_DIR / "step3_crispor_output.tsv"


# ──────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────

def get_local_python39(project_root: Path) -> str:
    """Return the python interpreter inside local venv (Linux/macOS/WSL + Windows)."""
    unix = project_root / "venv" / "bin" / "python"
    win  = project_root / "venv" / "Scripts" / "python.exe"

    if unix.exists():
        return str(unix)
    if win.exists():
        return str(win)

    raise RuntimeError(
        f"ERROR: Could not find Python inside venv/\n"
        f"Expected at:\n  {unix}\n  {win}\n"
        f"Did you run:  make onehdr  ?"
    )

PY39 = get_local_python39(PROJECT_ROOT)



def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def _load_sequences() -> Dict[str, Any]:
    if not SEQS_JSON.exists():
        raise FileNotFoundError(f"Missing Step 1 JSON: {SEQS_JSON}")
    return json.loads(SEQS_JSON.read_text())


def _find_edge_index(seqs: Dict[str, Any]) -> Optional[int]:
    """Find the exon edge (gene-sense)."""
    exon = seqs["seq_exon"]
    tmpl = seqs["seq_primer_template"]
    reverse = bool(seqs.get("reverseStrand", False))

    if not exon or not tmpl:
        return None

    # use 15-mer from correct side
    edge15 = exon[:15] if reverse else exon[-15:]
    idx = tmpl.find(edge15)

    if idx != -1:
        return idx if reverse else idx + 14

    # fallback shorter
    for k in (12, 10, 8):
        ek = exon[:k] if reverse else exon[-k:]
        idx = tmpl.find(ek)
        if idx != -1:
            return idx if reverse else idx + k - 1
    return None


def _guide_position_in_template(guide: str, orientation: str, tmpl: str) -> Optional[int]:
    query = guide if orientation == "forw" else reverse_complement(guide)
    pos = tmpl.find(query)
    return None if pos == -1 else pos


def _distance_from_edge(edge_idx: int, guide_start: int) -> int:
    return guide_start - edge_idx



# ──────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────

def run() -> Dict[str, Any]:

    # 0. Confirm genome installation
    if not GENOME_DIR.exists():
        raise RuntimeError(
            f"CRISPOR genome NOT FOUND:\n"
            f"  {GENOME_DIR}\n\n"
            f"Run:  make install-genomes\n"
        )

    # 1. Check CRISPOR script
    if not CRISPOR_SCRIPT.exists():
        raise FileNotFoundError(
            f"CRISPOR script not found:\n  {CRISPOR_SCRIPT}\n"
            f"Expected inside OneHDRSingleInstall/crisporWebsite/"
        )

    # 2. Run CRISPOR
    if not INPUT_FASTA.exists():
        raise FileNotFoundError(f"Missing CRISPOR input FASTA: {INPUT_FASTA}")

    print("\nUsing CRISPOR script:", CRISPOR_SCRIPT, "\n")
    print("Running CRISPOR:")
    cmd = [
        PY39,
        str(CRISPOR_SCRIPT),
        "hg38",
        str(INPUT_FASTA),
        str(CRISPOR_TSV),
        "--genomeDir", str(GENOME_DIR)
    ]

    print(" ", " ".join(str(x) for x in cmd), "\n")

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if proc.returncode != 0:
        raise RuntimeError(
            "[Step 3] ERROR: CRISPOR FAILED (exit {code})\n"
            "--- STDOUT ---\n{out}\n\n"
            "--- STDERR ---\n{err}".format(
                code=proc.returncode,
                out=proc.stdout,
                err=proc.stderr
            )
        )

    print("CRISPOR run complete.\n")

    if not CRISPOR_TSV.exists():
        raise FileNotFoundError("CRISPOR output TSV missing after run.")


    # 3. Load sequence context
    seqs = _load_sequences()
    primer_template = seqs.get("seq_primer_template", "")
    edge_idx = _find_edge_index(seqs)


    # 4. Parse CRISPOR TSV
    df = pd.read_csv(CRISPOR_TSV, sep="\t")

    needed = [
        "guideId", "targetSeq", "mitSpecScore", "cfdSpecScore", "Doench '16-Score",
        "Moreno-Mateos-Score", "Doench-RuleSet3-Score", "Out-of-Frame-Score",
        "Lindel-Score"
    ]
    for c in needed:
        if c not in df.columns:
            df[c] = pd.NA


    # Extract orientation and number from guideId
    loc_orient = df["guideId"].astype(str).str.extract(r"(\d+)(forw|rev)?")
    df["location"] = pd.to_numeric(loc_orient[0], errors="coerce")
    df["orientation"] = loc_orient[1].fillna("forw")

    # Distance from exon
    distances = []
    for _, row in df.iterrows():
        guide = str(row["targetSeq"])
        orient = row["orientation"]
        if not guide or not primer_template or edge_idx is None:
            distances.append(pd.NA)
            continue
        gpos = _guide_position_in_template(guide, orient, primer_template)
        distances.append(pd.NA if gpos is None else _distance_from_edge(edge_idx, gpos))

    df["DistanceFromExon"] = distances


    # 5. Rank guides
    df_sorted = df.copy()
    df_sorted["absD"] = df_sorted["DistanceFromExon"].apply(lambda x: abs(x) if pd.notna(x) else 1e9)
    df_sorted["Doench16_num"] = pd.to_numeric(df_sorted["Doench '16-Score"], errors="coerce")
    df_sorted = df_sorted.sort_values(by=["absD", "Doench16_num"], ascending=[True, False])

    best_row = df_sorted.iloc[0].to_dict()

    best = {
        "targetSeq": best_row.get("targetSeq"),
        "orientation": best_row.get("orientation"),
        "DistanceFromExon": best_row.get("DistanceFromExon"),
        "scores": {
            "mitSpecScore": best_row.get("mitSpecScore"),
            "cfdSpecScore": best_row.get("cfdSpecScore"),
            "Doench16": best_row.get("Doench '16-Score"),
            "MorenoMateos": best_row.get("Moreno-Mateos-Score"),
            "RuleSet3": best_row.get("Doench-RuleSet3-Score"),
            "OutOfFrame": best_row.get("Out-of-Frame-Score"),
            "Lindel": best_row.get("Lindel-Score"),
        }
    }

    (RESULTS_DIR / "step3_best_sgRNA.json").write_text(json.dumps(best, indent=2))
    print("Best sgRNA saved → results/step3_best_sgRNA.json\n")

    return {"best_guide": best}



if __name__ == "__main__":
    try:
        out = run()
        b = out["best_guide"]
        print("[Step 3 OK]")
        print(f"Best sgRNA: {b['targetSeq']}  | orientation={b['orientation']}  | distance={b['DistanceFromExon']}")
    except Exception as e:
        print(e)
        sys.exit(1)
