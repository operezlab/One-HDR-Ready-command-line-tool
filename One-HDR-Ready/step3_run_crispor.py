#!/usr/bin/env python3
"""
Step 3 — Run CRISPOR (NO efficiency models) and Select Best sgRNA
----------------------------------------------------------------

This version:
✔ Bulletproof project-root detection (crisporWebsite + crisporGenomes)
✔ Uses sys.executable
✔ Runs CRISPOR with --noEffScores (portable; avoids rs3/Azimuth issues)
✔ Parses CRISPOR TSV and ranks guides using off-target specificity + distance
✔ Ignores mismatch 0–3 off-target count fields entirely
✔ Keeps ONLY mismatch-4 off-target count as an optional tie-breaker (if present)

Outputs:
    results/step3_crispor_output.tsv
    results/step3_best_sgRNA.json
"""

from __future__ import annotations

import json
import sys
import subprocess
import re
from pathlib import Path
from typing import Dict, Any, Optional, List

import pandas as pd
from Bio.Seq import Seq

def _json_safe(x):
    """Convert pandas/numpy missing values to None for JSON."""
    try:
        # handles pd.NA, np.nan, NaT, etc.
        if pd.isna(x):
            return None
    except Exception:
        pass
    # convert numpy scalars to Python scalars
    if hasattr(x, "item"):
        try:
            return x.item()
        except Exception:
            pass
    return x


# ──────────────────────────────────────────────────────────────
# PATH DETECTION
# ──────────────────────────────────────────────────────────────

SCRIPT_DIR = Path(__file__).resolve().parent


def find_project_root(start: Path) -> Path:
    cur = start
    for _ in range(8):
        if (cur / "crisporWebsite").exists() and (cur / "crisporGenomes").exists():
            return cur
        cur = cur.parent

    raise RuntimeError(
        "ERROR: Could not locate project root.\n"
        "Expected to find:\n"
        "  crisporWebsite/\n"
        "  crisporGenomes/\n"
        f"Starting search from: {start}"
    )


PROJECT_ROOT = find_project_root(SCRIPT_DIR)

RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

CRISPOR_SCRIPT = PROJECT_ROOT / "crisporWebsite" / "crispor.py"
GENOME_DIR = PROJECT_ROOT / "crisporGenomes"

INPUT_FASTA = RESULTS_DIR / "step1_crispor_input.fa"
SEQS_JSON = RESULTS_DIR / "step1_sequences.json"
CRISPOR_TSV = RESULTS_DIR / "step3_crispor_output.tsv"

PYTHON = sys.executable


# ──────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────

def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def _load_sequences() -> Dict[str, Any]:
    if not SEQS_JSON.exists():
        raise FileNotFoundError(f"Missing Step 1 JSON: {SEQS_JSON}")
    return json.loads(SEQS_JSON.read_text())


def _find_edge_index(seqs: Dict[str, Any]) -> Optional[int]:
    exon = seqs.get("seq_exon", "")
    tmpl = seqs.get("seq_primer_template", "")
    reverse = bool(seqs.get("reverseStrand", False))

    if not exon or not tmpl:
        return None

    edge15 = exon[:15] if reverse else exon[-15:]
    idx = tmpl.find(edge15)
    if idx != -1:
        return idx if reverse else idx + 14

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


def _ensure_cols(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    for c in cols:
        if c not in df.columns:
            df[c] = pd.NA
    return df


def _find_mm4_column(df: pd.DataFrame) -> Optional[str]:
    """
    Try to find a mismatch-4 off-target count column.
    We intentionally ignore mm0–mm3 columns.
    """
    candidates = [
        "otCountMM4",
        "otCountMm4",
        "offTargetCountMM4",
        "offTargetsMM4",
        "otCount_4MM",
        "otCount_4",
    ]
    for c in candidates:
        if c in df.columns:
            return c

    rgx_list = [
        re.compile(r"^otCount.*MM[_-]?4$", re.IGNORECASE),
        re.compile(r".*MM[_-]?4.*Count.*", re.IGNORECASE),
        re.compile(r".*off.*MM[_-]?4.*", re.IGNORECASE),
    ]
    for col in df.columns:
        for rgx in rgx_list:
            if rgx.match(str(col)):
                return col

    return None


# ──────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────

def run() -> Dict[str, Any]:
    if not GENOME_DIR.exists():
        raise RuntimeError(
            f"CRISPOR genome NOT FOUND:\n  {GENOME_DIR}\n\n"
            f"Run:  make install-genomes"
        )

    if not CRISPOR_SCRIPT.exists():
        raise FileNotFoundError(
            f"CRISPOR script not found:\n  {CRISPOR_SCRIPT}\n"
            f"Expected inside: crisporWebsite/"
        )

    if not INPUT_FASTA.exists():
        raise FileNotFoundError(f"Missing CRISPOR input FASTA: {INPUT_FASTA}")

    print("\nUsing CRISPOR script:", CRISPOR_SCRIPT)
    print("Genome directory:", GENOME_DIR)
    print("\nRunning CRISPOR (no efficiency models)...\n")

    cmd = [
        PYTHON,
        str(CRISPOR_SCRIPT),
        "hg38",
        str(INPUT_FASTA),
        str(CRISPOR_TSV),
        "--genomeDir", str(GENOME_DIR),
        "--noEffScores",
    ]

    print("Command:\n ", " ".join(str(x) for x in cmd), "\n")

    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    if proc.returncode != 0:
        raise RuntimeError(
            f"[Step 3] ERROR: CRISPOR FAILED (exit {proc.returncode})\n"
            f"--- STDOUT ---\n{proc.stdout}\n\n"
            f"--- STDERR ---\n{proc.stderr}"
        )

    if not CRISPOR_TSV.exists():
        raise FileNotFoundError("CRISPOR output TSV missing after run.")

    print("CRISPOR run complete.\n")

    seqs = _load_sequences()
    primer_template = seqs.get("seq_primer_template", "")
    edge_idx = _find_edge_index(seqs)

    df = pd.read_csv(CRISPOR_TSV, sep="\t")

    # We only require these for ranking + output
    df = _ensure_cols(df, ["guideId", "targetSeq", "mitSpecScore", "cfdSpecScore"])

    # mm4 count only (optional)
    mm4_col = _find_mm4_column(df)
    if mm4_col is None:
        df["otCountMM4_only"] = pd.NA
        mm4_col = "otCountMM4_only"

    # Parse orientation + location (guideId often like "123forw"/"123rev")
    loc_orient = df["guideId"].astype(str).str.extract(r"(\d+)(forw|rev)?")
    df["location"] = pd.to_numeric(loc_orient[0], errors="coerce")
    df["orientation"] = loc_orient[1].fillna("forw")

    # DistanceFromExon
    distances = []
    for _, row in df.iterrows():
        guide = str(row["targetSeq"]) if pd.notna(row["targetSeq"]) else ""
        orient = str(row["orientation"]) if pd.notna(row["orientation"]) else "forw"

        if not guide or not primer_template or edge_idx is None:
            distances.append(pd.NA)
            continue

        gpos = _guide_position_in_template(guide, orient, primer_template)
        distances.append(pd.NA if gpos is None else _distance_from_edge(edge_idx, gpos))

    df["DistanceFromExon"] = distances

    # Rank metrics
    df["absD"] = df["DistanceFromExon"].apply(lambda x: abs(x) if pd.notna(x) else 1e9)
    df["cfd_num"] = pd.to_numeric(df["cfdSpecScore"], errors="coerce")
    df["mit_num"] = pd.to_numeric(df["mitSpecScore"], errors="coerce")

    # mm4 off-target count only (lower is better). Unknown -> worst.
    df["mm4_num"] = pd.to_numeric(df[mm4_col], errors="coerce")
    df["mm4_rank"] = df["mm4_num"].apply(lambda x: x if pd.notna(x) else 1e9)

    # Sort order: distance, CFD, MIT, mm4 count
    df_sorted = df.sort_values(
        by=["absD", "cfd_num", "mit_num", "mm4_rank"],
        ascending=[True, False, False, True],
        na_position="last",
    )

    best_row = df_sorted.iloc[0].to_dict()

    best = {
        "targetSeq": _json_safe(best_row.get("targetSeq")),
        "orientation": _json_safe(best_row.get("orientation")),
        "DistanceFromExon": _json_safe(best_row.get("DistanceFromExon")),
        "scores": {
            "cfdSpecScore": _json_safe(best_row.get("cfdSpecScore")),
            "mitSpecScore": _json_safe(best_row.get("mitSpecScore")),
            "offtargets_mm4_count": _json_safe(best_row.get(mm4_col)),
        },
        "meta": {
            "guideId": _json_safe(best_row.get("guideId")),
            "location": _json_safe(best_row.get("location")),
        },
    }


    (RESULTS_DIR / "step3_best_sgRNA.json").write_text(json.dumps(best, indent=2))
    print("Best sgRNA saved → results/step3_best_sgRNA.json\n")

    return {"best_guide": best}


if __name__ == "__main__":
    try:
        out = run()
        b = out["best_guide"]
        print("[Step 3 OK]")
        print(
            f"Best sgRNA: {b['targetSeq']} | orientation={b['orientation']} | "
            f"distance={b['DistanceFromExon']} | cfd={b['scores'].get('cfdSpecScore')} | "
            f"mit={b['scores'].get('mitSpecScore')} | mm4={b['scores'].get('offtargets_mm4_count')}"
        )
    except Exception as e:
        print(e)
        sys.exit(1)
