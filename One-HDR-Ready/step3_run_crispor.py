#!/usr/bin/env python3
"""
Step 3 — Run CRISPOR and Select Best sgRNA (automatic)

Inputs
------
- Uses Step 1 outputs:
  - `results/step1_crispor_input.fa` (85 bp: 50u / 35d around stop edge)
  - `results/step1_sequences.json` (has exon + primer template strings)

What this step does
-------------------
1) Runs `crispor.py` (hg38) automatically via subprocess on the FA file.
2) Parses CRISPOR scored TSV.
3) Estimates distance of each sgRNA from the exon edge using the primer template as reference.
4) Picks a "best" sgRNA prioritizing:
   - minimal absolute distance from the exon edge (closer is better)
   - then highest Doench '16 score (if present)
5) Saves:
   - `results/step3_crispor_output.tsv` (raw CRISPOR output)
   - `results/step3_best_sgRNA.json` (selected guide + scores + distance)

Returns
-------
- Dictionary containing `best_guide` and a small summary.
"""

from __future__ import annotations
import os
import sys
import json
import subprocess
from pathlib import Path
from typing import Dict, Any, Optional

import pandas as pd
from Bio.Seq import Seq

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Relative CRISPOR path (confirmed by user)
SCRIPT_DIR = Path(__file__).resolve().parent
CRISPOR_SCRIPT = SCRIPT_DIR / ".." / "crisporWebsite" / "crispor.py"

INPUT_FASTA = RESULTS_DIR / "step1_crispor_input.fa"
SEQS_JSON   = RESULTS_DIR / "step1_sequences.json"
CRISPOR_TSV = RESULTS_DIR / "step3_crispor_output.tsv"

# ----------------- helpers -----------------

def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def _load_sequences() -> Dict[str, Any]:
    if not SEQS_JSON.exists():
        raise FileNotFoundError(f"Missing Step 1 JSON: {SEQS_JSON}")
    data = json.loads(SEQS_JSON.read_text())
    # expect keys: seq_exon, seq_crispor_input, seq_guide_window, seq_primer_template, reverseStrand
    return data


def _find_edge_index(seqs: Dict[str, Any]) -> Optional[int]:
    """Find the exon edge within the primer template.
    For forward (reverseStrand=False): edge ~ last base of exon (3' end) within primer template.
    For reverse (reverseStrand=True): edge ~ first base of exon within primer template.
    Returns index (0-based) within the primer template string.
    """
    exon = seqs["seq_exon"]
    tmpl = seqs["seq_primer_template"]
    reverse = bool(seqs.get("reverseStrand", False))

    if not exon or not tmpl:
        return None

    # Use a 15-mer at the relevant edge for robust matching
    if reverse:
        # edge at exon start (5' in gene sense)
        edge_15 = exon[:15]
    else:
        # edge at exon end (3' in gene sense)
        edge_15 = exon[-15:]

    idx = tmpl.find(edge_15)
    if idx == -1:
        # try longer search/fallbacks
        for k in (12, 10, 8):
            if reverse:
                edge_k = exon[:k]
            else:
                edge_k = exon[-k:]
            idx = tmpl.find(edge_k)
            if idx != -1:
                # point to the last base of that partial edge for forward; first base for reverse
                return idx if reverse else idx + k - 1
        return None
    else:
        # if forward: return last base index; reverse: return first base index
        return idx if reverse else idx + 15 - 1


def _guide_position_in_template(guide_seq: str, orientation: str, tmpl: str) -> Optional[int]:
    """Return 0-based index of the guide *start* in template (gene sense), handling orientation.
    orientation: 'forw' means guide is already in gene-sense; 'rev' means RC first.
    """
    query = guide_seq if orientation == "forw" else reverse_complement(guide_seq)
    pos = tmpl.find(query)
    return None if pos == -1 else pos


def _distance_from_edge(edge_idx: int, guide_start: int) -> int:
    """Simple signed distance (bp) from exon edge to guide start.
    Positive when guide is downstream (to the right) of the edge in primer template string.
    """
    return guide_start - edge_idx


# ----------------- main functionality -----------------

def run() -> Dict[str, Any]:
    # 1) Check inputs and run CRISPOR
    if not INPUT_FASTA.exists():
        raise FileNotFoundError(f"Missing CRISPOR input FASTA: {INPUT_FASTA}")
    if not CRISPOR_SCRIPT.exists():
        raise FileNotFoundError(f"crispor.py not found at {CRISPOR_SCRIPT}")

    # Run crispor.py with hg38
    cmd = [sys.executable, str(CRISPOR_SCRIPT), "hg38", str(INPUT_FASTA), str(CRISPOR_TSV)]
    print("Running CRISPOR:", " ".join(map(str, cmd)))
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"CRISPOR failed (code {proc.returncode})\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    print("CRISPOR finished.")

    if not CRISPOR_TSV.exists():
        raise FileNotFoundError("Expected CRISPOR TSV not produced.")

    # 2) Load sequences context
    seqs = _load_sequences()
    primer_template = seqs.get("seq_primer_template", "")
    edge_idx = _find_edge_index(seqs)

    # 3) Parse CRISPOR TSV
    df = pd.read_csv(CRISPOR_TSV, sep="\t")
    # attempt to keep only needed columns, but guard if missing
    needed_cols = ["guideId", "targetSeq", "mitSpecScore", "cfdSpecScore", "Doench '16-Score",
                   "Moreno-Mateos-Score", "Doench-RuleSet3-Score", "Out-of-Frame-Score", "Lindel-Score"]
    for c in needed_cols:
        if c not in df.columns:
            df[c] = pd.NA

    # Extract location/orientation if present in guideId like 123forw / 456rev
    if "guideId" in df.columns:
        loc_orient = df["guideId"].astype(str).str.extract(r"(\d+)(forw|rev)?")
        if 0 in loc_orient.columns:
            df["location"] = pd.to_numeric(loc_orient[0], errors="coerce")
        if 1 in loc_orient.columns:
            df["orientation"] = loc_orient[1].fillna("forw")
        else:
            df["orientation"] = "forw"
    else:
        df["orientation"] = "forw"

    # 4) Distance estimation relative to exon edge (if we could place the edge)
    distances = []
    for _, row in df.iterrows():
        guide = str(row.get("targetSeq", ""))
        orient = str(row.get("orientation", "forw"))
        if not guide or not primer_template or edge_idx is None:
            distances.append(pd.NA)
            continue
        gpos = _guide_position_in_template(guide, orient, primer_template)
        if gpos is None:
            distances.append(pd.NA)
        else:
            distances.append(_distance_from_edge(edge_idx, gpos))
    df["DistanceFromExon"] = distances

    # 5) Choose best guide: minimal |distance| first, then highest Doench '16
    df_sorted = df.copy()
    # abs distance with NA treated as large
    df_sorted["absD"] = df_sorted["DistanceFromExon"].apply(lambda x: abs(x) if pd.notna(x) else 1_000_000)
    # Doench might be string, ensure numeric
    df_sorted["Doench16_num"] = pd.to_numeric(df_sorted["Doench '16-Score"], errors="coerce")
    df_sorted = df_sorted.sort_values(by=["absD", "Doench16_num"], ascending=[True, False])

    best_row = df_sorted.iloc[0].to_dict() if not df_sorted.empty else {}
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
        },
    }

    (RESULTS_DIR / "step3_best_sgRNA.json").write_text(json.dumps(best, indent=2))
    print("Best sgRNA written to results/step3_best_sgRNA.json")

    return {"best_guide": best}


if __name__ == "__main__":
    try:
        out = run()
        b = out.get("best_guide", {})
        print("[Step 3] CRISPOR run complete.")
        print(f"Best sgRNA: {b.get('targetSeq')}  orientation={b.get('orientation')}  d={b.get('DistanceFromExon')}")
    except Exception as e:
        print(f"[Step 3] ERROR: {e}")
        sys.exit(1)
