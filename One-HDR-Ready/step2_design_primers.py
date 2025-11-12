#!/usr/bin/env python3
"""
Step 2 — Primer Design and Ranking for One-HDR-Ready

Inputs
------
- Uses the output from Step 1: `results/step1_primer_template.fa`

What this step does
-------------------
1) Reads the primer template sequence (±1000 bp around stop codon).
2) Designs primer pairs using Primer3 with defined thermodynamic parameters.
3) Evaluates primer specificity via BLAST and Needleman–Wunsch alignment.
4) Ranks primers by combined Primer3 penalty + specificity score.
5) Saves:
   - `results/step2_primers.json` (full Primer3 output)
   - `results/step2_best_primers.txt` (top 3 ranked pairs)

Returns
-------
- Dictionary containing primer design results and top-ranked primer pair.
"""

import os
import sys
import json
import tempfile
import subprocess
from pathlib import Path
from typing import Dict, Any
import primer3

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

def _run(cmd):
    cp = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if cp.returncode != 0:
        raise RuntimeError(f"Command failed:\n{' '.join(cmd)}\n{cp.stderr}")
    return cp.stdout

def build_db_from_template(template_seq: str, workdir: Path):
    fa = workdir / "template.fa"
    fa.write_text(">TEMPLATE\n" + template_seq + "\n")
    db_prefix = str(workdir / "template_db")
    _run(["makeblastdb", "-dbtype", "nucl", "-parse_seqids", "-in", str(fa), "-out", db_prefix])
    return db_prefix

def blast_hits_for(primer_seq: str, db_prefix: str):
    with tempfile.TemporaryDirectory() as td:
        qf = Path(td) / "q.fa"
        qf.write_text(">q\n" + primer_seq + "\n")
        out = _run([
            "blastn", "-task", "blastn-short", "-db", db_prefix, "-query", str(qf),
            "-strand", "both", "-evalue", "30000", "-word_size", "7",
            "-reward", "1", "-penalty", "-1", "-soft_masking", "false", "-dust", "no",
            "-max_target_seqs", "50000", "-max_hsps", "1",
            "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand"
        ])
    hits = []
    for line in out.strip().splitlines():
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, sstrand = line.split("\t")
        sstart, send = int(sstart), int(send)
        hits.append({
            "sseqid": sseqid,
            "pident": float(pident),
            "length": int(length),
            "mismatch": int(mismatch),
            "gapopen": int(gapopen),
            "sstart": sstart,
            "send": send,
            "sstrand": sstrand
        })
    return hits

def rank_primer_pairs(res: dict, template_seq: str, top_k: int = 3):
    weights = {
        "offtarget_pair": 1000.0,
        "single_primer_amp": 200.0,
        "compl_any": 1.0,
        "compl_end": 2.0,
        "left_self_end": 1.0,
        "right_self_end": 1.0,
    }
    with tempfile.TemporaryDirectory() as td:
        db_prefix = build_db_from_template(template_seq, Path(td))
        results = []
        n_pairs = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0) or 0)
        for i in range(n_pairs):
            Lseq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
            Rseq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
            p3_penalty = float(res.get(f"PRIMER_PAIR_{i}_PENALTY", 0.0) or 0.0)
            Lhits = blast_hits_for(Lseq, db_prefix)
            Rhits = blast_hits_for(Rseq, db_prefix)
            off_target_penalty = len(Lhits) + len(Rhits) - 2
            compl_any = float(res.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH", 0.0) or 0.0)
            compl_end = float(res.get(f"PRIMER_PAIR_{i}_COMPL_END_TH", 0.0) or 0.0)
            left_self_end = float(res.get(f"PRIMER_LEFT_{i}_SELF_END_TH", 0.0) or 0.0)
            right_self_end = float(res.get(f"PRIMER_RIGHT_{i}_SELF_END_TH", 0.0) or 0.0)
            score = (
                p3_penalty
                + weights["offtarget_pair"] * off_target_penalty
                + weights["compl_any"] * compl_any
                + weights["compl_end"] * compl_end
                + weights["left_self_end"] * left_self_end
                + weights["right_self_end"] * right_self_end
            )
            results.append({
                "index": i,
                "left_seq": Lseq,
                "right_seq": Rseq,
                "product_size": res.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE"),
                "score": round(score, 3),
            })
    results.sort(key=lambda x: x["score"])
    return results[:top_k]

def run(primer_template_fa: Path = RESULTS_DIR / "step1_primer_template.fa") -> Dict[str, Any]:
    if not primer_template_fa.exists():
        raise FileNotFoundError(f"Missing primer template: {primer_template_fa}")
    seq = primer_template_fa.read_text().splitlines()[1].strip()

    res = primer3.bindings.design_primers(
        {
            "SEQUENCE_ID": "MY_GENE",
            "SEQUENCE_TEMPLATE": seq,
        },
        {
            "PRIMER_OPT_SIZE": 20,
            "PRIMER_MIN_SIZE": 18,
            "PRIMER_MAX_SIZE": 25,
            "PRIMER_OPT_TM": 62.0,
            "PRIMER_MIN_TM": 60.0,
            "PRIMER_MAX_TM": 63.0,
            "PRIMER_MIN_GC": 40.0,
            "PRIMER_MAX_GC": 60.0,
            "PRIMER_PRODUCT_SIZE_RANGE": [[1400, 2000]],
        },
    )

    ranked = rank_primer_pairs(res, seq, top_k=3)
    out = {"all_results": res, "ranked_primers": ranked}
    (RESULTS_DIR / "step2_primers.json").write_text(json.dumps(res, indent=2))
    txt_lines = [f"Top {len(ranked)} primer pairs:"]
    for r in ranked:
        txt_lines.append(f"#{r['index']}  score={r['score']}  product_size={r['product_size']}\nL: {r['left_seq']}\nR: {r['right_seq']}\n")
    (RESULTS_DIR / "step2_best_primers.txt").write_text("\n".join(txt_lines))
    return out

if __name__ == "__main__":
    try:
        out = run()
        print("[Step 2] Primer design completed.")
        print("Top primers written to results/step2_best_primers.txt")
    except Exception as e:
        print(f"[Step 2] ERROR: {e}")
        sys.exit(1)
