#!/usr/bin/env python3
"""
Step 5 — Final Output Summary for One-HDR-Ready

Outputs:
  - <geneID>_summary.txt
  - <geneID>_summary.json
  - final_summary.txt
  - final_summary.json

Uses inputs:
  - step1_sequences.json
  - step1_primer_template.fa
  - step2_best_primers.txt
  - step3_best_sgRNA.json
  - step4_hdr_arms.json
"""

from __future__ import annotations
import json, sys, re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# MAIN
# ============================================================

def run():

    # Required step outputs
    seq_json = RESULTS_DIR / "step1_sequences.json"
    sg_json  = RESULTS_DIR / "step3_best_sgRNA.json"
    hdr_json = RESULTS_DIR / "step4_hdr_arms.json"

    # Correct primer source file
    primers_txt = RESULTS_DIR / "step2_best_primers.txt"
    primer_template_fa = RESULTS_DIR / "step1_primer_template.fa"

    # Ensure required files exist
    if not seq_json.exists() or not sg_json.exists() or not hdr_json.exists():
        raise FileNotFoundError("Missing required step files (1, 3, or 4).")

    # Load Step 1, 3, 4 JSON data
    seqs = json.loads(seq_json.read_text())
    sg   = json.loads(sg_json.read_text())
    hdr  = json.loads(hdr_json.read_text())

    seqs = seqs.get("sequences", seqs)
    sg   = sg.get("best_guide", sg)

    # ============================================================
    # LOAD PRIMER TEMPLATE (step1_primer_template.fa)
    # ============================================================
    primer_template = ""
    if primer_template_fa.exists():
        lines = primer_template_fa.read_text().splitlines()
        if len(lines) > 1:
            primer_template = lines[1].strip()

    # ============================================================
    # LOAD ACTUAL PRIMER PAIRS (step2_best_primers.txt)
    # ============================================================
    top_left = top_right = None
    product_size = None

    if primers_txt.exists():
        lines = primers_txt.read_text().splitlines()

        for line in lines:
            if "product_size=" in line:
               product_size = int(line.split("product_size=")[1].strip())

            if line.startswith("L: "):
                top_left = line.replace("L: ", "").strip()

            if line.startswith("R: "):
                top_right = line.replace("R: ", "").strip()

    else:
        print("[Step 5] WARNING: No step2_best_primers.txt found.")
        top_left = top_right = "(missing)"

    # ============================================================
    # Extract metadata from Steps 1, 3, 4
    # ============================================================
    gene_id    = seqs.get("gene_id") or "unknown_gene"
    uniprot_id = seqs.get("uniprot_id") or "N/A"
    chrom      = seqs.get("chromosome")
    exon_start = seqs.get("exon_start")
    exon_end   = seqs.get("exon_end")

    sg_seq = sg.get("targetSeq")
    orientation = sg.get("orientation")
    dist_from_exon = sg.get("DistanceFromExon")
    scores = sg.get("scores", {})

    left_arm  = hdr.get("left_arm")
    right_arm = hdr.get("right_arm")
    left_len  = len(left_arm) if left_arm else 0
    right_len = len(right_arm) if right_arm else 0

    verif_ok = hdr.get("verification_passed")

    if verif_ok is True:
        verif_note = "Verification successful: PAM/codon edits preserve amino acid sequence."
    elif verif_ok is False:
        verif_note = "No codon edits required or attempted."
    else:
        verif_note = "Verification inconclusive."

    # ============================================================
    # SAFE FILENAME BASE
    # ============================================================
    safe_id = gene_id if gene_id not in ("unknown_gene", None, "") else uniprot_id
    safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", str(safe_id))

    txt_gene_path  = RESULTS_DIR / f"{safe_id}_summary.txt"
    json_gene_path = RESULTS_DIR / f"{safe_id}_summary.json"
    txt_default    = RESULTS_DIR / "final_summary.txt"
    json_default   = RESULTS_DIR / "final_summary.json"

    # ============================================================
    # BUILD FULL TEXT SUMMARY
    # ============================================================
    summary = f"""
One-HDR-Ready Final Summary
===========================

Gene ID:        {gene_id}
UniProt ID:     {uniprot_id}
Chromosome:     {chrom}
Exon region:    {exon_start}-{exon_end}
Reverse strand: {seqs.get('reverseStrand')}

Best sgRNA (Step 3)
-------------------
Sequence:              {sg_seq}
Orientation:           {orientation}
Distance from exon:    {dist_from_exon} bp

CRISPOR Scores:
  MIT Spec Score:      {scores.get('mitSpecScore')}
  CFD Spec Score:      {scores.get('cfdSpecScore')}
  Doench '16 Score:    {scores.get('Doench16')}
  Moreno-Mateos Score: {scores.get('MorenoMateos')}
  RuleSet3 Score:      {scores.get('RuleSet3')}
  Out-of-Frame Score:  {scores.get('OutOfFrame')}
  Lindel Score:        {scores.get('Lindel')}

HDR Arms (Step 4)
-----------------
Left Arm Length:   {left_len} bp
Right Arm Length:  {right_len} bp
Verification OK:   {verif_ok}
{verif_note}

Primers (Step 2)
----------------
Top Left Primer:     {top_left}
Top Right Primer:    {top_right}
Product Size:        {product_size}

Primer Template Length: {len(primer_template) if primer_template else 'N/A'} bp
Primer Template (first 120 bp):
  {primer_template[:120] if primer_template else '(missing)'}

Output Files:
  - {safe_id}_summary.txt
  - {safe_id}_summary.json
  - final_summary.txt
  - final_summary.json
""".strip() + "\n"

    # ============================================================
    # WRITE SUMMARY FILES
    # ============================================================
    txt_gene_path.write_text(summary)
    txt_default.write_text(summary)

    combined = {
        "gene_id": gene_id,
        "uniprot_id": uniprot_id,
        "chromosome": chrom,
        "exon_start": exon_start,
        "exon_end": exon_end,
        "sgRNA": sg,
        "hdr_arms": hdr,
        "primers": {
            "left_primer": top_left,
            "right_primer": top_right,
            "product_size": product_size,
            "template": primer_template,
        }
    }

    json_gene_path.write_text(json.dumps(combined, indent=2))
    json_default.write_text(json.dumps(combined, indent=2))

    print(summary)
    print(f"[Step 5] Summary written to {txt_gene_path.name} and {json_gene_path.name}")


# ============================================================
# EXECUTION
# ============================================================

if __name__ == "__main__":
    try:
        run()
    except Exception as e:
        print(f"[Step 5] ERROR: {e}")
        sys.exit(1)
