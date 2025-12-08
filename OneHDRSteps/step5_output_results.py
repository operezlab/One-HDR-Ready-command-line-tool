#!/usr/bin/env python3
"""
Step 5 — Final Output Summary for One-HDR-Ready

Outputs:
  - <geneID>_summary.txt
  - <geneID>_summary.json
  - final_summary.txt
  - final_summary.json

Pulls data from:
  - step1_sequences.json
  - step2_primers.json
  - step1_primer_template.fa
  - step3_best_sgRNA.json
  - step4_hdr_arms.json
"""

from __future__ import annotations
import json, sys, re
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = PROJECT_ROOT / "results"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def run():
    # === Required input files ===
    seq_json = RESULTS_DIR / "step1_sequences.json"
    sg_json  = RESULTS_DIR / "step3_best_sgRNA.json"
    hdr_json = RESULTS_DIR / "step4_hdr_arms.json"

    # === Optional (but recommended) primer files ===
    primers_json = RESULTS_DIR / "step2_primers.json"
    primer_template_fa = RESULTS_DIR / "step1_primer_template.fa"

    # ---- ERRORS if required files missing ----
    if not seq_json.exists() or not sg_json.exists() or not hdr_json.exists():
        raise FileNotFoundError(
            "Missing required step files (step1, step3, or step4)."
        )

    # ------------------------------------------------------------
    # Load required JSONs
    # ------------------------------------------------------------
    seqs = json.loads(seq_json.read_text())
    sg   = json.loads(sg_json.read_text())
    hdr  = json.loads(hdr_json.read_text())

    # Handle nested keys
    seqs = seqs.get("sequences", seqs)
    sg   = sg.get("best_guide", sg)

    # ------------------------------------------------------------
    # Load optional Step 2 data (primers)
    # ------------------------------------------------------------
    primers = {}
    primer_template = ""
    ranked_primers = None
    top_left = top_right = product_size = None

    if primers_json.exists():
        try:
            primers = json.loads(primers_json.read_text())

            # Step 2 stores top primers under "ranked_primers"
            ranked_primers = primers.get("ranked_primers")
            if ranked_primers and len(ranked_primers) > 0:
                best = ranked_primers[0]
                top_left = best.get("left_seq")
                top_right = best.get("right_seq")
                product_size = best.get("product_size")
        except Exception:
            primers = {}

    # Load primer template
    if primer_template_fa.exists():
        lines = primer_template_fa.read_text().splitlines()
        if len(lines) >= 2:
            primer_template = lines[1].strip()

    # ------------------------------------------------------------
    # Extract normal Step 1/3/4 metadata
    # ------------------------------------------------------------
    gene_id    = seqs.get("gene_id") or "unknown_gene"
    uniprot_id = seqs.get("uniprot_id") or "N/A"
    chrom      = seqs.get("chromosome")
    exon_start = seqs.get("exon_start")
    exon_end   = seqs.get("exon_end")

    sg_seq     = sg.get("targetSeq")
    orientation = sg.get("orientation")
    dist_from_exon = sg.get("DistanceFromExon")
    scores = sg.get("scores", {})

    left_arm  = hdr.get("left_arm")
    right_arm = hdr.get("right_arm")
    verif_ok  = hdr.get("verification_passed")
    left_len  = len(left_arm) if left_arm else 0
    right_len = len(right_arm) if right_arm else 0

    # Human-readable interpretation of verification
    if verif_ok is True:
        verif_note = "Verification successful: PAM/codon edits preserve amino acid sequence."
    elif verif_ok is False:
        verif_note = "No edits were required or attempted for this sgRNA."
    else:
        verif_note = "Verification inconclusive — check sgRNA and exon overlap manually."

    # ------------------------------------------------------------
    # Choose safe gene ID for filenames
    # ------------------------------------------------------------
    safe_id = gene_id if gene_id not in ("", None, "unknown_gene") else uniprot_id
    if not safe_id:
        safe_id = "gene"
    safe_id = re.sub(r"[^A-Za-z0-9_.-]", "_", str(safe_id))

    txt_gene_path  = RESULTS_DIR / f"{safe_id}_summary.txt"
    json_gene_path = RESULTS_DIR / f"{safe_id}_summary.json"
    txt_default    = RESULTS_DIR / "final_summary.txt"
    json_default   = RESULTS_DIR / "final_summary.json"

    # ------------------------------------------------------------
    # Build the long text summary
    # ------------------------------------------------------------
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
Primer Template Length: {len(primer_template) if primer_template else 'N/A'} bp

Top Primer Pair:
  Left Primer:   {top_left}
  Right Primer:  {top_right}
  Product Size:  {product_size}

Primer Template (first 120 bp):
  {primer_template[:120] if primer_template else '(missing)'}

Output Files:
  - {safe_id}_summary.txt
  - {safe_id}_summary.json
  - final_summary.txt
  - final_summary.json
""".strip() + "\n"

    # ------------------------------------------------------------
    # Write text summary files
    # ------------------------------------------------------------
    txt_gene_path.write_text(summary)
    txt_default.write_text(summary)

    # ------------------------------------------------------------
    # Build JSON structured output
    # ------------------------------------------------------------
    combined = {
        "gene_id": gene_id,
        "uniprot_id": uniprot_id,
        "chromosome": chrom,
        "exon_start": exon_start,
        "exon_end": exon_end,
        "sgRNA": sg,
        "hdr_arms": hdr,
        "primers": {
            "template": primer_template,
            "top_pair": {
                "left_primer": top_left,
                "right_primer": top_right,
                "product_size": product_size,
            },
            "all_results": primers,
        },
    }

    json_data = json.dumps(combined, indent=2)
    json_gene_path.write_text(json_data)
    json_default.write_text(json_data)

    print(summary)
    print(f"[Step 5] Final summary written to:")
    print(f"  - {txt_gene_path.name}")
    print(f"  - {json_gene_path.name}")
    print(f"  - final_summary.txt")
    print(f"  - final_summary.json")


# ------------------------------------------------------------
# Main Execution
# ------------------------------------------------------------

if __name__ == "__main__":
    try:
        run()
    except Exception as e:
        print(f"[Step 5] ERROR: {e}")
        sys.exit(1)
