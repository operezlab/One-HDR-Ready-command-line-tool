#!/usr/bin/env python3
"""
Step 5 — Final Output Summary for One-HDR-Ready

Combines results from previous steps and produces:
- results/final_summary.txt  (human-readable summary)
- results/final_summary.json (structured version)

Inputs:
  - step1_sequences.json
  - step3_best_sgRNA.json
  - step4_hdr_arms.json
"""

from __future__ import annotations
import json, sys
from pathlib import Path

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

def run():
    seq_json = RESULTS_DIR / "step1_sequences.json"
    sg_json  = RESULTS_DIR / "step3_best_sgRNA.json"
    hdr_json = RESULTS_DIR / "step4_hdr_arms.json"

    if not (seq_json.exists() and sg_json.exists() and hdr_json.exists()):
        raise FileNotFoundError("One or more required step files missing (1, 3, or 4).")

    seqs = json.loads(seq_json.read_text())
    sg   = json.loads(sg_json.read_text())
    hdr  = json.loads(hdr_json.read_text())

    # Handle flat or nested keys gracefully
    seqs = seqs.get("sequences", seqs)
    sg   = sg.get("best_guide", sg)

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

    # Interpret verification status for readability
    if verif_ok is True:
        verif_note = "Verification successful: PAM/codon edits preserved amino acid sequence."
    elif verif_ok is False and hdr.get("distance_from_exon", 0) is not None:
        # no PAM overlap or no edit performed
        verif_note = "No codon/PAM edits were required or attempted for this sgRNA."
    else:
        verif_note = "Verification inconclusive — check sgRNA and exon overlap manually."


    summary = f"""
One-HDR-Ready Final Summary
===========================

Gene ID:        {gene_id}
UniProt ID:     {uniprot_id}
Chromosome:     {chrom}
Exon region:    {exon_start}-{exon_end}
Reverse strand: {seqs.get('reverseStrand')}

Best sgRNA
----------
Sequence:       {sg_seq}
Orientation:    {orientation}
Distance from exon edge: {dist_from_exon} bp

CRISPOR Scores:
  MIT Spec Score:        {scores.get('mitSpecScore')}
  CFD Spec Score:        {scores.get('cfdSpecScore')}
  Doench '16 Score:      {scores.get('Doench16')}
  Moreno-Mateos Score:   {scores.get('MorenoMateos')}
  RuleSet3 Score:        {scores.get('RuleSet3')}
  Out-of-Frame Score:    {scores.get('OutOfFrame')}
  Lindel Score:          {scores.get('Lindel')}

HDR Arms
--------
Left Arm Length:  {left_len} bp
Right Arm Length: {right_len} bp
Verification OK:  {verif_ok}
{verif_note}

Output Files:
  - step1_sequences.json
  - step3_best_sgRNA.json
  - step4_hdr_arms.json
  - step4_hdr_arms.fa
  - final_summary.txt
"""

    # Write summary files
    (RESULTS_DIR / "final_summary.txt").write_text(summary.strip() + "\n")
    combined = {
        "gene_id": gene_id,
        "uniprot_id": uniprot_id,
        "chromosome": chrom,
        "exon_start": exon_start,
        "exon_end": exon_end,
        "sgRNA": sg,
        "hdr_arms": hdr,
    }
    (RESULTS_DIR / "final_summary.json").write_text(json.dumps(combined, indent=2))

    print(summary)
    print("[Step 5] Final summary written to results/final_summary.txt and .json")

    # --- Text-based visualization section ---
    primer_seq = seqs.get("seq_primer_template", "")
    sg = sg_seq or ""
    visual = ""

    if primer_seq and sg:
        idx = primer_seq.find(sg)
        if idx != -1:
            window = 60
            start = max(0, idx - window)
            end = min(len(primer_seq), idx + len(sg) + window)
            sub = list(primer_seq[start:end])

            # Mark sgRNA
            for i in range(idx - start, idx - start + len(sg)):
                if 0 <= i < len(sub):
                    sub[i] = ">"
            # Mark PAM if present immediately downstream
            pam_start = idx + len(sg)
            if pam_start + 2 <= len(primer_seq) and primer_seq[pam_start:pam_start+2] == "GG":
                pam_rel = pam_start - start
                for j in range(pam_rel, min(pam_rel + 2, len(sub))):
                    sub[j] = "P"
            visual = "".join(sub)
        else:
            visual = f"(sgRNA sequence not found in primer template of length {len(primer_seq)})"
    else:
        visual = "(primer template or sgRNA missing)"


if __name__ == "__main__":
    try:
        run()
    except Exception as e:
        print(f"[Step 5] ERROR: {e}")
        sys.exit(1)
