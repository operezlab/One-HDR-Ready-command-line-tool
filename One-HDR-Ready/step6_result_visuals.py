#!/usr/bin/env python3
"""
Step 6 — Final Fixed Annotated Primer Template Visualization
-------------------------------------------------------------
Readable, wrapped sequence visualization (100 bp per line)
with correctly aligned and labeled components:
Left/Right HDR arms, sgRNA, PAM, exon, and primers.

Output:
  results/final_sequence_visualization.png
"""

from __future__ import annotations
import json, sys
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio.Seq import Seq

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)


def run(line_length: int = 100, spacing_x: float = 0.25, spacing_y: float = 0.8):
    # --- Load data ---
    seq_json = RESULTS_DIR / "step1_sequences.json"
    sg_json = RESULTS_DIR / "step3_best_sgRNA.json"
    hdr_json = RESULTS_DIR / "step4_hdr_arms.json"
    primers_json = RESULTS_DIR / "step2_primers.json"

    seqs = json.loads(seq_json.read_text())
    sg = json.loads(sg_json.read_text())
    hdr = json.loads(hdr_json.read_text())
    primers = json.loads(primers_json.read_text()) if primers_json.exists() else {}

    seqs = seqs.get("sequences", seqs)
    sg = sg.get("best_guide", sg)

    primer_seq = seqs.get("seq_primer_template", "")
    exon_seq = seqs.get("seq_exon", "")
    sg_seq = sg.get("targetSeq", "")
    orientation = sg.get("orientation", "forw")
    gene_id = seqs.get("gene_id", "target")

    if not primer_seq or not sg_seq:
        raise ValueError("Missing primer template or sgRNA sequence.")

    # --- Handle sgRNA orientation ---
    sg_to_find = sg_seq
    if orientation.lower().startswith("rev"):
        sg_to_find = str(Seq(sg_seq).reverse_complement())

    sg_start = primer_seq.find(sg_to_find)
    if sg_start == -1:
        print("[Step 6] sgRNA not found in primer template — cannot visualize.")
        sys.exit(0)
    sg_end = sg_start + len(sg_to_find)
    pam_start = sg_end
    pam = (pam_start, pam_start + 2) if primer_seq[pam_start:pam_start + 2] == "GG" else None

    # --- Helper: robust sequence search ---
    def _find_span(template: str, query: str):
        if not template or not query:
            return None
        i = template.find(query)
        if i != -1:
            return (i, i + len(query))
        rc = str(Seq(query).reverse_complement())
        i = template.find(rc)
        if i != -1:
            return (i, i + len(query))
        return None

    # --- Locate features ---
    exon_span = _find_span(primer_seq, exon_seq)
    left_span = _find_span(primer_seq, hdr.get("left_arm", ""))
    right_span = _find_span(primer_seq, hdr.get("right_arm", ""))

    # fallback for arms
    if not left_span and hdr.get("left_arm", ""):
        L = len(hdr["left_arm"])
        left_span = (max(0, sg_start - L), sg_start)
    if not right_span and hdr.get("right_arm", ""):
        R = len(hdr["right_arm"])
        right_span = (sg_end, min(len(primer_seq), sg_end + R))

    # primers
    fwd_prim = primers.get("primers", {}).get("forward", {})
    rev_prim = primers.get("primers", {}).get("reverse", {})
    fwd_span = (fwd_prim.get("start", 0),
                fwd_prim.get("start", 0) + fwd_prim.get("length", 0)) if fwd_prim else None
    rev_span = (rev_prim.get("start", 0),
                rev_prim.get("start", 0) + rev_prim.get("length", 0)) if rev_prim else None

    # --- Prepare layout ---
    num_lines = (len(primer_seq) + line_length - 1) // line_length
    fig_height = max(6, num_lines * 0.5 + 2)
    fig, ax = plt.subplots(figsize=(18, fig_height))
    ax.axis("off")

    # --- draw colored bars and labels ---
    y_bar = 1.3
    def draw_bar(span, color, label):
        if not span:
            return
        s, e = span
        if e > s:
            ax.add_patch(mpatches.Rectangle((s * spacing_x / line_length, y_bar),
                                            (e - s) * spacing_x / line_length, 0.15,
                                            color=color, alpha=0.7))
            ax.text(((s + e) / 2) * spacing_x / line_length, y_bar + 0.25,
                    label, color=color, fontsize=8, ha="center", va="bottom")

    draw_bar(left_span, "skyblue", "Left HDR Arm")
    draw_bar((sg_start, sg_end), "crimson", "sgRNA")
    if pam:
        draw_bar(pam, "orange", "PAM")
    draw_bar(right_span, "lightgreen", "Right HDR Arm")
    if exon_span:
        draw_bar(exon_span, "navy", "Exon")
    draw_bar(fwd_span, "purple", "Forward Primer")
    draw_bar(rev_span, "violet", "Reverse Primer")

    # --- Render sequence (wrapped 100 bp per line) ---
    y_offset = 0
    for line_idx in range(num_lines):
        start = line_idx * line_length
        end = min(start + line_length, len(primer_seq))
        segment = primer_seq[start:end]

        # base text
        for j, base in enumerate(segment):
            pos = start + j
            color = "black"
            if sg_start <= pos < sg_end:
                color = "crimson"
            elif pam and pam[0] <= pos < pam[1]:
                color = "orange"
            ax.text(j * spacing_x, -line_idx * spacing_y, base,
                    fontsize=7, color=color, ha="center", va="center", fontfamily="monospace")

        # left coordinate
        ax.text(-6, -line_idx * spacing_y, f"{start:>4}",
                fontsize=7, color="gray", ha="right", va="center", fontfamily="monospace")

    ax.set_xlim(-10, line_length * spacing_x)
    ax.set_ylim(-num_lines * spacing_y - 0.5, y_bar + 0.8)
    plt.title(f"{gene_id} Annotated Primer Template (5'→3')", fontsize=11, fontweight="bold")

    # --- Legend ---
    handles = [
        mpatches.Patch(color="skyblue", label="Left HDR Arm"),
        mpatches.Patch(color="lightgreen", label="Right HDR Arm"),
        mpatches.Patch(color="navy", label="Exon"),
        mpatches.Patch(color="crimson", label="sgRNA"),
        mpatches.Patch(color="orange", label="PAM"),
        mpatches.Patch(color="purple", label="Forward Primer"),
        mpatches.Patch(color="violet", label="Reverse Primer"),
    ]
    ax.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8)

    plt.tight_layout()
    out_png = RESULTS_DIR / "final_sequence_visualization.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[Step 6] Readable annotated visualization saved → {out_png}")


if __name__ == "__main__":
    try:
        run()
    except Exception as e:
        print(f"[Step 6] ERROR: {e}")
        sys.exit(1)
