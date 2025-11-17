#!/usr/bin/env python3
"""
main.py — One-HDR-Ready Modular Pipeline Orchestrator
-----------------------------------------------------
Executes all pipeline steps sequentially:

  Step 1 → Fetch sequences
  Step 2 → Design and rank primers
  Step 3 → Run CRISPOR automatically
  Step 4 → Build HDR arms
  Step 5 → Verify and summarize
  Step 6 → Generate final annotated visualization (PDF + PNG)

Usage:
    python3 main.py <GENE_SYMBOL or UNIPROT_ID>

All outputs are written to ./results/
"""

import subprocess, sys, os
from pathlib import Path

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)

# --------------------------------------------------------------------------- #
# helper to run each step safely
def run_step(step_name: str, cmd: list[str]):
    print(f"\n{'='*80}")
    print(f"[ One-HDR-Ready ] Running {step_name}")
    print(f"{'='*80}")
    try:
        result = subprocess.run(cmd, check=True, text=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        print(result.stdout)
        if result.stderr.strip():
            print("stderr:", result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"\n[ERROR] {step_name} failed with exit code {e.returncode}")
        print(e.stdout)
        print(e.stderr)
        sys.exit(e.returncode)

# --------------------------------------------------------------------------- #
def main():
    if len(sys.argv) < 2:
        print("Usage: python3 main.py <GENE_SYMBOL or UNIPROT_ID>")
        sys.exit(1)

    gene_id = sys.argv[1]

    steps = [
        ("Step 1 — Fetch Sequences",
         ["python3", "step1_fetch_sequences.py", gene_id]),
        ("Step 2 — Primer Design + Ranking",
         ["python3", "step2_design_primers.py"]),
        ("Step 3 — Run CRISPOR",
         ["python3", "step3_run_crispor.py"]),
        ("Step 4 — Build HDR Arms",
         ["python3", "step4_build_hdr_arms.py"]),
        ("Step 5 — Verification / Summary",
         ["python3", "step5_verify_summary.py"]),
        ("Step 6 — Annotated Visualization",
         ["python3", "step6_result_visuals.py"]),
    ]

    for name, cmd in steps:
        run_step(name, cmd)

    print(f"\n✅ Pipeline completed successfully.")
    print(f"Results and plots saved in: {RESULTS_DIR.resolve()}")
    print("Key outputs:")
    print(" ├─ step1_sequences.json / .fa")
    print(" ├─ step3_crispor_output.tsv")
    print(" ├─ step4_hdr_arms.json / .fa")
    print(" ├─ step5_final_summary.json / .txt")
    print(" └─ final_sequence_visualization.pdf / .png")

# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
