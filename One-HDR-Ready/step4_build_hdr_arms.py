#!/usr/bin/env python3
"""
Step 4 — Advanced HDR Arm Builder (codon-aware, PAM-sensitive)

Inputs
------
- Step 1: results/step1_sequences.json   (expects keys incl. seq_exon, seq_guide_window, reverseStrand)
- Step 3: results/step3_best_sgRNA.json   (expects keys: best_guide.targetSeq, orientation, DistanceFromExon)

Behavior
--------
- Extracts left/right HDR arms from seq_guide_window around the chosen sgRNA.
- Handles strand and sgRNA orientation (rev/forw), codon-frame, and PAM edits.
- If PAM overlaps coding bases, performs synonymous codon edits to break PAM while
  preserving amino-acid sequence; verifies AA identity post‑edit.
- Emits detailed logs to stdout (as requested).

Outputs
-------
- results/step4_hdr_arms.json
- results/step4_hdr_arms.fa
"""
from __future__ import annotations
import sys, json
from pathlib import Path
from typing import Dict, Any, Tuple, Optional
from Bio.Seq import Seq
from Bio.Data import CodonTable

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ---------------- Basics ----------------

def reverse_complement(s: str) -> str:
    return str(Seq(s).reverse_complement())

def translate_dna_to_protein(dna_sequence: str) -> str:
    return str(Seq(dna_sequence).translate(to_stop=False))

def trim_to_multiple_of_three(sequence: str) -> str:
    start_index = len(sequence) % 3
    return sequence[start_index:] if start_index else sequence

def codon_order(sequence: str) -> str:
    return ''.join(str((i % 3) + 1) for i in range(len(sequence)))

# ---------------- Sequence search helpers ----------------

def find_and_extract_sequence(gdna: str, target: str, upstream: int, downstream: int) -> Tuple[Optional[str], Optional[str]]:
    """Return (upstream_sequence, downstream_sequence) around *end* of target in gdna.
    Upstream length ~= upstream; downstream length ~= downstream.
    """
    target_index = gdna.find(target)
    if target_index == -1:
        return None, None
    u_seq = gdna[max(0, target_index + len(target) - upstream): target_index + len(target)]
    d_seq = gdna[target_index + len(target) - 1: target_index + len(target) - 1 + downstream]
    return u_seq, d_seq


def find_sgrna_coordinates_in_exon(exon_sequence: str, sgrna_sequence: str) -> Tuple[Optional[int], Optional[int]]:
    start_index = exon_sequence.find(sgrna_sequence)
    if start_index == -1:
        return None, None
    return start_index, start_index + len(sgrna_sequence)


def find_partial_sgrna_coordinates_in_exon(exon_sequence: str, sgrna_sequence: str, min_match_length: int) -> Tuple[Optional[int], Optional[int]]:
    for i in range(len(sgrna_sequence) - min_match_length + 1):
        partial = sgrna_sequence[i:i + min_match_length]
        start_index = exon_sequence.find(partial)
        if start_index != -1:
            return start_index, start_index + min_match_length
    return None, None

# ---------------- Codon editing helpers ----------------

def get_alternative_codon(original_codon: str, amino_acid: str) -> Optional[str]:
    table = CodonTable.unambiguous_dna_by_id[1]  # Standard table
    syn = [codon for codon, aa in table.forward_table.items() if aa == amino_acid and codon != original_codon]
    return syn[0] if syn else None


def modify_sgrna_codons(dna_sequence_trimmed: str, sgrna_start: int, sgrna_end: int,
                         codon_order_sequence: str, left_arm: str, sgRNA_sequence: str) -> Tuple[str, bool]:
    """Synonymously modify codons that span the sgRNA within left_arm; verify AA preserved."""
    # locate sgRNA in left arm
    l_start, l_end = find_sgrna_coordinates_in_exon(left_arm, sgRNA_sequence)
    if l_start is None or l_end is None:
        print("sgRNA sequence not found in the left_arm.")
        return left_arm, False
    print(f"sgRNA sequence found in left_arm at: start={l_start}, end={l_end}")

    sg_in_left = left_arm[l_start:l_end]
    edited = list(left_arm)

    for i in range(0, len(sg_in_left), 3):
        c_start = l_start + i
        c_end = c_start + 3
        codon = left_arm[c_start:c_end]
        if len(codon) != 3:
            print(f"Incomplete codon near end of sgRNA window: '{codon}'")
            continue
        aa = translate_dna_to_protein(codon)
        alt = get_alternative_codon(codon, aa)
        if alt and alt != codon:
            print(f"  Synonymous change: {codon}->{alt} (AA {aa}) at {c_start}-{c_end}")
            edited[c_start:c_end] = list(alt)
        else:
            print(f"  No synonymous alternative for {codon} (AA {aa}); leaving as-is.")

    edited_left = ''.join(edited)

    # verify AA preservation on frame-trimmed sequences
    orig_trim = trim_to_multiple_of_three(left_arm)
    edit_trim = trim_to_multiple_of_three(edited_left)
    if translate_dna_to_protein(orig_trim) == translate_dna_to_protein(edit_trim):
        print("Verification successful: amino-acid sequence unchanged after edit.")
        return edited_left, True
    print("Verification failed: amino-acid sequence changed.")
    return edited_left, False


def modify_sgrna_PAM(left_arm: str, last_three_start: int, last_three_bases: str,
                      dna_sequence_trimmed: str, sgrna_start: int, sgrna_end: int,
                      codon_order_sequence: str, sgRNA_sequence: str) -> Tuple[str, bool]:
    """Try a direct codon swap for the PAM-adjacent codon; verify AA preserved.
    If verification fails, fall back to codon-by-codon edits across the sgRNA.
    """
    aa = translate_dna_to_protein(last_three_bases)
    alt = get_alternative_codon(last_three_bases, aa)
    if alt:
        edited = left_arm[:last_three_start] + alt + left_arm[last_three_start + 3:]
        orig_trim = trim_to_multiple_of_three(left_arm)
        edit_trim = trim_to_multiple_of_three(edited)
        if translate_dna_to_protein(orig_trim) == translate_dna_to_protein(edit_trim):
            print("PAM-adjacent codon swap preserved AA; using edited left arm.")
            return edited, True
        print("PAM swap changed AA; falling back to codon-by-codon edits across sgRNA.")
        return modify_sgrna_codons(dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, left_arm, sgRNA_sequence)
    print("No synonymous alternative for PAM-adjacent codon; falling back to codon-by-codon edits.")
    return modify_sgrna_codons(dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, left_arm, sgRNA_sequence)

# ---------------- Main ----------------

def run() -> Dict[str, Any]:
    seq_json = RESULTS_DIR / "step1_sequences.json"
    sg_json  = RESULTS_DIR / "step3_best_sgRNA.json"
    if not seq_json.exists() or not sg_json.exists():
        raise FileNotFoundError("Missing required Step 1 or Step 3 output files.")

    data = json.loads(seq_json.read_text())
    seqs = data.get("sequences", data) 
    gene_name = seqs.get("gene_id") or seqs.get("uniprot_id")


    best_data = json.loads(sg_json.read_text())
    best = best_data.get("best_guide", best_data)


    # Inputs
    last_exon_seq: str = seqs.get("seq_exon", "")
    gdna_sequence: str = seqs.get("seq_guide_window", "")  # per user: use seq_guide_window
    reverse: bool = bool(seqs.get("reverseStrand", False))

    original_sgRNA: str = best.get("targetSeq", "")
    selected_orientation: str = best.get("orientation", "forw")
    selected_distance_from_exon = best.get("DistanceFromExon", 0)

    if not gdna_sequence or not original_sgRNA or not last_exon_seq:
        raise ValueError("Missing required sequences or guide (seq_guide_window / targetSeq / seq_exon).")

    is_sgRNA_inverted = "Y" if selected_orientation == "rev" else "N"

    # build sgRNA_core (strip PAM depending on orientation)
    if is_sgRNA_inverted == "Y":
        sgRNA_core = reverse_complement(original_sgRNA)[3:]  # remove PAM at RC-front
    else:
        sgRNA_core = original_sgRNA[:-3]  # remove PAM at tail

    # Build initial left arm around last 20 bases of last exon (protein coding edge proxy)
    protein_coding_exon = last_exon_seq[-20:] if len(last_exon_seq) >= 20 else last_exon_seq
    left_arm, _ = find_and_extract_sequence(gdna_sequence, protein_coding_exon, upstream=417, downstream=0)

    # Detailed logs
    print(f"Orientation: {selected_orientation}  inverted={is_sgRNA_inverted}")
    print(f"DistanceFromExon: {selected_distance_from_exon}")
    print(f"Initial left_arm len={len(left_arm) if left_arm else None}")

    # Prepare trimmed exon frame/codon map
    dna_sequence_trimmed = trim_to_multiple_of_three(last_exon_seq)
    protein_sequence = translate_dna_to_protein(dna_sequence_trimmed)
    codon_order_sequence = codon_order(dna_sequence_trimmed)

    # Locate sgRNA in trimmed exon (partial allowed)
    min_match_length = max(8, len(original_sgRNA) // 2)
    sgrna_start, sgrna_end = find_partial_sgrna_coordinates_in_exon(dna_sequence_trimmed, original_sgRNA, min_match_length)
    print(f"sgRNA in trimmed exon (partial) start={sgrna_start} end={sgrna_end}")

    # Branching logic from user’s workflow
    modified_left_arm = None
    verification_status = False

    if selected_distance_from_exon >= -3 or (selected_distance_from_exon <= -3 and sgrna_end is None):
        # Pathway 1: Standard/no special PAM overlap handling needed
        print("Pathway 1: Standard processing (no PAM/codon-3 overlap detected or sgRNA not in trimmed exon)")
    else:
        print("Handling negative distance from exon with PAM-aware logic…")
        # Rebuild sgRNA_sequence to include PAM in consistent orientation for left-arm search
        if is_sgRNA_inverted == "Y":
            sgRNA_sequence = reverse_complement(original_sgRNA[-3:]) + sgRNA_core
        else:
            sgRNA_sequence = sgRNA_core + original_sgRNA[-3:]

        print(f"sgRNA_sequence (with PAM): {sgRNA_sequence}")
        print(f"Last Exon (trimmed): {dna_sequence_trimmed}")

        # Re-evaluate sgRNA in trimmed exon (partial)
        sgrna_start, sgrna_end = find_partial_sgrna_coordinates_in_exon(dna_sequence_trimmed, sgRNA_sequence, min_match_length)
        print(f"sgRNA Coord (trimmed exon): {sgrna_start} {sgrna_end}")

        if is_sgRNA_inverted == "Y":
            position_1 = codon_order_sequence[sgrna_start] if sgrna_start is not None else None
            position_2 = codon_order_sequence[sgrna_start + 1] if sgrna_start is not None else None
            if position_1 == '3' or position_2 == '3':
                print("Pathway: One or both of the first two bases are in codon position 3 (rev orientation).")
                # locate sgRNA_core in left arm (partial)
                la_start, la_end = find_partial_sgrna_coordinates_in_exon(left_arm or '', sgRNA_core, min_match_length)
                print(f"sgRNA in left_arm (partial): {la_start}, {la_end}")
                if la_start is not None:
                    # choose last three bases window relative to partial match
                    last_three_start = la_start - 4
                    last_three_bases = (left_arm or '')[last_three_start: la_start - 1]
                    print(f"Last codon near sgRNA (rev path): {last_three_bases}")
                    if last_three_bases == "TGG":
                        modified_left_arm, verification_status = modify_sgrna_codons(
                            dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, left_arm, sgRNA_sequence
                        )
                    else:
                        modified_left_arm, verification_status = modify_sgrna_PAM(
                            left_arm, last_three_start, last_three_bases, dna_sequence_trimmed,
                            sgrna_start, sgrna_end, codon_order_sequence, sgRNA_sequence
                        )
            else:
                print("Pathway: Neither of the first two bases is position 3 (rev orientation).")
                la_start, la_end = find_partial_sgrna_coordinates_in_exon(left_arm or '', sgRNA_core, min_match_length)
                modified_left_arm, verification_status = modify_sgrna_codons(
                    dna_sequence_trimmed, la_start, la_end, codon_order_sequence, left_arm, sgRNA_sequence
                )
        else:
            # forward orientation path: check last two bases positions
            if sgrna_end is None:
                print("sgRNA not found in trimmed exon; skipping codon-3 logic.")
            else:
                position_2 = codon_order_sequence[sgrna_end - 2]
                position_3 = codon_order_sequence[sgrna_end - 1]
                print(f"Codon positions of last two sgRNA bases: {position_2}, {position_3}")
                if position_2 == '3' or position_3 == '3':
                    print("Pathway: One of the last two bases is in codon position 3 (forward orientation).")
                    la_start, la_end = find_sgrna_coordinates_in_exon(left_arm or '', sgRNA_sequence)
                    if la_start is not None:
                        # choose last-three bases window depending on which base is 3
                        if position_2 == '3':
                            last_three_start = la_start + (la_end - la_start) - 4
                            last_three_bases = (left_arm or '')[last_three_start: la_start + (la_end - la_start) - 1]
                        else:
                            last_three_start = la_start + (la_end - la_start) - 3
                            last_three_bases = (left_arm or '')[last_three_start: la_start + (la_end - la_start)]
                        print(f"Last codon near sgRNA (fwd path): {last_three_bases}")
                        if last_three_bases == "TGG":
                            modified_left_arm, verification_status = modify_sgrna_codons(
                                dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, left_arm, sgRNA_sequence
                            )
                        else:
                            modified_left_arm, verification_status = modify_sgrna_PAM(
                                left_arm, last_three_start, last_three_bases, dna_sequence_trimmed,
                                sgrna_start, sgrna_end, codon_order_sequence, sgRNA_sequence
                            )
                    else:
                        print("sgRNA sequence not found in the left arm.")
                else:
                    print("Pathway: Neither of the last two bases is position 3 (forward orientation).")
                    la_start, la_end = find_sgrna_coordinates_in_exon(left_arm or '', sgRNA_sequence)
                    modified_left_arm, verification_status = modify_sgrna_codons(
                        dna_sequence_trimmed, sgrna_start, sgrna_end, codon_order_sequence, left_arm, sgRNA_sequence
                    )

    # Decide final left arm sequence
    left_arm_seq = modified_left_arm if modified_left_arm else (left_arm or '')

    # Build right arm using sgRNA prefix in gdna (orientation-aware)
    if is_sgRNA_inverted == "Y" and selected_orientation == "rev":
        rev_sg = reverse_complement(sgRNA_sequence if 'sgRNA_sequence' in locals() else original_sgRNA)
        _, right_arm = find_and_extract_sequence(gdna_sequence, rev_sg[:7], upstream=0, downstream=321)
        if right_arm is None:
            _, right_arm = find_and_extract_sequence(gdna_sequence, (sgRNA_sequence if 'sgRNA_sequence' in locals() else original_sgRNA)[:7], 0, 321)
    elif is_sgRNA_inverted == "Y":
        _, right_arm = find_and_extract_sequence(gdna_sequence, (sgRNA_sequence if 'sgRNA_sequence' in locals() else original_sgRNA)[:18], 0, 321)
    else:
        _, right_arm = find_and_extract_sequence(gdna_sequence, (sgRNA_sequence if 'sgRNA_sequence' in locals() else original_sgRNA)[:18], 0, 321)

    print(f"Right_arm length: {len(right_arm) if right_arm else None}")
    print(f"Left_arm final length: {len(left_arm_seq)}  verification={verification_status}")

    # AI1 templates (placeholders copied from user logic)
    left_ai1_template = "tgctggccttttgctcaggatccsnggatccCaaggcggtggaCTCGA"
    if is_sgRNA_inverted == "Y":
        left_ai1_seq = left_ai1_template.replace("s", original_sgRNA).replace("n", left_arm_seq)
    else:
        left_ai1_seq = left_ai1_template.replace("s", reverse_complement(original_sgRNA)).replace("n", left_arm_seq)

    right_ai1_template = (
        "CCTGCGGTGTCTTTGCTTrycatgtGGTTCCATGGTGTAATGGTTAGCACTCTGGACTCTGAATCCAGCGATCCGAGTTCAAATCTCGGTGGAACCTxGTTTTAGAGCTAGAAATAGCAA"
    )
    right_ai1_seq = (
        right_ai1_template.replace("r", right_arm or '')
        .replace("y", original_sgRNA)
        .replace("x", original_sgRNA[:20])
    )

    # Save outputs
    out = {
        "guide_seq": original_sgRNA,
        "orientation": selected_orientation,
        "is_inverted": is_sgRNA_inverted,
        "distance_from_exon": selected_distance_from_exon,
        "left_arm": left_arm_seq,
        "right_arm": right_arm or '',
        "verification_passed": bool(verification_status),
        "left_ai1_seq": left_ai1_seq,
        "right_ai1_seq": right_ai1_seq,
    }

    (RESULTS_DIR / "step4_hdr_arms.json").write_text(json.dumps(out, indent=2))

    fa_lines = [
        f">{gene_name}-Left-Arm\n{left_arm_seq}",
        f">{gene_name}-Right-Arm\n{right_arm or ''}",
        f">{gene_name}-Left-AI1\n{left_ai1_seq}",
        f">{gene_name}-Right-AI1\n{right_ai1_seq}",
    ]
    (RESULTS_DIR / "step4_hdr_arms.fa").write_text("\n".join(fa_lines) + "\n")

    print("[Step 4] HDR arms + AI1 written to results/step4_hdr_arms.fa and step4_hdr_arms.json")
    return out


if __name__ == "__main__":
    try:
        run()
    except Exception as e:
        print(f"[Step 4] ERROR: {e}")
        sys.exit(1)
