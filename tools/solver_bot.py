#!/usr/bin/env python3
"""
BioWordle Solver-Bot & Difficulty Calibration Harness (P0-2 & Appendix A)
=========================================================================
Evaluates the greedy max-entropy information-theoretic solver across the full
puzzle catalog under both ``ranked`` and ``learn`` profiles, prints summary
statistics, and enforces the build-time solver difficulty gates defined in
``data/feedback_spec.json``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict

ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from chem_engine import (  # noqa: E402
    assert_solver_gate,
    compute_seqle_information_metrics,
    run_chemle_greedy_solver,
)


def run_all_gates(puzzles_path: Path, spec_path: Path) -> Dict[str, Any]:
    data = json.loads(puzzles_path.read_text(encoding="utf-8"))
    spec = json.loads(spec_path.read_text(encoding="utf-8"))

    chemle = data["chemle"]
    sim_matrix = data.get("chemle_similarity") or data.get("similarity_matrix")
    seqle = data["seqle"]
    aa_props = data["aa_properties"]

    report: Dict[str, Any] = {"chemle": {}, "seqle": {}}

    for profile_name, profile_spec in spec["profiles"].items():
        stats = run_chemle_greedy_solver(chemle, sim_matrix, profile_spec)
        assert_solver_gate(stats, profile_spec["solver_gate"], profile_name)
        report["chemle"][profile_name] = stats
        print(
            f"[Solver-Bot · Chemle · {profile_name.upper():6s}] "
            f"mean={stats['mean']:.2f}  median={stats['median']:.1f}  "
            f"p95={stats['p95']}  max={stats['max']}  "
            f"<=2-guess={stats['pct_solved_le2']:.2f}%  "
            f"<=6-guess={stats['pct_solved_le6']:.1f}%  "
            f"min_opening_exp_rem={stats['min_opening_expected_remaining']:.2f}  "
            f"dist={stats['distribution']}"
        )

    # Also compute Seqle opening information metrics per profile
    for profile_name, profile_spec in spec["profiles"].items():
        eligible_seq = [
            m for m in seqle if profile_name in m.get("game_eligibility", ["learn"])
        ]
        metrics = compute_seqle_information_metrics(eligible_seq, aa_props)
        min_exp_rem = min(m["expected_remaining_candidates"] for m in metrics.values())
        max_ig = max(m["information_gain_bits"] for m in metrics.values())
        report["seqle"][profile_name] = {
            "eligible_count": len(eligible_seq),
            "min_opening_expected_remaining": round(min_exp_rem, 2),
            "max_opening_ig_bits": round(max_ig, 2),
        }
        print(
            f"[Solver-Bot · Seqle  · {profile_name.upper():6s}] "
            f"eligible={len(eligible_seq)}  "
            f"min_opening_exp_rem={min_exp_rem:.2f}  "
            f"max_opening_ig={max_ig:.2f} bits"
        )

    return report


def main() -> None:
    puzzles_path = ROOT_DIR / "data" / "puzzles.json"
    spec_path = ROOT_DIR / "data" / "feedback_spec.json"
    if not puzzles_path.exists():
        print(f"ERROR: {puzzles_path} not found. Run build_puzzles.py first.", file=sys.stderr)
        sys.exit(1)
    run_all_gates(puzzles_path, spec_path)
    print("[Solver-Bot] All profile difficulty & k-anonymity gates PASSED.")


if __name__ == "__main__":
    main()
