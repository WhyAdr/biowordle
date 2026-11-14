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
    assert_seqle_solver_gate,
    assert_solver_gate,
    run_chemle_greedy_solver,
    run_seqle_greedy_solver,
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
        stats = run_chemle_greedy_solver(
            chemle, sim_matrix, profile_spec, profile_id=profile_name, condition_on_clues=True
        )
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

    for profile_name, profile_spec in spec["profiles"].items():
        seq_stats = run_seqle_greedy_solver(
            seqle, aa_props, profile_spec, profile_id=profile_name, condition_on_clues=True
        )
        assert_seqle_solver_gate(seq_stats, profile_spec["seqle_solver_gate"], profile_name)
        report["seqle"][profile_name] = seq_stats
        print(
            f"[Solver-Bot · Seqle  · {profile_name.upper():6s}] "
            f"eligible={seq_stats['eligible_count']}  "
            f"mean={seq_stats['mean']:.2f}  median={seq_stats['median']:.1f}  max={seq_stats['max']}  "
            f"<=2-guess={seq_stats['pct_solved_le2']:.2f}%  "
            f"yigsr_distinct={seq_stats['yigsr_distinct_signatures']}  "
            f"max_probe_distinct={seq_stats['max_single_probe_distinct_targets']}  "
            f"min_opening_exp_rem={seq_stats['min_opening_expected_remaining']:.2f}  "
            f"dist={seq_stats['distribution']}"
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
