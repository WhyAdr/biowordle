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
    generate_legal_seqle_probes,
    run_chemle_greedy_solver,
    run_seqle_greedy_solver,
)


def _format_solver_summary(label: str, stats: Dict[str, Any]) -> str:
    mean_s = f"{stats['mean']:.2f}" if stats.get("mean") is not None else "n/a"
    med_s = f"{stats['median']:.1f}" if stats.get("median") is not None else "n/a"
    probe_n = stats.get("probe_count", stats["target_count"])
    return (
        f"[{label}] "
        f"targets={stats['target_count']} probes={probe_n} solved={stats['solved_count']} unsolved={stats['unsolved_count']} "
        f"cap={stats['mode_cap']} ({stats['pct_solved_by_cap']:.1f}% solved) "
        f"mean={mean_s} median={med_s} p95={stats['p95']} max={stats['max']} "
        f"<=2-guess={stats['pct_solved_le2']:.2f}% <=6-guess={stats['pct_solved_le6']:.1f}% "
        f"<=8-guess={stats['pct_solved_le8']:.1f}% dist={stats['distribution']}"
    )


def run_all_gates(
    puzzles_path: Path,
    spec_path: Path,
    include_legal_probe_sample: bool = True,
) -> Dict[str, Any]:
    data = json.loads(puzzles_path.read_text(encoding="utf-8"))
    spec = json.loads(spec_path.read_text(encoding="utf-8"))

    chemle = data["chemle"]
    sim_matrix = data.get("chemle_similarity") or data.get("similarity_matrix")
    seqle = data["seqle"]
    aa_props = data["aa_properties"]

    report: Dict[str, Any] = {"chemle": {}, "seqle": {}, "seqle_legal_probe": {}}

    for profile_name, profile_spec in spec["profiles"].items():
        pool_profile = "learn" if profile_spec.get("chemle_only") else profile_name
        ans_ids = [
            r["id"]
            for r in chemle
            if pool_profile in r.get("game_eligibility", ["learn", "ranked"])
        ]
        probe_ids = [
            r["id"]
            for r in chemle
            if pool_profile in r.get("guess_eligibility", r.get("game_eligibility", ["learn", "ranked"]))
        ]
        stats = run_chemle_greedy_solver(
            chemle,
            sim_matrix,
            profile_spec,
            profile_id=profile_name,
            condition_on_clues=True,
            answer_ids=ans_ids,
            probe_ids=probe_ids,
        )
        assert_solver_gate(stats, profile_spec["solver_gate"], profile_name)
        report["chemle"][profile_name] = stats
        max_sing = (stats.get("opening_probe_audit") or {}).get("max_singleton_fraction")
        sing_note = f"  opening_max_singleton={max_sing:.4f}" if max_sing is not None else ""
        print(
            _format_solver_summary(f"Solver-Bot · Chemle · {profile_name.upper():9s}", stats)
            + f"  whole_pool_min_exp_rem={stats['min_opening_expected_remaining']:.2f}{sing_note}"
        )

    catalog_seqs = [r["seq"] for r in seqle]
    legal_probes = (
        sorted(set(catalog_seqs) | set(generate_legal_seqle_probes(aa_props, seed=20261114, count=10000)))
        if include_legal_probe_sample
        else catalog_seqs
    )

    for profile_name, profile_spec in spec["profiles"].items():
        if profile_spec.get("chemle_only"):
            continue
        seq_stats = run_seqle_greedy_solver(
            seqle, aa_props, profile_spec, profile_id=profile_name, condition_on_clues=True
        )
        assert_seqle_solver_gate(seq_stats, profile_spec["seqle_solver_gate"], profile_name)
        report["seqle"][profile_name] = seq_stats
        print(
            _format_solver_summary(
                f"Solver-Bot · Seqle · {profile_name.upper():6s} · {seq_stats['search_method']}",
                seq_stats,
            )
            + f"  yigsr_distinct={seq_stats['yigsr_distinct_signatures']}"
            + f"  max_probe_distinct={seq_stats['max_single_probe_distinct_targets']}"
            + f"  whole_pool_min_exp_rem={seq_stats['min_opening_expected_remaining']:.2f}"
        )
        for og in seq_stats["opening_group_metrics"]:
            print(
                f"    OpeningGroup[{og['group_size']} targets | {og['probe_count']} probes | {og['group_label']}] -> "
                f"best_probe={og['best_probe']} exp_surv={og['best_expected_survivors']:.3f} "
                f"singletons={og['singleton_count']}/{og['group_size']} ({og['singleton_coverage_pct']:.1f}%) "
                f"2-step_P={og['two_step_success_prob']:.3f}"
            )

        if include_legal_probe_sample:
            lp_stats = run_seqle_greedy_solver(
                seqle,
                aa_props,
                profile_spec,
                profile_id=profile_name,
                condition_on_clues=True,
                probe_seqs=legal_probes,
                strategy="greedy_entropy",
                strategy_label="greedy_entropy_sampled_upper_bound_seed20261114",
            )
            lp_two_step = run_seqle_greedy_solver(
                seqle,
                aa_props,
                profile_spec,
                profile_id=profile_name,
                condition_on_clues=True,
                probe_seqs=legal_probes,
                strategy="two_step",
                strategy_label="two_step_sampled_upper_bound_seed20261114",
            )
            report["seqle_legal_probe"][profile_name] = {
                "greedy_entropy": lp_stats,
                "two_step": lp_two_step,
            }
            print(_format_solver_summary(f"LegalProbe · Seqle · {profile_name.upper():6s} · {lp_stats['search_method']}", lp_stats))
            print(_format_solver_summary(f"LegalProbe · Seqle · {profile_name.upper():6s} · {lp_two_step['search_method']}", lp_two_step))
            for og in lp_stats["opening_group_metrics"]:
                print(
                    f"    OpeningGroup-LegalProbe[{og['group_size']} targets | {og['probe_count']} probes | {og['group_label']}] -> "
                    f"best_probe={og['best_probe']} exp_surv={og['best_expected_survivors']:.3f} "
                    f"singletons={og['singleton_count']}/{og['group_size']} ({og['singleton_coverage_pct']:.1f}%) "
                    f"2-step_P={og['two_step_success_prob']:.3f}"
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
