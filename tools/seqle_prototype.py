#!/usr/bin/env python3
"""
Phase 3 Seqle Deeper-Pool Prototype Harness (60 / 90 / 120 Ranked Targets)
===========================================================================
Evaluates candidate corpus sizes (28 baseline, 60, 90, 120), clue schedules
(`opening_charge` vs `charge_after_guess_2`), probe spaces (`catalog_only`,
deterministic legal-probe sample seed=20261114 N=10000, holdout seed=20261115
N=5000), and search strategies (`greedy_entropy` vs `two_step`) against the
Phase 3 promotion criteria:
  - Legal-probe median >= 3.0
  - Legal-probe mean >= 2.6
  - Legal-probe <=2-guess rate <= 45.0%
  - >=98.0% solved by Ranked mode cap (<=6 guesses)
"""

from __future__ import annotations

import copy
import json
import statistics
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple

ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from chem_engine import (  # noqa: E402
    compute_peptide_biochem,
    generate_legal_seqle_probes,
    run_seqle_greedy_solver,
)

PROMOTION_CRITERIA = {
    "legal_probe_median_min": 3.0,
    "legal_probe_mean_min": 2.6,
    "legal_probe_le2_max_pct": 45.0,
    "solved_by_cap_min_pct": 98.0,
}


def _format_charge_band(q: float) -> str:
    if q < -0.5:
        return "Charge Band (pH 7.4): Acidic (Q < -0.5)"
    if q > 0.5:
        return "Charge Band (pH 7.4): Basic (Q > +0.5)"
    return "Charge Band (pH 7.4): Near-Neutral (-0.5 ≤ Q ≤ +0.5)"


def build_prototype_records(
    raw_targets: Sequence[Dict[str, Any]],
    aa_props: Dict[str, Dict[str, Any]],
    schedule_mode: str,
) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    for raw in raw_targets:
        rec = copy.deepcopy(raw)
        bio = compute_peptide_biochem(rec["seq"], aa_props)
        rec.update(bio)
        charge_text = _format_charge_band(bio["net_charge_ph74"])
        rule = rec.get("functional_rule") or rec.get("consensus") or rec.get("pattern") or rec["seq"]
        if schedule_mode == "opening_charge":
            ranked_ladder = [
                {
                    "level": 0,
                    "after_guesses": 0,
                    "label": "Opening Band · Net Charge (pH 7.4)",
                    "text": charge_text,
                },
                {
                    "level": 1,
                    "after_guesses": 2,
                    "label": "Functional Family",
                    "text": f"Family: {rec['family']}",
                },
                {
                    "level": 2,
                    "after_guesses": 4,
                    "label": "Biological Context Requirement",
                    "text": f"Context: {rec['context_requirement']}",
                },
                {
                    "level": 3,
                    "after_guesses": 5,
                    "label": "Motif Rule & Biochemical Role",
                    "text": f"Rule: {rule} — {rec['hint_role']}",
                },
            ]
        elif schedule_mode == "charge_after_guess_2":
            ranked_ladder = [
                {
                    "level": 1,
                    "after_guesses": 2,
                    "label": "Net Charge Band (pH 7.4)",
                    "text": charge_text,
                },
                {
                    "level": 2,
                    "after_guesses": 4,
                    "label": "Functional Family & Biological Context",
                    "text": f"Family: {rec['family']} · Context: {rec['context_requirement']}",
                },
                {
                    "level": 3,
                    "after_guesses": 5,
                    "label": "Motif Rule & Biochemical Role",
                    "text": f"Rule: {rule} — {rec['hint_role']}",
                },
            ]
        else:
            raise ValueError(f"Unknown schedule_mode: {schedule_mode}")

        rec["clue_ladder"] = {
            "ranked": ranked_ladder,
            "learn": rec.get("clue_ladder", {}).get("learn") or ranked_ladder,
        }
        records.append(rec)
    return records


def check_ranked_distinguishability(records: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    seqs = [r["seq"] for r in records]
    dup_seqs = len(seqs) - len(set(seqs))
    # Ranked feedback marks self-guess exact ('CCCCC'), while pairwise Ranked distinguishability
    # across non-self probes holds as long as sequences are distinct. Also track exact Level-0
    # + Level-1 + Level-2 + Level-3 clue collisions to ensure clues alone never leave zero separation.
    clue_sigs: Dict[Tuple[str, ...], List[str]] = defaultdict(list)
    for r in records:
        sig = tuple(c["text"] for c in r["clue_ladder"]["ranked"])
        clue_sigs[sig].append(r["seq"])
    clue_collisions = {str(k): v for k, v in clue_sigs.items() if len(v) > 1}
    return {
        "target_count": len(records),
        "unique_sequences": len(set(seqs)),
        "duplicate_sequences": dup_seqs,
        "distinguishable_under_ranked": dup_seqs == 0,
        "full_ladder_clue_collision_groups": len(clue_collisions),
    }


def summarize_per_family(
    records: Sequence[Dict[str, Any]],
    per_target_results: Dict[str, Tuple[int, bool]],
) -> Dict[str, Any]:
    by_fam: Dict[str, List[Tuple[str, int, bool]]] = defaultdict(list)
    early_wins: List[Dict[str, Any]] = []
    slow_or_unsolved: List[Dict[str, Any]] = []
    for r in records:
        seq = r["seq"]
        steps, solved = per_target_results[seq]
        by_fam[r["family"]].append((seq, steps, solved))
        entry = {"seq": seq, "family": r["family"], "steps": steps, "solved": solved}
        if solved and steps <= 2:
            early_wins.append(entry)
        if (not solved) or steps >= 5:
            slow_or_unsolved.append(entry)

    family_stats: Dict[str, Dict[str, Any]] = {}
    for fam, items in sorted(by_fam.items()):
        solved_steps = [s for _, s, ok in items if ok]
        family_stats[fam] = {
            "count": len(items),
            "solved": len(solved_steps),
            "unsolved": len(items) - len(solved_steps),
            "mean_steps": round(statistics.mean(solved_steps), 2) if solved_steps else None,
            "max_steps": max(solved_steps) if solved_steps else None,
        }
    return {
        "per_family": family_stats,
        "early_le2_count": len(early_wins),
        "early_le2_sample": early_wins[:12],
        "slow_ge5_or_unsolved_count": len(slow_or_unsolved),
        "slow_ge5_or_unsolved_targets": slow_or_unsolved,
    }


def passes_promotion_criteria(stats: Dict[str, Any]) -> Dict[str, Any]:
    med = stats.get("median") or 0.0
    mean = stats.get("mean") or 0.0
    le2 = stats.get("pct_solved_le2", 100.0)
    by_cap = stats.get("pct_solved_by_cap", 0.0)
    checks = {
        "median_ge_3": med >= PROMOTION_CRITERIA["legal_probe_median_min"],
        "mean_ge_2_6": mean >= PROMOTION_CRITERIA["legal_probe_mean_min"],
        "le2_le_45pct": le2 <= PROMOTION_CRITERIA["legal_probe_le2_max_pct"],
        "solved_by_cap_ge_98pct": by_cap >= PROMOTION_CRITERIA["solved_by_cap_min_pct"],
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
    }


def evaluate_prototype_corpus(
    corpus_path: Path,
    puzzles_path: Path,
    spec_path: Path,
    out_report_path: Path | None = None,
) -> Dict[str, Any]:
    corpus = json.loads(corpus_path.read_text(encoding="utf-8"))
    puzzles = json.loads(puzzles_path.read_text(encoding="utf-8"))
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    aa_props = puzzles["aa_properties"]
    ranked_spec = spec["profiles"]["ranked"]

    raw_targets = corpus["targets"]
    legal_sample_10k = generate_legal_seqle_probes(aa_props, seed=20261114, count=10000)
    legal_holdout_5k = generate_legal_seqle_probes(aa_props, seed=20261115, count=5000)

    sizes = [28, 60, 90, 120]
    schedules = ["opening_charge", "charge_after_guess_2"]
    strategies = ["greedy_entropy", "two_step"]

    runs: List[Dict[str, Any]] = []
    any_promoted = False

    for size in sizes:
        slice_raw = raw_targets[:size]
        catalog_seqs = [r["seq"] for r in slice_raw]
        sample_probes = sorted(set(catalog_seqs) | set(legal_sample_10k))
        holdout_probes = sorted(set(catalog_seqs) | set(legal_holdout_5k))

        for schedule_mode in schedules:
            recs = build_prototype_records(slice_raw, aa_props, schedule_mode)
            dist_check = check_ranked_distinguishability(recs)
            assert dist_check["distinguishable_under_ranked"], f"Duplicate targets in slice {size}"

            probe_configs = [
                ("catalog_only", catalog_seqs),
                ("legal_probe_sample_10k", sample_probes),
                ("legal_probe_holdout_5k", holdout_probes),
            ]
            for probe_label, probes in probe_configs:
                for strategy in strategies:
                    t0 = time.perf_counter()
                    res = run_seqle_greedy_solver(
                        recs,
                        aa_props,
                        ranked_spec,
                        profile_id="ranked",
                        condition_on_clues=True,
                        answer_records=recs,
                        probe_seqs=probes,
                        strategy=strategy,
                        strategy_label=f"{strategy}_{probe_label}",
                        include_unconditioned=False,
                    )
                    elapsed = round(time.perf_counter() - t0, 2)
                    fam_diag = summarize_per_family(recs, res["per_target_results"])
                    promo = (
                        passes_promotion_criteria(res)
                        if probe_label.startswith("legal_probe")
                        else None
                    )
                    run_entry = {
                        "pool_size": size,
                        "schedule": schedule_mode,
                        "probe_mode": probe_label,
                        "probe_count": len(probes),
                        "strategy": strategy,
                        "elapsed_sec": elapsed,
                        "distinguishability": dist_check,
                        "target_count": res["target_count"],
                        "solved_count": res["solved_count"],
                        "unsolved_count": res["unsolved_count"],
                        "mode_cap": res["mode_cap"],
                        "pct_solved_by_cap": res["pct_solved_by_cap"],
                        "mean": res["mean"],
                        "median": res["median"],
                        "p95": res["p95"],
                        "max": res["max"],
                        "pct_solved_le2": res["pct_solved_le2"],
                        "pct_solved_le6": res["pct_solved_le6"],
                        "pct_solved_le8": res["pct_solved_le8"],
                        "distribution": res["distribution"],
                        "opening_group_metrics": res["opening_group_metrics"],
                        "family_diagnostics": fam_diag,
                        "promotion_gate": promo,
                    }
                    runs.append(run_entry)
                    promo_str = (
                        f"PROMO={'PASS' if promo['passed'] else 'FAIL'}"
                        if promo is not None
                        else "PROMO=n/a(catalog)"
                    )
                    print(
                        f"[Proto N={size:3d} | {schedule_mode:20s} | {probe_label:22s} | {strategy:14s}] "
                        f"mean={res['mean']:.2f} med={res['median']:.1f} <=2={res['pct_solved_le2']:5.1f}% "
                        f"<=6={res['pct_solved_by_cap']:5.1f}% unsolved={res['unsolved_count']} "
                        f"dist={res['distribution']} {promo_str} ({elapsed:.1f}s)"
                    )

    # Check if any candidate (size >= 60, schedule) passes legal-probe promotion on BOTH sample & holdout under both strategies (or greedy + two_step)
    promoted_configs = []
    for size in [60, 90, 120]:
        for sched in schedules:
            legal_runs = [
                r
                for r in runs
                if r["pool_size"] == size
                and r["schedule"] == sched
                and r["probe_mode"].startswith("legal_probe")
            ]
            if legal_runs and all(r["promotion_gate"]["passed"] for r in legal_runs):
                promoted_configs.append({"pool_size": size, "schedule": sched})
                any_promoted = True

    report = {
        "corpus_total_count": corpus["total_count"],
        "uniprot_release": corpus["uniprot_release"],
        "promotion_criteria": PROMOTION_CRITERIA,
        "any_config_promoted": any_promoted,
        "promoted_configs": promoted_configs,
        "shipped_ranked_decision": (
            "Retain v0.2.2 N=28 Ranked pool with KKMPQ->KEKSD replacement and document legal-probe ceiling as known limitation; blind human playtesting on N>=60 required before Ranked pool expansion."
            if not any_promoted
            else "Candidate passed quantitative solver gates; awaiting human blind playtest."
        ),
        "runs": runs,
    }
    if out_report_path is not None:
        out_report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    return report


def main() -> int:
    corpus_path = ROOT_DIR / "data" / "seqle_prototype_corpus.json"
    puzzles_path = ROOT_DIR / "data" / "puzzles.json"
    spec_path = ROOT_DIR / "data" / "feedback_spec.json"
    out_path = ROOT_DIR / "data" / "seqle_prototype_report.json"

    report = evaluate_prototype_corpus(corpus_path, puzzles_path, spec_path, out_report_path=out_path)
    print(f"\n[Prototype Decision] any_config_promoted={report['any_config_promoted']} | {report['shipped_ranked_decision']}")
    print(f"Saved report -> {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
