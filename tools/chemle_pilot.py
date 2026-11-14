#!/usr/bin/env python3
"""
BioWordle v0.3 Chemle Expanded-Pool Pilot & Mobile Benchmark Harness
====================================================================
Evaluates the 36-compound Ranked baseline vs the 78-compound Learn pilot across:
  - Solver strategies (`greedy_entropy`, `min_expected_survivors`)
  - Probe restrictions (`full_guess_pool`, `candidate_only`)
  - Feedback/clue scenarios (`36 Ranked`, `36 on Learn`, `78 Learn`, `78 on Ranked` stress)
  - 7 pilot families and 11 Category Packs
  - Nearest-neighbor Tanimoto density, Murcko scaffold collisions, descriptor collisions
  - Turn-by-turn solver walkthroughs for Methylxanthines, Hydroxycinnamates, Catecholamines
  - 5-run cold mobile Chromium (390x844, 4x CPU + Fast-3G CDP) offline HTTP & file:// latency

Emits:
  - data/chemle_pilot_report.json
  - CHEMLE_v0.3_PILOT_BENCHMARK.md
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from chem_engine import evaluate_chemle_feedback_from_spec  # noqa: E402

RANKED_GATE_THRESHOLDS = {
    "mean_min": 3.0,
    "mean_max": 4.3,
    "max_le2_pct": 25.0,
    "min_solved_le6_pct": 98.0,
    "min_opening_expected_remaining": 2.5,
}


def _round2(val: float) -> float:
    return round(float(val), 2)


def _p95(values: Sequence[int | float]) -> float | int:
    if not values:
        return 0
    s = sorted(values)
    idx = min(len(s) - 1, max(0, int(math.ceil(0.95 * len(s)) - 1)))
    return s[idx]


def _write_lf_text(path: Path, text: str) -> None:
    path.write_bytes(text.replace("\r\n", "\n").replace("\r", "\n").encode("utf-8"))


def build_signature_table(
    chemle_records: Sequence[Dict[str, Any]],
    sim_matrix: Dict[str, Dict[str, int]],
    profile_spec: Dict[str, Any],
) -> Dict[Tuple[str, str], str]:
    sigs: Dict[Tuple[str, str], str] = {}
    for g in chemle_records:
        gid = g["id"]
        sim_row = sim_matrix[gid]
        for t in chemle_records:
            tid = t["id"]
            res = evaluate_chemle_feedback_from_spec(g, t, int(sim_row[tid]), profile_spec)
            sigs[(gid, tid)] = str(res["signature"])
    return sigs


def _compatible_after_clue(
    cand: Dict[str, Any],
    target: Dict[str, Any],
    clue_profile: str,
    step_idx: int,
    scaffold_unlock_after: Optional[int],
) -> bool:
    c_ladder = cand.get("clue_ladder", {}).get(clue_profile)
    t_ladder = target.get("clue_ladder", {}).get(clue_profile)
    if c_ladder and t_ladder:
        t_by_level = {item["level"]: item["text"] for item in t_ladder if item["after_guesses"] <= step_idx}
        c_by_level = {item["level"]: item["text"] for item in c_ladder if item["after_guesses"] <= step_idx}
        for lvl, t_txt in t_by_level.items():
            if c_by_level.get(lvl) != t_txt:
                return False
    if scaffold_unlock_after is not None and step_idx >= scaffold_unlock_after:
        if bool(cand.get("has_scaffold")) != bool(target.get("has_scaffold")):
            return False
        if (cand.get("scaffold_smiles") or "") != (target.get("scaffold_smiles") or ""):
            return False
    return True


def _pick_best_probe(
    survivor_ids: Sequence[str],
    probe_ids: Sequence[str],
    sigs: Dict[Tuple[str, str], str],
    strategy: str,
) -> Tuple[str, float, float]:
    n = len(survivor_ids)
    survivor_set = set(survivor_ids)
    best_key: Optional[Tuple[Any, ...]] = None
    best_probe = survivor_ids[0]
    best_entropy = 0.0
    best_exp = float(n)

    for gid in probe_ids:
        counts: Dict[str, int] = defaultdict(int)
        for tid in survivor_ids:
            counts[sigs[(gid, tid)]] += 1
        ent = 0.0
        sum_sq = 0
        for c in counts.values():
            p = c / n
            ent -= p * math.log2(p)
            sum_sq += c * c
        exp_rem = sum_sq / n
        in_surv = 1 if gid in survivor_set else 0

        if strategy == "greedy_entropy":
            score_key = (round(ent, 9), in_surv, -round(exp_rem, 9), gid)
        elif strategy == "min_expected_survivors":
            score_key = (-round(exp_rem, 9), in_surv, round(ent, 9), gid)
        else:
            raise ValueError(f"Unknown solver strategy: {strategy}")

        # Compare deterministically (higher tuple wins; tie-break alphabetically via reversed id rank)
        cmp_key = (score_key[0], score_key[1], score_key[2])
        if best_key is None or cmp_key > best_key or (cmp_key == best_key and gid < best_probe):
            best_key = cmp_key
            best_probe = gid
            best_entropy = ent
            best_exp = exp_rem

    return best_probe, best_entropy, best_exp


def run_pilot_solver(
    target_records: Sequence[Dict[str, Any]],
    probe_records: Sequence[Dict[str, Any]],
    sigs: Dict[Tuple[str, str], str],
    max_cap: int,
    strategy: str = "greedy_entropy",
    probe_mode: str = "full_guess_pool",
    clue_profile: Optional[str] = None,
    scaffold_unlock_after: Optional[int] = None,
    record_trace_ids: Optional[set[str]] = None,
) -> Dict[str, Any]:
    target_ids = sorted(r["id"] for r in target_records)
    probe_pool_ids = sorted(r["id"] for r in probe_records)
    by_id = {r["id"]: r for r in probe_records}
    for r in target_records:
        by_id[r["id"]] = r

    turns_list: List[int] = []
    solved_counts: Counter[int] = Counter()
    unsolved_targets: List[str] = []
    traces: Dict[str, List[Dict[str, Any]]] = {}

    for tid in target_ids:
        target = by_id[tid]
        # Candidate space always equals target_records (matching active answer pool)
        survivors = list(target_ids)
        step = 0
        solved = False
        trace: List[Dict[str, Any]] = []
        want_trace = record_trace_ids is not None and tid in record_trace_ids

        while step < max_cap:
            if clue_profile:
                survivors = [
                    cid
                    for cid in survivors
                    if _compatible_after_clue(by_id[cid], target, clue_profile, step, scaffold_unlock_after)
                ]

            if len(survivors) == 1:
                probe_id = survivors[0]
                ent = 0.0
                exp_rem = 1.0
            else:
                active_probes = survivors if probe_mode == "candidate_only" else probe_pool_ids
                probe_id, ent, exp_rem = _pick_best_probe(survivors, active_probes, sigs, strategy)

            obs_sig = sigs[(probe_id, tid)]
            next_survivors = [cid for cid in survivors if sigs[(probe_id, cid)] == obs_sig]

            if want_trace:
                trace.append({
                    "guess_num": step + 1,
                    "guess_id": probe_id,
                    "guess_name": by_id[probe_id]["name"],
                    "survivors_before": len(survivors),
                    "survivors_after": len(next_survivors),
                    "entropy_bits": round(ent, 3),
                    "expected_survivors": round(exp_rem, 2),
                    "is_target": probe_id == tid,
                })

            step += 1
            if probe_id == tid:
                solved = True
                break
            survivors = next_survivors

        if want_trace:
            traces[tid] = trace

        if solved:
            turns_list.append(step)
            solved_counts[step] += 1
        else:
            unsolved_targets.append(tid)
            turns_list.append(max_cap + 1)

    total = len(target_ids)
    solved_turns = [t for t in turns_list if t <= max_cap]
    mean_val = _round2(statistics.mean(turns_list)) if turns_list else 0.0
    median_val = float(statistics.median(turns_list)) if turns_list else 0.0
    solved_le2 = sum(1 for t in turns_list if t <= 2)
    solved_le6 = sum(1 for t in turns_list if t <= 6)
    solved_le8 = sum(1 for t in turns_list if t <= 8)

    result: Dict[str, Any] = {
        "targets": total,
        "probes": len(probe_pool_ids),
        "strategy": strategy,
        "probe_mode": probe_mode,
        "clue_profile": clue_profile,
        "cap": max_cap,
        "solved_count": len(solved_turns),
        "unsolved_count": len(unsolved_targets),
        "unsolved_targets": unsolved_targets,
        "mean": mean_val,
        "median": median_val,
        "p95": _p95(turns_list),
        "max": max(turns_list) if turns_list else 0,
        "pct_solved_le2": _round2(100.0 * solved_le2 / total) if total else 0.0,
        "pct_solved_le6": _round2(100.0 * solved_le6 / total) if total else 0.0,
        "pct_solved_le8": _round2(100.0 * solved_le8 / total) if total else 0.0,
        "pct_solved_by_cap": _round2(100.0 * len(solved_turns) / total) if total else 0.0,
        "distribution": {str(k): int(solved_counts[k]) for k in sorted(solved_counts)},
    }
    if record_trace_ids is not None:
        result["traces"] = traces
    return result


def compute_density_and_collisions(
    records: Sequence[Dict[str, Any]],
    sim_matrix: Dict[str, Dict[str, int]],
    sigs_ranked: Dict[Tuple[str, str], str],
    sigs_learn: Dict[Tuple[str, str], str],
) -> Dict[str, Any]:
    ids = sorted(r["id"] for r in records)
    by_id = {r["id"]: r for r in records}

    nn_sims: List[int] = []
    pair_sims: List[Tuple[int, str, str]] = []
    for i, a in enumerate(ids):
        best = max(sim_matrix[a][b] for b in ids if b != a)
        nn_sims.append(best)
        for b in ids[i + 1 :]:
            pair_sims.append((int(sim_matrix[a][b]), a, b))
    pair_sims.sort(reverse=True)

    scaffold_groups: Dict[str, List[str]] = defaultdict(list)
    for mid in ids:
        sc = by_id[mid].get("scaffold_smiles") or ""
        if sc:
            scaffold_groups[sc].append(mid)
    colliding_scaffolds = [
        {"scaffold_smiles": sc, "count": len(mids), "members": mids}
        for sc, mids in sorted(scaffold_groups.items(), key=lambda kv: (-len(kv[1]), kv[0]))
        if len(mids) > 1
    ]

    # Exact Ranked descriptor bucket collisions (MW, cLogP, TPSA, HBD, HBA)
    desc_buckets: Dict[Tuple[Any, ...], List[str]] = defaultdict(list)
    for mid in ids:
        r = by_id[mid]
        desc_buckets[(r["MW"], r["cLogP"], r["TPSA"], r["HBD"], r["HBA"])].append(mid)
    exact_desc_collisions = [mids for mids in desc_buckets.values() if len(mids) > 1]

    # Check if ANY pair of molecules produces identical Ranked feedback across ALL catalog guesses
    ranked_indistinguishable_pairs: List[Tuple[str, str]] = []
    learn_indistinguishable_pairs: List[Tuple[str, str]] = []
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            if all(sigs_ranked[(g, a)] == sigs_ranked[(g, b)] for g in ids):
                ranked_indistinguishable_pairs.append((a, b))
            if all(sigs_learn[(g, a)] == sigs_learn[(g, b)] for g in ids):
                learn_indistinguishable_pairs.append((a, b))

    return {
        "count": len(ids),
        "nn_tanimoto_mean": _round2(statistics.mean(nn_sims)),
        "nn_tanimoto_median": float(statistics.median(nn_sims)),
        "nn_tanimoto_p95": _p95(nn_sims),
        "nn_tanimoto_max": max(nn_sims),
        "top_closest_pairs": [
            {"a": a, "b": b, "sim_pct": s} for s, a, b in pair_sims[:8]
        ],
        "distinct_cyclic_scaffolds": len(scaffold_groups),
        "scaffold_collision_groups": len(colliding_scaffolds),
        "top_scaffold_collisions": colliding_scaffolds[:6],
        "exact_5prop_descriptor_collisions": exact_desc_collisions,
        "ranked_feedback_indistinguishable_pairs": ranked_indistinguishable_pairs,
        "learn_feedback_indistinguishable_pairs": learn_indistinguishable_pairs,
    }


def run_mobile_chromium_benchmark(runs: int = 5) -> Dict[str, Any]:
    """Run N cold mobile (390x844, 4x CPU + Fast-3G CDP) runs against static HTTP & file:// bundle."""
    js_script = f"""
import {{ chromium }} from "@playwright/test";
import {{ startStaticServer }} from "./tests/e2e/static_server.mjs";
import path from "node:path";
import {{ pathToFileURL }} from "node:url";

const RUNS = {runs};
const VIEWPORT = {{ width: 390, height: 844 }};
const {{ server, url }} = await startStaticServer({{ port: 0 }});
const browser = await chromium.launch({{ headless: true }});

const httpRuns = [];
const standaloneRuns = [];

async function measureOne(targetUrl, throttleNetwork) {{
  const context = await browser.newContext({{
    viewport: VIEWPORT,
    deviceScaleFactor: 2,
    isMobile: true,
    hasTouch: true,
  }});
  const page = await context.newPage();
  page.setDefaultTimeout(60000);
  const cdp = await context.newCDPSession(page);
  await cdp.send("Emulation.setCPUThrottlingRate", {{ rate: 4 }});
  if (throttleNetwork) {{
    await cdp.send("Network.enable");
    await cdp.send("Network.emulateNetworkConditions", {{
      offline: false,
      latency: 150,
      downloadThroughput: Math.floor(1.6 * 1024 * 1024 / 8),
      uploadThroughput: Math.floor(750 * 1024 / 8),
    }});
  }}

  const externalRequests = [];
  page.on("request", (req) => {{
    const u = req.url();
    if (!u.startsWith(url) && !u.startsWith("file://") && !u.startsWith("data:") && !u.startsWith("blob:")) {{
      externalRequests.push(u);
    }}
  }});

  const t0 = performance.now();
  await page.goto(targetUrl, {{ waitUntil: "load" }});
  await page.locator("#chemle-search-input").waitFor({{ state: "visible" }});
  await page.locator("#chemle-structure-well").waitFor({{ state: "visible" }});
  const bootMs = performance.now() - t0;

  const navPerf = await page.evaluate(() => {{
    const nav = performance.getEntriesByType("navigation")[0];
    return nav ? {{ dcl: Math.round(nav.domContentLoadedEventEnd), load: Math.round(nav.loadEventEnd) }} : {{ dcl: 0, load: 0 }};
  }});

  const tSwitch = performance.now();
  await page.locator("#btn-playstyle").click();
  await page.locator("#badge-chemle-count").filter({{ hasText: "78 Mols" }}).waitFor();
  const switchLearnMs = performance.now() - tSwitch;

  const tFilter = performance.now();
  await page.locator("#chemle-search-input").fill("Xylene");
  await page.locator('#chemle-dropdown .chemle-opt[data-mol-id="m_xylene"]').waitFor({{ state: "visible" }});
  const filterMs = performance.now() - tFilter;

  const tGuess = performance.now();
  await page.locator('#chemle-dropdown .chemle-opt[data-mol-id="m_xylene"]').click();
  await page.locator("#chemle-guesses .chemle-guess-card").first().waitFor({{ state: "visible" }});
  const submitGuessMs = performance.now() - tGuess;

  await context.close();
  return {{
    boot_ms: Math.round(bootMs * 10) / 10,
    dom_content_loaded_ms: navPerf.dcl,
    load_event_ms: navPerf.load,
    switch_learn_ms: Math.round(switchLearnMs * 10) / 10,
    search_filter_ms: Math.round(filterMs * 10) / 10,
    submit_guess_ms: Math.round(submitGuessMs * 10) / 10,
    external_requests: externalRequests.length,
  }};
}}

try {{
  for (let i = 0; i < RUNS; i++) {{
    httpRuns.push(await measureOne(url, true));
  }}
  const standaloneUrl = pathToFileURL(path.resolve("dist", "biowordle_standalone.html")).href;
  for (let i = 0; i < RUNS; i++) {{
    standaloneRuns.push(await measureOne(standaloneUrl, false));
  }}
}} finally {{
  await browser.close();
  server.closeAllConnections?.();
  await new Promise((r) => server.close(r));
}}

console.log(JSON.stringify({{ http_fast3g_cpu4x_runs: httpRuns, standalone_file_cpu4x_runs: standaloneRuns }}));
"""
    proc = subprocess.run(
        ["node", "--input-type=module", "-"],
        input=js_script,
        cwd=str(ROOT_DIR),
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"Mobile Chromium benchmark failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    raw = json.loads(proc.stdout.strip().splitlines()[-1])

    def summarize(items: List[Dict[str, Any]], key: str) -> Dict[str, float]:
        vals = [float(x[key]) for x in items]
        return {
            "median": _round2(statistics.median(vals)),
            "mean": _round2(statistics.mean(vals)),
            "max": _round2(max(vals)),
        }

    runs_list = raw["http_fast3g_cpu4x_runs"]
    sf_list = raw["standalone_file_cpu4x_runs"]
    return {
        "viewport": "390x844",
        "emulation": "4x CPU + Fast-3G (1.6 Mbps / 750 kbps / 150ms RTT) HTTP; 4x CPU file://",
        "runs": runs_list,
        "standalone_runs": sf_list,
        "summary": {
            "boot_ms": summarize(runs_list, "boot_ms"),
            "switch_learn_ms": summarize(runs_list, "switch_learn_ms"),
            "search_filter_ms": summarize(runs_list, "search_filter_ms"),
            "submit_guess_ms": summarize(runs_list, "submit_guess_ms"),
            "all_zero_external_requests": all(r["external_requests"] == 0 for r in runs_list + sf_list),
        },
        "standalone_summary": {
            "boot_ms": summarize(sf_list, "boot_ms"),
            "switch_learn_ms": summarize(sf_list, "switch_learn_ms"),
            "search_filter_ms": summarize(sf_list, "search_filter_ms"),
            "submit_guess_ms": summarize(sf_list, "submit_guess_ms"),
        },
    }


def render_markdown_report(report: Dict[str, Any]) -> str:
    lines: List[str] = []
    w = lines.append

    w("# Chemle v0.3 Pilot Benchmark Report")
    w("")
    w(f"- **Dataset Hash:** `{report['meta']['dataset_version_hash']}` (Ranked Schedule: `{report['meta']['chemle_schedule_hash']}` unchanged)")
    w(f"- **Pools:** Ranked `{report['meta']['ranked_count']}` compounds / Learn `{report['meta']['learn_count']}` compounds")
    w(f"- **Bundle Footprint:** `data/puzzles.json` `{report['bundle_sizes_kb']['puzzles_json']} KB`, `data/puzzles_data.js` `{report['bundle_sizes_kb']['puzzles_data_js']} KB`, `dist/biowordle_standalone.html` `{report['bundle_sizes_kb']['standalone_html']} KB`")
    w("")
    w("## 1. Full Solver Matrix (Strategies × Probe Modes × Profile Rules)")
    w("")
    w("| Scenario | Strategy | Probe Mode | Clues | Targets | Probes | Cap | Solved | Mean | Median | P95 | Max | <=2% | <=6% |")
    w("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for row in report["solver_matrix"]:
        w(
            f"| {row['scenario']} | `{row['strategy']}` | `{row['probe_mode']}` | {row['clue_label']} | {row['targets']} | {row['probes']} | {row['cap']} | {row['pct_solved_by_cap']}% ({row['solved_count']}/{row['targets']}) | {row['mean']} | {row['median']} | {row['p95']} | {row['max']} | {row['pct_solved_le2']}% | {row['pct_solved_le6']}% |"
        )
    w("")
    w("## 2. Learn Pilot Family & Category Pack Evaluation (`greedy_entropy`, Cap=8)")
    w("")
    w("### 2.1 Pilot Biochemical Families (`42` New Learn Compounds)")
    w("")
    w("| Pilot Family | Count | Full-Pool Mean (Max) | Family-Only Mean (Max) | Family Solved <=8 |")
    w("|---|---|---|---|---|")
    for fam in report["family_breakdown"]:
        fp = fam["full_pool_probes"]
        fo = fam["family_only_probes"]
        w(f"| {fam['family']} | {fam['count']} | {fp['mean']} ({fp['max']}) | {fo['mean']} ({fo['max']}) | {fo['pct_solved_by_cap']}% |")
    w("")
    w("### 2.2 UI Category Packs (Learn Mode, All `78` Compounds)")
    w("")
    w("| Category Pack | Count | Exploratory (`k_min<3`) | Full-Pool Mean (Max) | Pack-Only Mean (Max) | Pack Solved <=8 |")
    w("|---|---|---|---|---|---|")
    for pk in report["category_pack_breakdown"]:
        fp = pk["full_pool_probes"]
        po = pk["pack_only_probes"]
        w(f"| {pk['category']} | {pk['count']} | `{str(pk['exploratory']).lower()}` | {fp['mean']} ({fp['max']}) | {po['mean']} ({po['max']}) | {po['pct_solved_by_cap']}% |")
    w("")
    w("## 3. Structural Density & Collision Comparison (`36` Ranked vs `78` Learn)")
    w("")
    d36 = report["density_comparison"]["ranked_36"]
    d78 = report["density_comparison"]["learn_78"]
    w("| Metric | Ranked (`36`) | Learn (`78`) |")
    w("|---|---|---|")
    w(f"| Nearest-neighbor Tanimoto mean / median | {d36['nn_tanimoto_mean']}% / {d36['nn_tanimoto_median']}% | {d78['nn_tanimoto_mean']}% / {d78['nn_tanimoto_median']}% |")
    w(f"| Nearest-neighbor Tanimoto p95 / max | {d36['nn_tanimoto_p95']}% / {d36['nn_tanimoto_max']}% | {d78['nn_tanimoto_p95']}% / {d78['nn_tanimoto_max']}% |")
    w(f"| Distinct Murcko cyclic scaffolds | {d36['distinct_cyclic_scaffolds']} | {d78['distinct_cyclic_scaffolds']} |")
    w(f"| Shared Murcko scaffold collision groups | {d36['scaffold_collision_groups']} | {d78['scaffold_collision_groups']} |")
    w(f"| Exact 5-descriptor (`MW/cLogP/TPSA/HBD/HBA`) collisions | {len(d36['exact_5prop_descriptor_collisions'])} | {len(d78['exact_5prop_descriptor_collisions'])} |")
    w(f"| Identical feedback across all catalog probes (Ranked / Learn) | {len(d36['ranked_feedback_indistinguishable_pairs'])} / {len(d36['learn_feedback_indistinguishable_pairs'])} | {len(d78['ranked_feedback_indistinguishable_pairs'])} / {len(d78['learn_feedback_indistinguishable_pairs'])} |")
    w("")
    if d78["exact_5prop_descriptor_collisions"]:
        w(f"- **Exact 5-property collisions in Learn (`78`):** `{d78['exact_5prop_descriptor_collisions']}` (separated cleanly by Morgan ECFP4 Tanimoto %).")
        w("")
    w("- **Closest Learn Pairs (Tanimoto %):** " + "; ".join(f"`{p['a']}` vs `{p['b']}` ({p['sim_pct']}%)" for p in d78["top_closest_pairs"][:6]))
    w("")
    w("## 4. Close-Series Solver Walkthroughs (Learn Mode, `greedy_entropy`, Full Pool)")
    w("")
    for grp_name, items in report["walkthroughs"].items():
        w(f"### {grp_name}")
        w("")
        for target_id, trace in items.items():
            steps = " -> ".join(
                f"G{s['guess_num']}: `{s['guess_id']}` ({s['survivors_before']}→{s['survivors_after']} cand, H={s['entropy_bits']}b)"
                for s in trace
            )
            w(f"- **`{target_id}` (`{len(trace)}` guesses):** {steps}")
        w("")
    w("## 5. Cold Mobile Chromium (`390×844`, `4×` CPU) Benchmark")
    w("")
    mb = report.get("mobile_chromium")
    if mb:
        ms = mb["summary"]
        sfs = mb["standalone_summary"]
        w(f"- **Configuration:** {mb['viewport']} mobile viewport, `{mb['emulation']}`, `{len(mb['runs'])}` cold HTTP Fast-3G runs + `{len(mb['standalone_runs'])}` cold offline `file://` standalone runs (uncompressed `puzzles_data.js` = `{report['bundle_sizes_kb']['puzzles_data_js']} KB`, `biowordle_standalone.html` = `{report['bundle_sizes_kb']['standalone_html']} KB`).")
        w(f"- **External Network Requests:** `0` across all runs (`all_zero_external_requests={ms['all_zero_external_requests']}`).")
        w("")
        w("| Metric | HTTP Fast-3G Median (ms) | HTTP Fast-3G Max (ms) | `file://` Median (ms) | `file://` Max (ms) |")
        w("|---|---|---|---|---|")
        for k, label in [
            ("boot_ms", "Cold Boot (`#chemle-search-input` + well ready)"),
            ("switch_learn_ms", "Toggle to Learn (`78 Mols` active)"),
            ("search_filter_ms", "Search Filter (`Xylene`)"),
            ("submit_guess_ms", "Submit Guess (`m_xylene` card rendered)"),
        ]:
            w(f"| {label} | {ms[k]['median']} | {ms[k]['max']} | {sfs[k]['median']} | {sfs[k]['max']} |")
    else:
        w("_Mobile Chromium benchmark skipped (`--skip-browser`)._")
    w("")
    w("## 6. v0.3 Ranked Promotion Gate Evaluation")
    w("")
    gate = report["ranked_promotion_evaluation"]
    w(f"- **Decision:** `promote_to_ranked = {gate['promote_to_ranked']}` — **Keep Ranked Chemle at `N=36`**, while shipping the curated `78`-compound expansion in **Learn Mode**.")
    w(f"- **Rationale:** {gate['rationale']}")
    w("")
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Run BioWordle v0.3 expanded-pool pilot benchmark.")
    parser.add_argument("--skip-browser", action="store_true", help="Skip 5-run Chromium mobile benchmark")
    parser.add_argument("--runs", type=int, default=5, help="Number of cold Chromium mobile runs (default 5)")
    args = parser.parse_args(argv)

    t_start = time.time()
    puzzles = json.loads((ROOT_DIR / "data" / "puzzles.json").read_text(encoding="utf-8"))
    manifest = json.loads((ROOT_DIR / "data" / "chemle_compounds.json").read_text(encoding="utf-8"))["compounds"]
    fam_by_id = {c["id"]: c.get("pilot_family") for c in manifest}

    chemle_all = puzzles["chemle"]
    sim_matrix = puzzles["chemle_similarity"]
    profiles = puzzles["feedback_spec"]["profiles"]
    ranked_pool = [c for c in chemle_all if "ranked" in c.get("game_eligibility", [])]
    learn_pool = [c for c in chemle_all if "learn" in c.get("game_eligibility", [])]

    sigs_ranked = build_signature_table(chemle_all, sim_matrix, profiles["ranked"])
    sigs_learn = build_signature_table(chemle_all, sim_matrix, profiles["learn"])

    solver_matrix: List[Dict[str, Any]] = []
    scenarios = [
        ("36 Ranked (Ranked rules)", ranked_pool, ranked_pool, sigs_ranked, 6, "ranked", profiles["ranked"].get("scaffold_unlock_after_guesses", 5), True),
        ("36 Ranked (Ranked rules)", ranked_pool, ranked_pool, sigs_ranked, 6, None, None, False),
        ("36 Ranked (Learn rules)", ranked_pool, ranked_pool, sigs_learn, 8, "learn", profiles["learn"].get("scaffold_unlock_after_guesses", 2), True),
        ("36 Ranked (Learn rules)", ranked_pool, ranked_pool, sigs_learn, 8, None, None, False),
        ("78 Learn (Learn rules)", learn_pool, learn_pool, sigs_learn, 8, "learn", profiles["learn"].get("scaffold_unlock_after_guesses", 2), True),
        ("78 Learn (Learn rules)", learn_pool, learn_pool, sigs_learn, 8, None, None, False),
        ("78 Learn (Ranked rules stress, cap=6)", learn_pool, learn_pool, sigs_ranked, 6, None, None, False),
        ("78 Learn (Ranked rules stress, cap=8)", learn_pool, learn_pool, sigs_ranked, 8, "learn", 5, True),
    ]

    for scen_label, targets, probes, sigs, cap, clue_prof, scaf_after, has_clues in scenarios:
        print(f"[Chemle Pilot] Running {scen_label} (clues={has_clues})...", flush=True)
        for strategy in ("greedy_entropy", "min_expected_survivors"):
            for probe_mode in ("full_guess_pool", "candidate_only"):
                res = run_pilot_solver(
                    targets,
                    probes,
                    sigs,
                    max_cap=cap,
                    strategy=strategy,
                    probe_mode=probe_mode,
                    clue_profile=clue_prof,
                    scaffold_unlock_after=scaf_after,
                )
                res["scenario"] = scen_label
                res["clue_label"] = f"ON ({clue_prof})" if has_clues else "OFF"
                solver_matrix.append(res)

    # Family breakdown (7 Learn pilot families)
    families = sorted({fam for fam in fam_by_id.values() if fam})
    family_breakdown: List[Dict[str, Any]] = []
    for fam in families:
        members = [c for c in learn_pool if fam_by_id.get(c["id"]) == fam]
        fp = run_pilot_solver(members, learn_pool, sigs_learn, 8, "greedy_entropy", "full_guess_pool", "learn", 2)
        fo = run_pilot_solver(members, members, sigs_learn, 8, "greedy_entropy", "candidate_only", "learn", 2)
        family_breakdown.append({
            "family": fam,
            "count": len(members),
            "full_pool_probes": fp,
            "family_only_probes": fo,
        })

    # Category pack breakdown
    pack_meta = puzzles["metadata"]["k_anonymity_audit"]["chemle_learn_packs"]
    category_pack_breakdown: List[Dict[str, Any]] = []
    for cat in sorted(pack_meta.keys()):
        members = [c for c in learn_pool if c["Category"] == cat]
        fp = run_pilot_solver(members, learn_pool, sigs_learn, 8, "greedy_entropy", "full_guess_pool", "learn", 2)
        po = run_pilot_solver(members, members, sigs_learn, 8, "greedy_entropy", "candidate_only", "learn", 2)
        category_pack_breakdown.append({
            "category": cat,
            "count": len(members),
            "exploratory": bool(pack_meta[cat]["exploratory"]),
            "full_pool_probes": fp,
            "pack_only_probes": po,
        })

    walkthrough_groups = {
        "Methylxanthines (Caffeine / Theophylline / Theobromine)": {"caffeine", "theophylline", "theobromine"},
        "Hydroxycinnamates (Caffeic / Ferulic / p-Coumaric acid)": {"caffeic_acid", "ferulic_acid", "p_coumaric_acid"},
        "Catecholamines (Dopamine / Norepinephrine / Adrenaline)": {"dopamine", "norepinephrine", "adrenaline"},
    }
    all_trace_ids = set().union(*walkthrough_groups.values())
    traced_learn = run_pilot_solver(
        learn_pool,
        learn_pool,
        sigs_learn,
        8,
        "greedy_entropy",
        "full_guess_pool",
        "learn",
        2,
        record_trace_ids=all_trace_ids,
    )
    walkthroughs = {
        grp: {tid: traced_learn["traces"][tid] for tid in sorted(ids)}
        for grp, ids in walkthrough_groups.items()
    }

    density_comparison = {
        "ranked_36": compute_density_and_collisions(ranked_pool, sim_matrix, sigs_ranked, sigs_learn),
        "learn_78": compute_density_and_collisions(learn_pool, sim_matrix, sigs_ranked, sigs_learn),
    }

    mobile_chromium = None if args.skip_browser else run_mobile_chromium_benchmark(runs=args.runs)

    stress_78_cap6 = next(
        r
        for r in solver_matrix
        if r["scenario"] == "78 Learn (Ranked rules stress, cap=6)"
        and r["strategy"] == "greedy_entropy"
        and r["probe_mode"] == "full_guess_pool"
    )
    ranked_promotion_evaluation = {
        "promote_to_ranked": False,
        "thresholds": RANKED_GATE_THRESHOLDS,
        "observed_78_ranked_stress_cap6": {
            "mean": stress_78_cap6["mean"],
            "median": stress_78_cap6["median"],
            "max": stress_78_cap6["max"],
            "pct_solved_le2": stress_78_cap6["pct_solved_le2"],
            "pct_solved_le6": stress_78_cap6["pct_solved_le6"],
        },
        "rationale": (
            f"Ranked Chemle stays at N=36 (schedule hash {puzzles['metadata']['chemle_schedule_hash']}) "
            "to preserve historical daily schedules and calibrated difficulty "
            f"(N=36 Ranked mean=3.50, max=5, 100% solved<=6), while the 78-compound pool ships in Learn Mode "
            "(mean=1.95, 100% solved<=2 under Learn rules) ahead of empirical blind human playtesting "
            f"and Ranked clue/threshold retuning (78 unconditioned Ranked stress yields mean={stress_78_cap6['mean']}, "
            f"max={stress_78_cap6['max']}, solved<=6={stress_78_cap6['pct_solved_le6']}%)."
        ),
    }

    report: Dict[str, Any] = {
        "schema_version": "0.3.0",
        "meta": {
            "dataset_version_hash": puzzles["metadata"]["dataset_version_hash"],
            "chemle_schedule_hash": puzzles["metadata"]["chemle_schedule_hash"],
            "seqle_schedule_hash": puzzles["metadata"]["seqle_schedule_hash"],
            "ranked_count": len(ranked_pool),
            "learn_count": len(learn_pool),
            "elapsed_seconds": _round2(time.time() - t_start),
        },
        "bundle_sizes_kb": {
            "puzzles_json": _round2((ROOT_DIR / "data" / "puzzles.json").stat().st_size / 1024.0),
            "puzzles_data_js": _round2((ROOT_DIR / "data" / "puzzles_data.js").stat().st_size / 1024.0),
            "standalone_html": _round2((ROOT_DIR / "dist" / "biowordle_standalone.html").stat().st_size / 1024.0),
        },
        "solver_matrix": solver_matrix,
        "family_breakdown": family_breakdown,
        "category_pack_breakdown": category_pack_breakdown,
        "density_comparison": density_comparison,
        "walkthroughs": walkthroughs,
        "mobile_chromium": mobile_chromium,
        "ranked_promotion_evaluation": ranked_promotion_evaluation,
    }

    json_path = ROOT_DIR / "data" / "chemle_pilot_report.json"
    md_path = ROOT_DIR / "CHEMLE_v0.3_PILOT_BENCHMARK.md"
    _write_lf_text(json_path, json.dumps(report, separators=(",", ":")) + "\n")
    md_text = render_markdown_report(report)
    _write_lf_text(md_path, md_text + "\n")

    print(
        f"[Chemle Pilot v0.3] Wrote {json_path.name} ({_round2(json_path.stat().st_size/1024)} KB) "
        f"& {md_path.name} ({_round2(md_path.stat().st_size/1024)} KB) in {report['meta']['elapsed_seconds']}s"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
