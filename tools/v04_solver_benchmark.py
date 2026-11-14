#!/usr/bin/env python3
"""Experimental BioWordle v0.4 Learn/Challenge solver comparison.

This does not write production data. It calls the real compiler with an in-memory
manifest and records, rather than suppresses, every compiler difficulty-gate
failure. It never changes eligibility or reference structures. Use with either a
42-entry additions manifest or a complete 120-entry manifest in {"compounds": [...]}
format. See --help. Probe tie-breaking and Learn hint policy are reported in JSON.
"""
import argparse
import collections
import contextlib
import copy
import functools
import hashlib
import io
import json
import math
import statistics
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import build_puzzles as builder
import chem_engine as engine


def read_json(path):
    raw = Path(path).read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def ladder_text(record, profile, level):
    return next((x['text'] for x in record['clue_ladder'][profile] if x['level'] == level), None)


def solve(data, profile, strategy):
    spec = data['feedback_spec']['profiles'][profile]
    ans = sorted(r['id'] for r in data['chemle'] if 'learn' in r.get('game_eligibility', []))
    probes = sorted(r['id'] for r in data['chemle'] if 'learn' in r.get('guess_eligibility', []))
    by = {r['id']: r for r in data['chemle']}
    if not set(ans) <= set(probes):
        raise ValueError('answer set must be a subset of legal probes')
    sim = data['chemle_similarity']
    sig = {p: {a: engine.evaluate_chemle_feedback_from_spec(by[p], by[a], sim[p][a], spec)['signature'] for a in ans} for p in probes}
    opening = {a: ladder_text(by[a], profile, 0) for a in ans}
    groups = collections.Counter(opening.values())
    cap = int(spec['max_guesses'])
    assert cap == 8, f'Expected 8 guesses; got {cap}'
    def partition(p, candidates):
        out = collections.defaultdict(list)
        for a in candidates:
            out[sig[p][a]].append(a)
        return out
    @functools.lru_cache(None)
    def choose(candidates, used):
        if len(candidates) == 1:
            return candidates[0]
        count = len(candidates)
        cset = set(candidates)
        scores = []
        for p in probes:
            if p in used:
                continue
            parts = partition(p, candidates)
            counts = [len(x) for x in parts.values()]
            entropy = -sum(n/count*math.log2(n/count) for n in counts)
            expected = sum(n*n for n in counts)/count
            scores.append((p, entropy, expected, max(counts), parts))
        if not scores:
            raise ValueError('no legal unused probes')
        if strategy == 'entropy':
            return min(scores, key=lambda x: (-round(x[1], 10), x[0] not in cset, x[0]))[0]
        if strategy == 'expected':
            return min(scores, key=lambda x: (round(x[2], 10), x[0] not in cset, x[0]))[0]
        if strategy == 'minimax':
            return min(scores, key=lambda x: (x[3], round(x[2], 10), x[0] not in cset, x[0]))[0]
        finalists = sorted(scores, key=lambda x: (x[2], x[0] not in cset, x[0]))[:8]
        def second(x):
            value = 0.0
            for key, branch in x[4].items():
                if key == 'WIN':
                    continue
                m = len(branch)
                legal = [h for h in probes if h not in used and h != x[0]]
                value += min(sum(len(v)**2 for v in partition(h, branch).values())/m for h in legal) * m/count
            return value
        return min(finalists, key=lambda x: (round(second(x), 10), x[0] not in cset, x[0]))[0]
    def scaffold(r):
        # UI shows an acyclic/no-scaffold message separately; exact scaffold SMILES
        # is only a proxy for what a player can infer from the rendered SVG.
        return (bool(r.get('has_scaffold')), r.get('scaffold_smiles') or '')
    results = {}
    for target in ans:
        candidates = tuple(a for a in ans if opening[a] == opening[target])
        history = []
        observations = []
        for step in range(1, cap+1):
            p = choose(candidates, tuple(sorted(history)))
            history.append(p)
            before = len(candidates)
            feedback = sig[p][target]
            if p == target:
                observations.append({'guess':p, 'feedback':'WIN','candidates_before':before, 'candidates_after':1})
                break
            candidates = tuple(a for a in candidates if sig[p][a] == feedback)
            # Learn auto unlocks actual ladder levels at 1/2/3 guesses; no manual
            # hints are pressed. In Learn the scaffold is visible at level 1.
            # Challenge gives only actual ladder clues plus scaffold after step 5.
            for clue in by[target]['clue_ladder'][profile]:
                if clue['after_guesses'] <= step:
                    candidates = tuple(a for a in candidates if ladder_text(by[a], profile, clue['level']) == clue['text'])
            scaffold_at = 1 if profile == 'learn' else int(spec.get('scaffold_unlock_after_guesses', 5))
            if step >= scaffold_at:
                candidates = tuple(a for a in candidates if scaffold(by[a]) == scaffold(by[target]))
            if target not in candidates:
                raise AssertionError(f'target eliminated: {target}')
            observations.append({'guess':p, 'feedback':feedback,'candidates_before':before, 'candidates_after':len(candidates)})
        won = history[-1] == target
        results[target] = {'solved':won,'guesses':len(history) if won else None,'history':history,
                           'opening_k':groups[opening[target]],'trace':observations}
    counts = collections.Counter(str(x['guesses']) if x['solved'] else 'unsolved' for x in results.values())
    solved = [x['guesses'] for x in results.values() if x['solved']]
    all_probe = []
    for clue in sorted(groups):
        group = [a for a in ans if opening[a] == clue]
        for p in probes:
            values = collections.Counter(sig[p][a] for a in group)
            singleton = sum(a != p and values[sig[p][a]] == 1 for a in group)
            all_probe.append({'opening_clue':clue,'k':len(group),'probe':p,
                              'nonwinning_singletons':singleton,'singleton_fraction':singleton/len(group),
                              'expected_survivors':sum(n*n for n in values.values())/len(group)})
    all_probe.sort(key=lambda x:(-x['singleton_fraction'],-x['nonwinning_singletons'],x['probe']))
    return {'profile':profile,'strategy':strategy,'targets':len(ans),'legal_probes':len(probes),
            'attempt_cap':cap,'opening_groups':dict(sorted(groups.items())), 'opening_k_min':min(groups.values()),
            'distribution':dict(sorted(counts.items())), 'unsolved_targets':sorted(k for k,v in results.items() if not v['solved']),
            'success_within_6':sum(v<=6 for v in solved)/len(ans), 'success_within_8':len(solved)/len(ans),
            'mean_solved':statistics.mean(solved) if solved else None,
            'median_solved':statistics.median(solved) if solved else None,
            'p95_solved':sorted(solved)[math.ceil(.95*len(solved))-1] if solved else None,
            'opening_probe_worst_10':all_probe[:10], 'per_target':results}


def ranked_comparison(base, expansion):
    a = {r['id']:r for r in base['chemle'] if 'ranked' in r.get('game_eligibility',[])}
    b = {r['id']:r for r in expansion['chemle'] if 'ranked' in r.get('game_eligibility',[])}
    ap = {r['id'] for r in base['chemle'] if 'ranked' in r.get('guess_eligibility',[])}
    bp = {r['id'] for r in expansion['chemle'] if 'ranked' in r.get('guess_eligibility',[])}
    fields = ['inchikey','canonical_smiles','MW','cLogP','TPSA','HBD','HBA','Rings','AromaticRings','HeavyAtoms','functional_groups','Category','super_class','scaffold_smiles','has_scaffold']
    changed = []
    for k in sorted(a.keys() & b.keys()):
        for f in fields:
            if a[k].get(f) != b[k].get(f): changed.append(f'{k}.{f}')
        if a[k]['clue_ladder']['ranked'] != b[k]['clue_ladder']['ranked']: changed.append(f'{k}.ranked_clues')
    diffs = []
    if ap == bp and a.keys() == b.keys():
        for p in sorted(ap):
            for t in sorted(a):
                old = engine.evaluate_chemle_feedback_from_spec(a[p],a[t],base['chemle_similarity'][p][t],base['feedback_spec']['profiles']['ranked'])['signature']
                new = engine.evaluate_chemle_feedback_from_spec(b[p],b[t],expansion['chemle_similarity'][p][t],expansion['feedback_spec']['profiles']['ranked'])['signature']
                if old != new: diffs.append([p,t,old,new])
    return {'answer_roster_equal':set(a)==set(b),'probe_roster_equal':ap==bp,
            'ranked_descriptor_or_clue_changes':changed,'ranked_feedback_differences':diffs,
            'schedule_hashes':[base['metadata']['chemle_schedule_hash'],expansion['metadata']['chemle_schedule_hash']],
            'rules_hashes':[base['metadata']['chemle_rules_hash'],expansion['metadata']['chemle_rules_hash']],
            'note':'Ranked session storage key includes the global Chemle rules hash; changed hash archives previous sessions even when Ranked feedback and schedule are identical.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', help='42 additions or full 120, JSON object with compounds array; omit for baseline only')
    parser.add_argument('--expected-additions', type=int, default=42, help='required unique additions; use 0 for a baseline check')
    parser.add_argument('--output', default='tools/v04_solver_results.json')
    parser.add_argument('--compiled-output', help='optional experimental compiler JSON output; not a production path')
    parser.add_argument('--research-preview', action='store_true', help='temporarily enable pending candidate rows in memory for research; never marks them scientifically verified in the source manifest')
    parser.add_argument('--challenge-mw', type=float, help='experimental challenge MW exact tolerance')
    parser.add_argument('--challenge-logp', type=float, help='experimental challenge cLogP exact tolerance')
    parser.add_argument('--challenge-edges', nargs=2, type=int, metavar=('LOW','HIGH'), help='experimental challenge similarity bucket edges')
    parser.add_argument('--challenge-super-after', type=int, help='experimental challenge superclass unlock')
    parser.add_argument('--challenge-scaffold-after', type=int, help='experimental challenge scaffold unlock')
    args = parser.parse_args()
    production = {str((ROOT/x).resolve()) for x in ['data/puzzles.json','data/puzzles_data.js','dist/biowordle_standalone.html','data/chemle_compounds.json']}
    for out in [args.output, args.compiled_output]:
        if out and str(Path(out).resolve()) in production:
            parser.error('output must not overwrite production data')
    base, base_sha = read_json(ROOT/'data/puzzles.json')
    spec, spec_sha = read_json(ROOT/'data/feedback_spec.json')
    source, source_sha = read_json(ROOT/'data/chemle_compounds.json')
    spec = copy.deepcopy(spec)
    challenge = spec['profiles']['challenge']
    changes = {}
    for name, value, index in [('MW',args.challenge_mw,0),('cLogP',args.challenge_logp,1)]:
        if value is not None:
            if value <= 0: parser.error(f'{name} tolerance must be positive')
            challenge['chemle_fields'][index]['exact_tol'] = value
            changes[name + '_exact_tol'] = value
    if args.challenge_edges is not None:
        if not 0 <= args.challenge_edges[0] < args.challenge_edges[1] <= 100: parser.error('challenge edges must satisfy 0 <= LOW < HIGH <= 100')
        challenge['tanimoto_bucket']['edges'] = args.challenge_edges
        changes['tanimoto_edges'] = args.challenge_edges
    for flag, field in [(args.challenge_super_after,'super_class_unlock_after_guesses'),(args.challenge_scaffold_after,'scaffold_unlock_after_guesses')]:
        if flag is not None:
            if not 1 <= flag <= 8: parser.error('clue timing must be in 1..8')
            challenge[field] = flag
            changes[field] = flag
    baseline = source['compounds']
    rows = baseline
    manifest_sha = None
    if args.manifest:
        manifest, manifest_sha = read_json(args.manifest)
        rows = manifest['compounds'] if isinstance(manifest,dict) else manifest
        base_by = {r['id']:r for r in baseline}
        incoming = {r['id']:r for r in rows}
        if len(incoming) != len(rows): raise ValueError('duplicate IDs in candidate manifest')
        if set(base_by) <= set(incoming):
            if any(incoming[k] != v for k,v in base_by.items()):
                raise ValueError('full manifest must preserve all 78 base records byte-for-byte by value')
        elif set(base_by).isdisjoint(incoming):
            rows = baseline + rows
        else:
            raise ValueError('manifest must be additions only or contain all base records unchanged')
        added = [r for r in rows if r['id'] not in base_by]
        if len(added) != args.expected_additions:
            raise ValueError(f'expected {args.expected_additions} additions, got {len(added)}')
        if args.research_preview:
            rows = copy.deepcopy(rows)
            for r in rows:
                if r['id'] not in base_by:
                    if r.get('review_status') not in ('verified', 'pending', 'draft'):
                        raise ValueError(f"{r['id']}: rejected or unknown status cannot be previewed")
                    r['review_status'] = 'verified'  # in-memory validator adapter only; NOT a scientific review
                    r['game_eligibility'] = ['learn']
                    r['guess_eligibility'] = ['learn']
        if any(set(r.get('game_eligibility',[])) != {'learn'} or set(r.get('guess_eligibility',[])) != {'learn'} for r in rows if r['id'] not in base_by):
            raise ValueError('new records must be Learn-only answers and probes (or explicitly use --research-preview)')
    elif args.expected_additions != 0:
        parser.error('provide --manifest or use --expected-additions 0 for a baseline-only check')
    # Capture failed gates, retaining their exact assertions for the report. Every
    # schema, identity, duplicate and other compiler check runs unmodified.
    gate_failures = []
    old_gate = builder.assert_solver_gate
    def audited_gate(stats, gate, profile):
        try: old_gate(stats, gate, profile)
        except AssertionError as exc: gate_failures.append({'profile':profile,'message':str(exc)})
    builder.assert_solver_gate = audited_gate
    log = io.StringIO()
    t0=time.perf_counter()
    try:
        with contextlib.redirect_stdout(log):
            expanded = builder.build_all(copy.deepcopy(spec), rows)
    finally:
        builder.assert_solver_gate = old_gate
    elapsed = time.perf_counter()-t0
    if args.compiled_output:
        Path(args.compiled_output).write_text(json.dumps(expanded,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    result = {'base_puzzles_sha256':base_sha,'base_manifest_sha256':source_sha,'spec_sha256':spec_sha,
              'candidate_manifest_sha256':manifest_sha,'compile_seconds':elapsed,'experimental_challenge_overrides':changes,
              'research_preview':args.research_preview,
              'research_preview_warning':'Candidate review status and Learn eligibility are temporarily adapted in memory only. This is not a scientific verification or permission to ship.' if args.research_preview else None,
              'compiler_gate_failures':gate_failures,
              'gate_capture_note':'Each profile reports its first failing compiler assertion; numerical statistics remain available in compiled output. Existing release gates are not changed.',
              'compiler_log':log.getvalue().replace('Solver Gate PASS','Solver Gate EVALUATED (see compiler_gate_failures)'),
              'ranked_invariance':ranked_comparison(base,expanded),
              'simulation_policy':{'target_prior':'uniform','learn_hints':'automatic ladder only; no manual hint presses',
               'clue_model':'exact published ladder text; scaffold smiles proxy for visible scaffold',
               'challenge':'v0.3.2 spec and clue timings unless explicitly overridden; no narrative or manual hints',
               'lookahead':'two-ply, best eight first probes by expected survivors; future clue reveals ignored during optimization',
               'tie_break':'prefer possible answer then lexically first ID; no duplicate guesses; differs from production solver tie-breaking',
               'limitations':'perfect knowledge of exact clue text and scaffold proxy; not human difficulty or experimental validation'},
              'scenarios':{}}
    for label, data in [('base_78',base),('expanded',expanded)]:
        result['scenarios'][label]=[]
        for profile in ('learn','challenge'):
            for strategy in ('entropy','expected','minimax','lookahead'):
                outcome=solve(data,profile,strategy)
                result['scenarios'][label].append(outcome)
                print(label,profile,strategy,outcome['distribution'],'k_min',outcome['opening_k_min'],flush=True)
    Path(args.output).write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')

if __name__ == '__main__':
    main()
