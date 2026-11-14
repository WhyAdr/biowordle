#!/usr/bin/env python3
"""Deterministic, exhaustive-target Chemle challenge simulations; no production changes.

Example: python tools/challenge_benchmark.py --base 83f05f96102abea5d2fcce0ed75dbf3535127f71 --output challenge_results.json
"""
import argparse, collections, copy, functools, hashlib, json, math, subprocess, sys, types
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def load(path, rev):
    raw = subprocess.check_output(['git', 'show', f'{rev}:{path}'], cwd=ROOT) if rev else (ROOT/path).read_bytes()
    return json.loads(raw), hashlib.sha256(raw).hexdigest()


def run(data, original, name, pool, strategy, progression, mw, logp, edges):
    records = {r['id']: r for r in data['chemle'] if pool in r.get('game_eligibility', [])}
    ids = tuple(sorted(records)); n = len(ids)
    if len({records[t]['inchikey'] for t in ids}) != n:
        raise ValueError('Benchmark assumes unique winning InChIKeys; update solver for duplicate identities')
    spec = copy.deepcopy(original)
    spec['chemle_fields'][0]['exact_tol'] = mw
    spec['chemle_fields'][1]['exact_tol'] = logp
    spec['tanimoto_bucket']['edges'] = edges
    sim = data.get('chemle_similarity') or data['similarity_matrix']
    table = {g: {t: evaluate_chemle_feedback_from_spec(records[g], records[t], sim[g][t], spec)['signature'] for t in ids} for g in ids}
    # The existing Learn level-0 size/polarity band covers the whole 78 pool.
    # Ranked baseline uses its own published level-0 band.
    clue_profile = 'ranked' if pool == 'ranked' else 'learn'
    def opening(t):
        ladder = records[t]['clue_ladder'][clue_profile]
        return next(e['text'] for e in ladder if e['level'] == 0)
    groups = collections.Counter(opening(t) for t in ids)
    def progressive(cands, target, step):
        if pool == 'ranked':
            # Production Ranked reveals its full narrative ladder, including
            # potentially identifying text, at each configured unlock.
            target_ladder = records[target]['clue_ladder']['ranked']
            for clue in target_ladder:
                if int(clue['after_guesses']) <= step:
                    level, text = clue['level'], clue['text']
                    cands = tuple(c for c in cands if any(e['level'] == level and e['text'] == text for e in records[c]['clue_ladder']['ranked']))
        # Coarse pre-existing classes and scaffold are the modeled reveal only;
        # target-specific narrative clues are deliberately not assumed equivalent.
        sc_step, scaffold_step = progression
        if step >= sc_step:
            cands = tuple(c for c in cands if records[c].get('super_class') == records[target].get('super_class'))
        if step >= scaffold_step:
            cands = tuple(c for c in cands if (records[c].get('scaffold_smiles') or '') == (records[target].get('scaffold_smiles') or ''))
        return cands
    def parts(g, cands):
        d = collections.defaultdict(list)
        for t in cands: d[table[g][t]].append(t)
        return d
    @functools.lru_cache(None)
    def choose(cands):
        if len(cands) == 1: return cands[0]
        count = len(cands); cset = set(cands); scores=[]
        for g in ids:
            p = parts(g,cands); sizes=[len(v) for v in p.values()]
            entropy=-sum((k/count)*math.log2(k/count) for k in sizes)
            expected=sum(k*k for k in sizes)/count
            scores.append((g,entropy,expected,max(sizes),p))
        # Stable tie breaking prefers a possible answer, then lexical ID.
        if strategy == 'entropy': return min(scores,key=lambda x:(-round(x[1],10),x[0] not in cset,x[0]))[0]
        if strategy == 'expected': return min(scores,key=lambda x:(round(x[2],10),x[0] not in cset,x[0]))[0]
        if strategy == 'minimax': return min(scores,key=lambda x:(x[3],round(x[2],10),x[0] not in cset,x[0]))[0]
        # Bounded 2-ply: 8 best expected-survivor first probes, each branch
        # evaluated against every legal second probe. Ignores future clue reveals.
        finalists=sorted(scores,key=lambda x:(x[2],x[0] not in cset,x[0]))[:8]
        def second(x):
            total=0
            for sig, branch in x[4].items():
                if sig == 'WIN': continue
                b=tuple(branch)
                total += min(sum(len(v)**2 for v in parts(h,b).values())/len(b) for h in ids)*len(b)/count
            return total
        return min(finalists,key=lambda x:(round(second(x),10),x[0] not in cset,x[0]))[0]
    outcomes={}
    for target in ids:
        cands=tuple(t for t in ids if opening(t)==opening(target)); history=[]
        for step in range(1,9):
            guess=choose(cands); history.append(guess)
            if guess==target: break
            cands=tuple(t for t in cands if table[guess][t]==table[guess][target])
            cands=progressive(cands,target,step)
            assert target in cands
        outcomes[target]={'guesses':len(history) if history[-1]==target else None,'history':history,'opening_k':groups[opening(target)]}
    dist=collections.Counter(str(v['guesses']) if v['guesses'] is not None else 'unsolved' for v in outcomes.values())
    probes=[]
    # Opening probes are scored inside each actual clue group; singleton coverage
    # excludes the probe itself, whose WIN is necessarily identifying.
    for clue in sorted(groups):
        cs=tuple(t for t in ids if opening(t)==clue)
        for g in ids:
            counts=collections.Counter(table[g][t] for t in cs)
            singleton=sum(1 for t in cs if t!=g and counts[table[g][t]]==1)
            probes.append({'clue':clue,'k':len(cs),'probe':g,'singleton_nonwin':singleton,'singleton_nonwin_fraction':round(singleton/len(cs),4),'expected_survivors':round(sum(x*x for x in counts.values())/len(cs),3)})
    over=sorted(probes,key=lambda x:(-x['singleton_nonwin_fraction'],-x['singleton_nonwin'],x['probe']))[:10]
    vals=[v['guesses'] for v in outcomes.values() if v['guesses'] is not None]
    # Nonuniform / novice prior sensitivity over the active target pool
    novice_weights = {}
    heavy_weights = {}
    for tid in ids:
        r = records[tid]
        # Novice prior upweights ranked/familiar clinical/core metabolites (3x) vs specialized natural products (1x).
        novice_weights[tid] = 3.0 if 'ranked' in r.get('game_eligibility', []) else 1.0
        # Heavy/lipophilic prior upweights larger scaffolds (HeavyAtoms >= 18 -> 3x).
        heavy_weights[tid] = 3.0 if int(r.get('HeavyAtoms', 0)) >= 18 else 1.0
    def _weighted_summary(wmap):
        tot = sum(wmap[t] for t in ids)
        w_le6 = sum(wmap[t] for t, v in outcomes.items() if v['guesses'] is not None and v['guesses'] <= 6) / tot
        w_le8 = sum(wmap[t] for t, v in outcomes.items() if v['guesses'] is not None) / tot
        w_mean = sum(wmap[t] * v['guesses'] for t, v in outcomes.items() if v['guesses'] is not None) / sum(wmap[t] for t, v in outcomes.items() if v['guesses'] is not None)
        return {'weighted_mean': round(w_mean, 4), 'weighted_le6': round(w_le6, 4), 'weighted_le8': round(w_le8, 4)}
    scaffold_mismatches = [
        t for t in ids
        if records[t].get('has_scaffold') and not (records[t].get('scaffold_smiles') and records[t].get('svg_scaffold'))
    ]
    return {
        'name': name,
        'pool': pool,
        'targets': n,
        'strategy': strategy,
        'rules': {'mw_exact_tol': mw, 'logp_exact_tol': logp, 'tanimoto_edges': edges, 'superclass_after': progression[0], 'scaffold_after': progression[1]},
        'opening_groups': dict(sorted(groups.items())),
        'opening_k_min': min(groups.values()),
        'opening_k_distribution': dict(sorted(collections.Counter(groups[opening(t)] for t in ids).items())),
        'distribution': dict(sorted(dist.items())),
        'success_le6': sum(v <= 6 for v in vals) / n,
        'success_le8': len(vals) / n,
        'mean_solved': sum(vals) / len(vals) if vals else None,
        'prior_sensitivity': {
            'novice_familiarity_3x': _weighted_summary(novice_weights),
            'heavy_scaffold_ge18_3x': _weighted_summary(heavy_weights),
        },
        'scaffold_svg_smiles_mismatches': len(scaffold_mismatches),
        'overinformative_probes': over,
        'per_target': outcomes,
    }


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--base',help='read immutable tracked input files at git revision'); parser.add_argument('--output',default='challenge_benchmark_results.json'); args=parser.parse_args()
    data, dh=load('data/puzzles.json',args.base); spec, sh=load('data/feedback_spec.json',args.base)
    engine_raw = subprocess.check_output(['git', 'show', f'{args.base}:chem_engine.py'], cwd=ROOT) if args.base else (ROOT/'chem_engine.py').read_bytes()
    engine = types.ModuleType('challenge_engine_snapshot')
    exec(compile(engine_raw, 'chem_engine.py', 'exec'), engine.__dict__)
    global evaluate_chemle_feedback_from_spec
    evaluate_chemle_feedback_from_spec = engine.evaluate_chemle_feedback_from_spec
    scenarios=[('v03_ranked_baseline','ranked',(4,5),180,4.2,[72,88]),('expanded_ranked_like','learn',(4,5),180,4.2,[72,88]),('expanded_earlier_clue','learn',(3,5),180,4.2,[72,88]),('expanded_tighter_feedback','learn',(4,5),140,3.5,[60,80])]
    results=[]
    for name,pool,progress,mw,logp,edges in scenarios:
        for strategy in ('entropy','expected','minimax','lookahead'):
            result=run(data,spec['profiles']['ranked'],name,pool,strategy,progress,mw,logp,edges); results.append(result)
            print(name,strategy,result['distribution'],'<=6',round(result['success_le6']*100,1),'<=8',round(result['success_le8']*100,1),'kmin',result['opening_k_min'],flush=True)
    out={'input_revision':args.base or 'working tree','puzzles_sha256':dh,'spec_sha256':sh,'engine_sha256':hashlib.sha256(engine_raw).hexdigest(),'uniform_target_prior':True,'results':results}
    with open(args.output, 'w', encoding='utf-8', newline='\n') as f:
        f.write(json.dumps(out, indent=2, ensure_ascii=False) + '\n')
if __name__=='__main__': main()
