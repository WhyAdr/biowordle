# BioWordle v0.3.2 — Challenge Practice design and calibration

## Decision and scope

Keep the 36-compound Ranked daily roster and streaks unchanged. Prototype an explicitly experimental, repeatable **Challenge Practice** using all 78 curated Chemle compounds as answers and guesses. It is a separate local practice session, never a Ranked result, daily completion, classroom score, or streak event. Default candidate: eight attempts, Ranked-style coarse feedback, opening size/polarity band, broad superclass after guess 4, and scaffold after guess 5. Reveal neither target-specific narrative nor an exact percentage during play. Explain that feedback derives from a defined molecular representation, not a complete chemical identification.

This document and the solver are a **simulation prototype**, not an interactive implementation or evidence of human difficulty. Do not promote the expanded pool to Ranked on these numbers. The base benchmark below uses the uncorrected base data; the corrected-data rerun and its consequences are recorded separately below. Scope of this experiment is Chemle only.

## Reproduce and interpret

From the repository root with its existing Python environment:

```powershell
.\.venv\Scripts\python.exe tools/challenge_benchmark.py --base 83f05f96102abea5d2fcce0ed75dbf3535127f71 --output tools/challenge_benchmark_base_results.json
```

Omit `--base` to test working-tree data; write to a different output filename to preserve this comparison. The output records SHA-256 of both data inputs and the evaluator source, settings, each target's entire sequence of guesses and its opening candidate count, and the 10 most identifying opening probe/group pairs. The script uses the Python feedback evaluator from the specified revision (or working tree when `--base` is omitted) and precomputed similarity matrix; it does not independently validate scientific structures or browser parity. It exhaustively simulates each answer with a uniform target prior and all catalog guesses available as probes. Ties prefer a remaining candidate and then lexical ID. Every run has eight opportunities so <=6 and <=8 are measured on the same trajectories. A real six-attempt session would terminate after six; late hints cannot retroactively improve <=6.

For expanded scenarios, the opening clue is the existing Learn level-0 band over all 78 entries. Only the superclass and scaffold are modeled later; any uniquely identifying narrative is intentionally omitted. The v0.3 Ranked baseline uses the Ranked opening clue and actual published narrative ladder at its configured turns, plus the existing superclass and scaffold. Consequently, this comparison changes both pool and clues: it does **not** isolate the effect of pool size. Scaffold equivalence uses `scaffold_smiles`; implementation must confirm that it matches the visible SVG clue (especially blank/acyclic scaffolds). Missing clues, visual differences, synonym cues, and imperfect human recognition can change the true information content.

Strategies: Shannon entropy; minimum expected survivors; minimax maximum remaining branch, breaking ties by expected survivors; and bounded two-ply expected survivors (top eight first probes by expected survivors, then all 78/36 legal second probes). Lookahead ignores future clue reveals while selecting its probe and is not an optimal solver or upper bound. Across all methods candidates are conditioned on revealed clues and feedback, and the output includes missed targets explicitly. Mean is conditional on solved targets, but all targets solved within eight here. This is a bot with perfect catalog knowledge, not a human success estimate.

## Observed results on base commit

Distributions list the count solved in each number of guesses. All denominators are 36 or 78 targets, respectively.

| Configuration | Strategy | Guess distribution | Mean | <=6 | <=8 |
|---|---|---|---:|---:|---:|
| v0.3 Ranked | entropy | 2:5, 3:15, 4:9, 5:7 | 3.50 | 36/36 | 36/36 |
| v0.3 Ranked | expected / minimax | 2:4, 3:16, 4:9, 5:7 | 3.53 | 36/36 | 36/36 |
| v0.3 Ranked | bounded lookahead | 1:1, 2:2, 3:14, 4:10, 5:9 | 3.67 | 36/36 | 36/36 |
| Expanded Ranked-like | entropy / expected / minimax | 1:1, 2:12, 3:29, 4:20, 5:12, 6:3, 7:1 | 3.55 | 77/78 | 78/78 |
| Expanded Ranked-like | bounded lookahead | 1:1, 2:11, 3:23, 4:24, 5:15, 6:3, 7:1 | 3.69 | 77/78 | 78/78 |
| Expanded, superclass after 3 | entropy / expected / minimax | 1:1, 2:12, 3:29, 4:28, 5:6, 6:1, 7:1 | 3.42 | 77/78 | 78/78 |
| Expanded, superclass after 3 | bounded lookahead | 1:1, 2:11, 3:23, 4:31, 5:9, 6:2, 7:1 | 3.59 | 77/78 | 78/78 |
| Expanded, tighter feedback | entropy | 1:1, 2:14, 3:28, 4:21, 5:13, 6:1 | 3.44 | 78/78 | 78/78 |
| Expanded, tighter feedback | expected | 1:2, 2:12, 3:29, 4:21, 5:13, 6:1 | 3.44 | 78/78 | 78/78 |
| Expanded, tighter feedback | minimax | 1:2, 2:9, 3:32, 4:20, 5:13, 6:2 | 3.50 | 78/78 | 78/78 |
| Expanded, tighter feedback | bounded lookahead | 1:1, 2:10, 3:29, 4:23, 5:13, 6:2 | 3.55 | 78/78 | 78/78 |

Ranked-like uses current Ranked thresholds (MW exact tolerance 180 Da, cLogP 4.2, Tanimoto edges 72/88) and superclass/scaffold after 4/5 guesses. Earlier-clue changes superclass to after 3 only. Tighter feedback keeps 4/5 but tests 140 Da, 3.5, and Tanimoto edges 60/80. These are coarse classifications, not experimentally measured tolerances. They were a small, predeclared three-configuration sensitivity comparison; they do not prove global optimization. The seven-guess target in each of the first two expanded scenarios is phenol for the entropy bot. The tighter configuration solves all within six in this sample, but increases information and does not establish a better human experience; recommend current Ranked thresholds as the first playtest default, retaining tighter feedback as a controlled alternative.

## Opening anonymity and probe leakage

The Ranked opening groups have sizes 10, 13 and 13. Expanded opening groups have sizes 6, 11, 13, 13, 16 and 19: minimum opening anonymity is **6**, and the six targets in the smallest group comprise 7.7% of the catalog. This is anonymity within the defined answer pool, not privacy or a guarantee that players lack outside knowledge. Each answer appears in exactly one opening group.

For the Ranked-like expanded scenario, the most identifying opening pairs in this simulation identify two *non-winning* targets out of an 11-member group from a single response (18.2%); examples include ascorbic acid, caffeine, and ferulic acid as probes in that group. Ascorbic acid leaves 5.0 expected candidates under the uniform prior; caffeine and ferulic acid each leave 7.545. The directly guessed answer is excluded from singleton counts. Exact group text and all top pairs are in the JSON; this is an ordinal diagnostic rather than a pre-agreed failure threshold. Inspect both expected survivors and singleton coverage, not only the mean solve count, before changing thresholds. Narrative clues can trivially identify an answer and require a separate text-level audit.

## Provisional gates and implementation plan

These are **proposed**, not existing or approved promotion criteria: opening group minimum >=4; no probe identifies >25% of any opening group as non-winning singletons; entropy and minimax bot solve >=95% within six and >=99% within eight; report all target paths and compare a defined bounded-lookahead adversary. Run sensitivity analysis using plausible nonuniform and novice priors, inspect all clue wording and accessible scaffold variants, and evaluate separate novice and chemistry-experienced human cohorts. Preregister sample size, completion/abandonment and satisfaction measures, acceptable failure rates, and confidence intervals before playtesting; revise gates only with documented rationale. Human accessibility, comprehension, perceived fairness, and adequate sample coverage must pass along with quantitative gates. No evidence of such human testing is available here.

1. Complete the scientific-data migration and freeze an immutable rules/structures version; regenerate 78-entry clues from corrected data and rerun this benchmark and browser parity checks.
2. Add a separate practice mode and storage namespace keyed by roster and scientific rules version, seed and session ID; persist or reconstruct immutable feedback for existing guesses. Keep all practice completion and reset actions out of daily streak and classroom scoring paths. Provide eight guesses and the 4/5 clue progression described above.
3. Add deterministic Python/JavaScript feedback vectors and browser tests for all clue unlock boundaries, reload and mid-session updates, mode isolation, reproducible seed and no streak changes. Test accessibility and mobile behavior, then run a blinded playtest with predeclared measures and configuration assignment.
4. Record results and issue a separate promotion decision only after gates, data integrity, persistence, performance, and human evidence pass. Any scientific correction, clue/rule change, or roster change requires a fresh calibrated benchmark; release numbers alone do not define puzzle identity.

## Corrected-data rerun (v0.3.1 working tree)

Reproduce with `.\.venv\Scripts\python.exe tools/challenge_benchmark.py --output tools/challenge_benchmark_corrected_results.json`. The output identifies exact inputs by SHA-256 and includes every target trajectory. On the corrected catalog the expanded Ranked-like entropy, expected-survivors and minimax solvers each solve 10 targets in two guesses, 31 in three, 24 in four, 12 in five and one in six; bounded lookahead solves 6, 32, 26, 13 and 1 respectively. All four solve 78/78 within six and eight, and minimum opening k remains six.

The proposed 25% singleton gate **fails**: ascorbic acid and shikimic acid each distinguish 3 of 11 non-winning targets (27.3%) in one opening group. The tighter-feedback candidate's worst opening result is 2/11 (18.2%), though other groups and human interpretation still need review. This small comparison suggests the tighter candidate is worth further study, not that it is globally optimal or approved. Do not implement the current Ranked-like configuration as an accepted challenge design without resolving the probe gate. The solver-only prototype remains experimental; no interactive Challenge Practice mode or human playtest has been completed.
