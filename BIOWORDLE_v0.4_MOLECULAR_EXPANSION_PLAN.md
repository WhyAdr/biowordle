# BioWordle v0.4 — Molecular expansion plan

**Base:** `03febdf790c39b6d14202833af671af469146ba2` · **Status:** research proposal; no production implementation authorized by this handoff · **Date:** 2026-11-14.

## Decision and scope

The v0.3 proposal's seven families and 42 names were already incorporated into the current 78. Following the user's decision, this plan proposes **42 different candidates**, six in each of the same seven teaching families, for a 120-entry research build. The [candidate intake](BIOWORDLE_v0.4_CANDIDATE_MANIFEST.json) has all 42 marked `pending` and with empty production answer/guess eligibility. It passed the current offline manifest validator **in combination with all 78 existing records**. That is validation under the existing code, not a declaration that 42 candidates are scientifically accepted. Generic placeholder hints must be replaced only after claim-level review.

The present recommendation is **do not implement or ship the 120 pool yet**. The user has reserved implementation for after scientific review of this plan and approval of the benchmark strategy. The expansion must not alter the 36-answer/36-probe Ranked Daily roster or silently reset its sessions. The 150-compound option below is a separate checkpoint, not approval to add those molecules.

## Candidate intake and evidence

| Family | Six new proposals |
| --- | --- |
| Aromatic degradation | Toluene; benzene; ethylbenzene; biphenyl; phenylacetic acid; phthalic acid |
| Phenylpropanoids and plant metabolites | Gallic acid; syringic acid; thymol; quercetin; kaempferol; apigenin |
| Therapeutic analogues | Isoniazid; metronidazole; sulfamethoxazole; trimethoprim; lidocaine; dapsone |
| Neuroactive metabolites | 2-Phenylethylamine; agmatine; putrescine; cadaverine; spermidine; taurine |
| Primary metabolism | L-Alanine; L-leucine; L-serine; L-glutamic acid; L-lactic acid; 2-oxoglutaric acid |
| Microbial and related natural products | Indole; 3-methylindole; kojic acid; indigo; itaconic acid; tropolone |
| Purines and nucleobases | Guanine; cytosine; thymine; uric acid; uridine; inosine |

Every row records a PubChem CID, ChEBI accession, selected input graph, molecular formula, standardized SMILES, InChIKey, stereo/protonation policy, source links and unresolved review status. The named biological categories are organizational labels, not a claim that every polyamine is a neurotransmitter or every compound is naturally produced in every organism. Seven study leads for the latter four families and 18 for the former three are attached to individual rows; **17 rows do not have a claim-bound biological study lead in the manifest itself**; their separately located leads are in the [literature ledger](BIOWORDLE_v0.4_LITERATURE_LEADS.md). All 42 leads need full-study review and precise mapping of each final clue to its support. The independent first batch's [methods and context](tools/v04_candidates_a_notes.md), and [saved source records](tools/v04_candidate_b_source_records.json) for the second batch, show what was queried. The fetch/preparation scripts are research tooling; network sources can change.

All 42 selected PubChem keys match the cited ChEBI standard InChIKeys, but matching keys do not settle tautomer or bond placement. For nine second-batch candidates the displayed ChEBI and PubChem isomeric graphs either differ or the ChEBI SMILES cannot be parsed by RDKit: agmatine, indole, 3-methylindole, guanine, cytosine, thymine, uric acid, uridine and inosine. These remain open, not silently treated as identical graphs. The selected neutral species is a reproducible game representation, not an assertion that it dominates at physiological pH. Indigo uses the specific E-configured PubChem CID 5318432 and ChEBI CHEBI:749117; a broad-name lookup returned a different CID and tautomer, demonstrating why names alone are insufficient. Uric acid uses a specific ChEBI tautomer alongside the broad parent concept. Verify the source graphs and intended tautomer with an expert before accepting either.

The additional [duplicate-screen report](tools/v04_manifest_audit_results.json) found no candidate-involving repeated InChI connectivity block or repeated RDKit canonical tautomer (ignoring stereo); these are screening results, not proof of distinct identities. Tautomer enumeration hit its transform limit for existing violacein and proposed quercetin, so those two screens are incomplete. One candidate-involving pair exceeded the Morgan 0.8 similarity threshold: putrescine–cadaverine (0.875); review both identity and feedback discrimination. Reproduce with `.\.venv\Scripts\python.exe tools/v04_manifest_audit.py`.

Metformin was considered and deferred: RDKit flags an unspecified potential amidine stereochemical element, and the current validator cannot express an unspecified non-racemate without an inaccurate label. Dapsone replaces it in the provisional 42; do not weaken the validator. Each selected source form still needs curator acceptance and an explicit salt/protonation and stereo decision. A single InChIKey, a database synonym or an accession resolving is not enough to establish an identity or an experimentally demonstrated biosynthesis.

The separate [literature ledger](BIOWORDLE_v0.4_LITERATURE_LEADS.md) identifies a primary-study lead for each of the 17 remaining candidates, with its limited experimental context. All leads are at the abstract or source-discovery stage; they are not approved biological claims and must be checked against the full studies. The 78 baseline records also require a claim-by-claim audit before their existing hints can be described as fully supported.

The base manifest has 71 `verified` and seven `baseline_legacy` entries, but those labels conflate structural and biological review. Prodigiosin and violacein retain review notes describing pre-v0.3.2 representations despite corrected current structures. Quinine's `unspecified_legacy` reflects RDKit's bridgehead-nitrogen stereo flag and the validator's coupling of status to legacy metadata. Resolve these contradictions through an explicit schema migration and evidence ledger, preserving the selected v0.3.2 structures unless new source evidence warrants a separately versioned correction. The historical v0.3 audit is a snapshot, not current confirmation. The current online verifier tests PubChem key equality for `verified` records and ChEBI accession existence; it does **not** establish matching ChEBI graphs, all stereochemistry, correct tautomers, or biological claims.

## Reproducible 78-versus-120 solver prototype

The prototype compiles a temporary copy through `build_all` and the normal structural validator. `--research-preview` **temporarily changes pending rows to playable and `verified` in memory only** to exercise the compiler; this is a test adapter, not scientific verification. It captures existing compiler gate failures rather than changing gate thresholds. No production manifest or generated bundle is written. Reproduce at the repository root:

```powershell
$env:PYTHONUTF8 = '1'
.\.venv\Scripts\python.exe tools/v04_solver_benchmark.py --manifest BIOWORDLE_v0.4_CANDIDATE_MANIFEST.json --research-preview --compiled-output tools/v04_compiled_preview.json --output tools/v04_solver_results.json
.\.venv\Scripts\python.exe tools/v04_solver_benchmark.py --manifest BIOWORDLE_v0.4_CANDIDATE_MANIFEST.json --research-preview --challenge-super-after 3 --output tools/v04_solver_early_clue_results.json
.\.venv\Scripts\python.exe tools/v04_space_audit.py --preview tools/v04_compiled_preview.json
```

The solver output has each target's path, misses, per-opening-group probe statistics, input hashes and compiler results. All targets and probes are evaluated under a uniform prior. Four policies are compared: greedy entropy, expected survivors, minimax, and bounded two-ply lookahead over the eight best first probes. The tie-break differs from the production compiler bot. The model conditions on automatic clue text and uses scaffold SMILES as a proxy for the displayed scaffold. Learn assumes no optional manual hint clicks; Challenge models the shipped clue timings and no narrative. Placeholder candidate hints are not final educational clues. Future clues, visible SVG equivalence, accessibility, nonuniform player familiarity and human strategy may change these results.

| Pool / profile | Entropy: distribution of solves by guess | Bounded lookahead: distribution | Worst opening non-winning singleton coverage | Min opening k |
| --- | --- | --- | ---: | ---: |
| 78 Learn | 1:3, 2:75 | 1:6, 2:71, 3:1 | 19/19, 100% | 6 |
| 120 Learn | 1:3, 2:116, 3:1 | 1:6, 2:112, 3:2 | 23/23, 100% | 13 |
| 78 Challenge | 2:15, 3:30, 4:19, 5:13, 6:1 | 2:14, 3:25, 4:25, 5:12, 6:2 | 2/11, 18.2% | 6 |
| 120 Challenge | 2:8, 3:44, 4:43, 5:24, 6:1 | 1:1, 2:7, 3:41, 4:45, 5:25, 6:1 | 2/17, 11.8% | 13 |

For the 120 pool, expected-survivor Challenge gives `2:9, 3:43, 4:43, 5:24, 6:1`; minimax gives `2:8, 3:45, 4:42, 5:23, 6:2`. All four policies report **120/120 within both six and eight**; no target is unsolved in this model. The captured production compiler gates report no failures for this provisional preview. These are model results, not evidence of human solve rate. Learn's 100% singleton examples reflect its intentionally informative feedback; do not apply the Challenge probe threshold to Learn without a separate learning-design decision.

The existing Challenge profile (`MW` 140, `cLogP` 3.5, similarity edges 60/80, superclass after four, scaffold after five, eight guesses) satisfies the existing simulated gates on this pool; keep it as the initial candidate rather than changing thresholds without evidence. As a sensitivity run, revealing superclass after three shifts the entropy distribution to `2:8, 3:44, 4:54, 5:13, 6:1`, reducing mean 3.717 to 3.625. It does not establish superior human difficulty and is not a proposed production change. Review the 11.8% worst opening probe against the actual visual clues and all candidate-specific outcomes; retain the eight-attempt budget, report any misses after future changes, and run human testing and nonuniform-prior sensitivity before a release decision.

## Ranked continuity: blocking engineering finding

The isolated 120 build has the same 36 Ranked answer and legal-probe IDs, unchanged Ranked clue ladders/descriptors and all 1,296 Ranked feedback signatures, and unchanged schedule hash `9257bd168fa6`. Nevertheless the current `chemle_rules_hash` covers all 78/120 compounds, full similarity matrix and all profile feedback rules. The prototype changes that hash from `7edd1b8864c9` to `2c81484214e1`. Ranked daily and classroom session keys include the global hash, so this would archive/restart in-progress or completed Ranked sessions despite a stable schedule. **This fails the requested historical-session acceptance criterion.**

Before implementation, design and test a profile-scoped Ranked rules identity containing all Ranked answers, legal probes, descriptors, visible clue content, relevant feedback rules, similarity submatrix and algorithm revision; keep Challenge identities independently sensitive to its expanded roster and rules. Do not simply reuse an old hash or exempt arbitrary differences. For the specific audited v0.3.2 global hash, allow an idempotent, atomic migration only after a frozen equivalence fixture proves Ranked roster, rules, clues and all pair signatures match. Preserve original payloads and ledger entries, including classroom seeds and completed results. Prove same-version reload, completed and in-progress migration, mismatch refusal/archive, stale old tabs writing during/after migration, simultaneous tabs, duplicate win credits, old-version reopening and rollback. Define conflict behavior explicitly; if equivalence is not established, preserve/archive without re-scoring and show an accurate notice. Challenge sessions need their own version policy and may not be reinterpreted under new catalog/rules. Keep old bundles available when historical feedback reconstruction is required.

## Distribution and engineering measurements

The [performance report](tools/v04_perf_notes.md) and [raw results](tools/v04_perf_results.json) compare generated temporary 78 and 120 bundles, not a deployment. Regenerate with:

```powershell
.\.venv\Scripts\python.exe tools/v04_perf.py --preview tools/v04_compiled_preview.json --output tools/v04_perf_results.json --browser --runs 2
```

| Measure | 78 baseline | 120 research preview |
| --- | ---: | ---: |
| `puzzles_data.js`, raw / gzip-9 | 3,133,800 / 342,731 B | 3,907,536 / 435,073 B |
| Standalone HTML, raw / gzip-9 | 3,284,629 / 376,466 B | 4,058,365 / 468,771 B |
| Chemle cross-language vectors | 18,252 | 43,200 |
| Golden file / Node verifier | 1,037,622 B / 0.631 s | 2,317,686 B / 0.721 s |
| Mobile-emulated gzip HTTP ready / transfer | 2,511 ms / 378,724 B | 2,941 ms / 471,066 B |
| Mobile-emulated ready JS heap | 4.47 MB | 5.09 MB |
| Local offline-file ready | 1,154 ms | 1,396 ms |

A single local base `build_all` observation was 3.68 s and the separately measured research preview 5.22 s; differing gate-capture paths and host noise make this exploratory, not a controlled speedup ratio. Browser medians use only two samples, a 390×844 Chromium viewport, 4× CPU slowdown and simulated Fast 3G over loopback. File navigation is unthrottled; JS heap is not total browser memory. Search and submit timings and raw samples are in the JSON. No page errors or external requests were observed. This is not a physical-device or hosted-HTTP test. The temporary 120 bundle passed browser smoke tests and the full Node parity verifier, **not** the complete Playwright E2E suite. The 120 preview is not release-certified.

The preview retains a self-contained standalone generator. Production implementation must run the complete portable Chromium HTTP and `file://` E2E suite, and must replace hard-coded Challenge `78` labels and 78-specific E2E expectations with current-count assertions while retaining explicit frozen v0.3.2 regression fixtures. Re-measure an actual deployed service with HTTP compression, cache/version headers, real mobile devices, interaction latency and process memory. Do not claim the byte projections for deferred SVG or game-asset splitting are functional implementations.

## Separate selection recommendation for 150

The 120 preview has 6 sulfur-containing and **zero phosphorus-containing** molecules, 18 empty Murcko scaffolds, 45 distinct Murcko scaffolds including empty, and nearest-neighbor Morgan similarity mean 0.472 (78 base: 3 sulfur, zero phosphorus, six empty scaffolds, 32 distinct and mean 0.457). These are representation-dependent coverage measurements, not a scientific quality score. The preview remains dominated by small compounds (median MW 168.13 Da) and the Learn model is already highly identifying after one probe. A 150 shortlist should fill gaps without relying solely on fingerprint distance or increasing near-duplicates.

Use the following **30 unverified nominations**, six per area, as a search queue, not a commitment to include them: (1) sulfur/thiols: L-cysteine, L-methionine, L-cystine, glutathione, ergothioneine, dimethyl sulfoxide; (2) phosphorylated metabolites: D-glucose 6-phosphate, D-glucose 1-phosphate, phosphoenolpyruvic acid, adenosine monophosphate, glycerol 3-phosphate, riboflavin 5-phosphate; (3) carbohydrates: D-glucose, D-fructose, D-ribose, D-galactose, sucrose, N-acetyl-D-glucosamine; (4) cofactors/vitamins: nicotinamide, riboflavin, thiamine, pyridoxal, pyridoxine, folic acid; (5) aliphatic metabolites: glycerol, D-mannitol, 1,3-propanediol, 4-hydroxybutyric acid, butyric acid, 3-hydroxybutyric acid. Names do not select charge, stereoisomer, tautomer or exact structure; several are challenging ions, mixtures or close isomers. Assess chemistry, biological relevance, stable atom mapping, safety of representations, source and licensing, and uniquely supported instructional claims before retaining them. Replace any failing candidate rather than filling a quota.

The distinct 150 checkpoint must repeat full 120-to-150 identity/source/parent/tautomer duplicate checks, source-mapped biological evidence, k-anonymity and all four clue-conditioned policies, per-target misses, adversarial probes, human playtest, compiler/runtime/size/memory benchmarks, browser/offline E2E and Ranked persistence equivalence. At 150, exhaustive three-profile Chemle parity would have 67,500 ordered vectors; measure it rather than extrapolating performance. No 150-compound promotion or implementation follows from this shortlist.

## Phases and acceptance tests for an executing agent

1. **Evidence freeze.** Resolve the nine source-graph cases, all 42 biological evidence ledgers and baseline contradictions. Bind every proposed mechanism/origin/hint to a specific study and organism/context, identify experimentally demonstrated versus predicted/generalized statements, and document salt, charge, stereochemistry and tautomer. Validate IDs, full keys, standardized uncharged parents, candidate/base duplicates, tautomer and close-isomer groups with human adjudication. An unresolved record stays non-playable; do not fabricate placeholders into verified claims.
2. **Compatibility before expansion.** Specify and test the profile-scoped rules-hash migration above using frozen v0.3.2 data and real legacy daily/classroom payloads. Require unchanged 36-answer schedule, all Ranked clues and full pair signatures, correct histories and streaks, cross-tab and rollback behavior. Seqle must remain isolated.
3. **Approved 120 implementation.** Only after scientific review and benchmark-strategy approval, add accepted records in isolated commits through the existing manifest/compiler. Regenerate pinned outputs and golden vectors, retain the standalone, use dynamic accessible counts, and run offline unit tests, full Python–JS parity, deterministic rebuild, full Chromium HTTP/file E2E, solver gates and CI on the actual commit. Keep Ranked unchanged; report explicitly if fewer than 42 candidates pass.
4. **Evaluation and decision.** Run predeclared 78/120 comparisons, all four policies, every target and probe, clue/scaffold visual audits, nonuniform priors and human playtests. Report unsolved targets within eight and all failed gates; do not weaken gates to force acceptance. Measure actual serving and physical mobile performance. Then make a separate review decision for the 120 Learn/Challenge experience. Ranked promotion is out of scope.
5. **150 checkpoint.** Only after the 120 decision, independently review and prioritize the 30 nominations and perform the full scientific, gameplay, persistence and performance evaluation above. Stop before implementing 150 until that checkpoint is explicitly accepted.

## Handoff validation

On the unchanged base, `PYTHONUTF8=1` with `python -B -m unittest discover -v -s . -p "test_*.py"` passed all 21 tests. The repository's portable Playwright command passed all 12 HTTP and standalone Chromium tests. The temporary 120 preview independently passed the generated full Node golden/DOM verifier and browser smoke tests described above; it did not receive the full Playwright suite. Candidate, base and preview SHA-256 identifiers were checked against the benchmark outputs. No production asset was changed and no implementation commit was created.

## Known limitations and reopen conditions

The candidate intake is structurally checked but not scientifically approved; nine cross-source graph cases, two incomplete tautomer enumerations and 17 manifest-to-ledger literature mappings remain, all 42 leads require full-study review, and existing biology claims need source mapping. The performance comparison is local and small-sample; the research preview does not have full E2E or deployed CI results. The solver assumes exact catalog knowledge and proxy clues. Any new chemical evidence, failed acceptance check, changed roster/rules/clues, or user request reopens the relevant phase; otherwise this document is a handoff for review, not permission to ship.
