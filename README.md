# BioWordle v0.3.0

An interactive chemical (**Chemle**) and peptide motif (**Seqle**) scientific deduction suite built with deterministic RDKit/UniProt compilation, manifest-driven molecular provenance (`data/chemle_compounds.json`), zero external runtime dependencies (`dist/biowordle_standalone.html`), cross-language golden-vector verification, and portable Playwright Chromium E2E tests.

## Game Modes, Profiles & Exact Thresholds

### Profiles & Routes
- **Ranked Daily (`?seed=` optional for classroom seeds)**: UTC-day deterministic puzzle selection (`selectDailyIndex`) locked per game schedule hash (`36` Ranked Chemle / `28` Ranked Seqle).
- **Learn Mode (`?profile=learn`)**: Open practice profile with `8` guesses, manual `💡 Hint` progression, live rounded Tanimoto %, live shared functional groups, Seqle positional class matching (`≈`), and Chemle Category Pack filtering (`AnsPool` scoped to category over the expanded `78`-compound Learn pool while `GuessPool` retains all `78` Learn-eligible compounds as probes).

### Clue & Unlock Thresholds
| Game | Profile | Max Guesses | Level 0 (Start) | Level 1 Unlock | Level 2 Unlock | Level 3 / Scaffold Unlock |
| :--- | :--- | :---: | :--- | :--- | :--- | :--- |
| **Chemle** | Ranked | `6` | k-anonymous size band | Super-class (Guess `4`) | Role & mechanism (Guess `4`) | Biosynthesis (Guess `5`), Bemis–Murcko scaffold / carbon skeleton SVG (`scaffold_unlock_after_guesses: 5`) |
| **Chemle** | Learn | `8` | k-anonymous size band | Guess `1` or Hint | Guess `2` or Hint (+ Scaffold SVG) | Guess `3` or Hint |
| **Seqle** | Ranked | `6` | k-anonymous charge band | Architecture & terminal (Guess `2`) | Context requirement (Guess `4`) | Motif rule & role (Guess `5`) |
| **Seqle** | Learn | `8` | Charge × Hydropathy band | Guess `1` or Hint | Guess `2` or Hint | Guess `3` or Hint |

All pre-unlock structure wells (`chemle-preunlock`, `seqle-preunlock`) render strictly target-independent placeholder geometry before reveal/scaffold unlock.

## Curated Pools & Provenance (`v0.3.0`)
- **Chemle (`78` curated compounds: `36` Ranked / `78` Learn)**:
  - Single source of truth in `data/chemle_compounds.json`, validated offline on build via `validate_chemle_manifest()` in `chem_engine.py` (SMILES parsing, single-component check, parent uncharged InChIKey deduplication, formula/stereochemistry/racemate checks) and online via `tools/verify_chemle_manifest_online.py` against PubChem PUG REST & EMBL-EBI OLS4 ChEBI (`78/78` resolved).
  - **Ranked baseline (`36` compounds, `chemle_schedule_hash: 9257bd168fa6`)**: Unchanged from v0.2.2 so Ranked daily schedules, Level-0 quantiles, and solver calibration stay strictly isolated from Learn additions (`29` `verified` + `7` `baseline_legacy` flagged transparently post-reveal in the UI provenance panel).
  - **Learn pilot expansion (`42` new Learn-only compounds across `7` biochemical families)**: Microbial & related natural products, aromatic degradation, neuroactive metabolites, purines & nucleobases, therapeutic analogues, phenylpropanoids & plant metabolites, and primary metabolism (`CHEMLE_v0.3_MOLECULAR_AUDIT.md`).
- **Seqle (`30` motifs: `28` Ranked / `30` Learn)**:
  - Retired invalid COPI example `KKMPQ` and replaced it in Ranked with verified UniProt `P39656` / PDB `4J77` (`KEKSD`), updating `seqle_schedule_hash` to `623f02623c18`.
  - Learn-only pedagogical abstractions (`CCHHC`, `LDLAS`) remain restricted to Learn mode (`game_eligibility: ["learn"]`).
  - Every Seqle entry is governed 1:1 by `SEQLE_REFERENCE_REGISTRY` in `chem_engine.py` with typed citations (`verified_uniprot` with verified residue intervals vs `consensus_exemplar` such as `NPVYK`/`NLTGY`/`KENIP` vs `context_only` such as `PPPVP`/`1GBR` vs `abstraction`), cleanly separating **Functional Rule** from **Teaching Sequence** / **Parent Sequence Window**.

## Persistence, Schedule Isolation & Concurrency
- **V2 Stats Ledger (`biowordle_stats_v2`)**: Preserves v1 cumulative wins (`chemWins`, `seqWins`) via `migrateV1Stats`, stores deduplicated per-UTC daily win dates per game (`rankedDailyWins`), non-daily win IDs (`nonDailyWinIds`, incrementing win totals without changing daily streaks), and computes display streaks dynamically via `getDisplayStreak` (dropping to `0` if a UTC day is missed).
- **Cross-Tab Serialization**: Uses `navigator.locks` on HTTP origins (falling back synchronously on `file://` standalone) plus a durable `biowordle_pending_wins_v2` write intent queue drained on startup, paired with a `storage` event listener and stale-move conflict guard (`Concurrent Tab Sync`).
- **Per-Game Legacy Key Compatibility**: `chemle_schedule_hash` (`9257bd168fa6`) is decoupled from `seqle_schedule_hash` (`623f02623c18`). Recognized legacy schedule keys (`["9257bd168fa6", "8bb96ac156cc", "919ad8afb190"]`) are restored cleanly only when the persisted `targetId` matches the active daily target; incompatible records (e.g., retired `KKMPQ` or a rotated target) are archived safely (`biowordle_legacy_archive_v2`) instead of re-scoring old guesses against a new answer.

## Reproducibility, Build & Test Commands

### Line-Ending Normalization
`build_puzzles.py` writes all generated tracked artifacts (`data/puzzles.json`, `data/puzzles_data.js`, `tests/golden_feedback.json`, and `dist/biowordle_standalone.html`) with strict `LF` (`\n`) newlines so cross-platform CI diffs remain clean on Windows and Linux.

### Python Compiler & Verification
```bash
# Regenerate puzzles.json, puzzles_data.js, golden_feedback.json, and dist/biowordle_standalone.html (LF-normalized)
python build_puzzles.py

# Run solver-bot calibration gates (Chemle & Seqle Ranked/Learn + deterministic 10k legal-probe run)
python tools/solver_bot.py

# Run Python unit & audit suite (provenance fixtures, solver clue/unsolved/probe tests, schedule hash isolation)
python -m unittest discover -v -s . -p "test_*.py"

# Verify live UniProtKB (2026_04) & RCSB PDB coordinates online
python tools/verify_provenance_online.py

# Verify live PubChem PUG REST & EMBL-EBI OLS4 ChEBI structures for all 78 manifest compounds
python tools/verify_chemle_manifest_online.py

# Run v0.3 Chemle expanded-pool solver matrix & 5-run cold mobile Chromium benchmark
python tools/chemle_pilot.py
```

### Node Parity & Portable Playwright Chromium E2E
```bash
# Cross-language golden vector & stats-ledger unit suite
node tests/verify_golden_and_dom.mjs

# Install Playwright test dependencies and Chromium (with Linux deps in CI via pw:install:ci)
npm ci
npm run pw:install

# Run full Playwright E2E suite against loopback HTTP server (tests/e2e/static_server.mjs) + file:// offline smoke
npm run test:e2e
```
On Windows sandboxed environments where loopback binding/Chromium launch is restricted inside the sandbox, run `npm.cmd run test:e2e` with `BypassSandbox: true`.

## Known Limitation & Phase 3 Seqle Prototype Findings
- **Known Balance Ceiling (`N=28` Ranked Seqle)**: Under `tools/solver_bot.py`, Ranked Seqle solves catalog-only at `mean=2.00` (`median=2.0`, `92.86% <=2`) and hits an algorithmic legal-probe upper bound of `mean=1.93` (`median=2.0`, `100.0% <=2`) across `10,000` legal 5-AA probes (`seed=20261114`), because five positional tiles across `20` amino acids yield enough entropy to split 28 targets in 2 guesses.
- **60 / 90 / 120 Prototype Evaluation (`tools/seqle_prototype.py`)**: Tested verified UniProt 5-mer expansion (`data/seqle_prototype_corpus.json`, `total_count=120`) across `opening_charge` vs delayed `charge_after_guess_2` schedules and `greedy_entropy` vs exact `two_step` solvers:
  - Even at `N=120` with `charge_after_guess_2`, catalog-only reaches `mean=2.54, median=3.0, <=2=45.0%`, while unrestricted legal probes (`seed=20261114` 10k sample & `seed=20261115` 5k holdout) achieve `mean=2.46–2.50, median=2.0–2.5, <=2=50.0%–54.2%, <=6=100.0%`, falling short of the promotion thresholds (`legal_probe median >= 3.0, mean >= 2.6, <=2 <= 45.0%`).
  - Therefore `any_config_promoted=False`: v0.2.2+ retains `N=28` Ranked Seqle pending blind human playtesting on an expanded pool (`data/seqle_prototype_report.json`).

## Chemle v0.3 Expanded-Pool Pilot Findings (`tools/chemle_pilot.py`)
- **Expanded Learn Pool (`36` Ranked → `78` Learn)**: All `42` curated additions across `7` biochemical families are separated from Ranked daily schedules (`chemle_schedule_hash: 9257bd168fa6`, `CHEM101-WEEK3 -> pyocyanin`, `2026-11-14 -> artemisinin`) via `game_eligibility` / `guess_eligibility` in `data/chemle_compounds.json`.
- **Solver & Density Evaluation (`data/chemle_pilot_report.json`, `CHEMLE_v0.3_PILOT_BENCHMARK.md`)**:
  - Nearest-neighbor Morgan ECFP4 similarity rises from `mean=30.75%` (`max=60%`) in `36` Ranked to `mean=45.64%` (`max=74%`) in `78` Learn, with `34` distinct Murcko scaffolds, `4` exact 5-descriptor collisions (all cleanly separated by Tanimoto % and feedback), and `0` indistinguishable pairs across all catalog probes.
  - With Learn feedback + clues ON (`cap=8`), both `greedy_entropy` and `min_expected_survivors` solve `78/78` (`100.0%`, `mean=1.94–1.95`, `max=2–3`). Under unconditioned Ranked feedback (`cap=6`) stress check, `78` compounds drop to `66.67%` solved (`mean=5.35`, `max=7`), confirming `promote_to_ranked=False` until empirical human playtesting and Ranked clue retuning.
- **Offline Mobile Chromium (`390×844`, `4×` CPU + Fast-3G CDP)**: Zero external network requests across all `5` cold HTTP Fast-3G runs and `5` cold `file://` standalone runs; `file://` boot median is `1013.5 ms` and interactive Learn switch / filter / guess submission medians stay `126.8–151.1 ms` (`16264.2 ms` cold HTTP boot over raw un-gzipped `1.6 Mbps` Fast-3G due to `3027.99 KB` `data/puzzles_data.js`).
