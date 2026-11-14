# Chemle v0.3 Pilot Benchmark Report

- **Dataset Hash:** `488378e66071` (Ranked Schedule: `9257bd168fa6` unchanged)
- **Pools:** Ranked `36` compounds / Learn `78` compounds
- **Bundle Footprint:** `data/puzzles.json` `3253.93 KB`, `data/puzzles_data.js` `3027.99 KB`, `dist/biowordle_standalone.html` `3159.91 KB`

## 1. Full Solver Matrix (Strategies × Probe Modes × Profile Rules)

| Scenario | Strategy | Probe Mode | Clues | Targets | Probes | Cap | Solved | Mean | Median | P95 | Max | <=2% | <=6% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 36 Ranked (Ranked rules) | `greedy_entropy` | `full_guess_pool` | ON (ranked) | 36 | 36 | 6 | 100.0% (36/36) | 3.5 | 3.0 | 5 | 5 | 13.89% | 100.0% |
| 36 Ranked (Ranked rules) | `greedy_entropy` | `candidate_only` | ON (ranked) | 36 | 36 | 6 | 100.0% (36/36) | 3.69 | 4.0 | 5 | 5 | 25.0% | 100.0% |
| 36 Ranked (Ranked rules) | `min_expected_survivors` | `full_guess_pool` | ON (ranked) | 36 | 36 | 6 | 100.0% (36/36) | 3.53 | 3.0 | 5 | 5 | 11.11% | 100.0% |
| 36 Ranked (Ranked rules) | `min_expected_survivors` | `candidate_only` | ON (ranked) | 36 | 36 | 6 | 100.0% (36/36) | 3.69 | 4.0 | 5 | 5 | 25.0% | 100.0% |
| 36 Ranked (Ranked rules) | `greedy_entropy` | `full_guess_pool` | OFF | 36 | 36 | 6 | 86.11% (31/36) | 4.69 | 5.0 | 7 | 7 | 11.11% | 86.11% |
| 36 Ranked (Ranked rules) | `greedy_entropy` | `candidate_only` | OFF | 36 | 36 | 6 | 58.33% (21/36) | 5.03 | 5.5 | 7 | 7 | 13.89% | 58.33% |
| 36 Ranked (Ranked rules) | `min_expected_survivors` | `full_guess_pool` | OFF | 36 | 36 | 6 | 86.11% (31/36) | 4.75 | 5.0 | 7 | 7 | 8.33% | 86.11% |
| 36 Ranked (Ranked rules) | `min_expected_survivors` | `candidate_only` | OFF | 36 | 36 | 6 | 58.33% (21/36) | 5.03 | 5.5 | 7 | 7 | 13.89% | 58.33% |
| 36 Ranked (Learn rules) | `greedy_entropy` | `full_guess_pool` | ON (learn) | 36 | 36 | 8 | 100.0% (36/36) | 1.83 | 2.0 | 2 | 2 | 100.0% | 100.0% |
| 36 Ranked (Learn rules) | `greedy_entropy` | `candidate_only` | ON (learn) | 36 | 36 | 8 | 100.0% (36/36) | 1.83 | 2.0 | 2 | 2 | 100.0% | 100.0% |
| 36 Ranked (Learn rules) | `min_expected_survivors` | `full_guess_pool` | ON (learn) | 36 | 36 | 8 | 100.0% (36/36) | 1.83 | 2.0 | 2 | 2 | 100.0% | 100.0% |
| 36 Ranked (Learn rules) | `min_expected_survivors` | `candidate_only` | ON (learn) | 36 | 36 | 8 | 100.0% (36/36) | 1.83 | 2.0 | 2 | 2 | 100.0% | 100.0% |
| 36 Ranked (Learn rules) | `greedy_entropy` | `full_guess_pool` | OFF | 36 | 36 | 8 | 100.0% (36/36) | 1.97 | 2.0 | 2 | 2 | 100.0% | 100.0% |
| 36 Ranked (Learn rules) | `greedy_entropy` | `candidate_only` | OFF | 36 | 36 | 8 | 100.0% (36/36) | 1.97 | 2.0 | 2 | 2 | 100.0% | 100.0% |
| 36 Ranked (Learn rules) | `min_expected_survivors` | `full_guess_pool` | OFF | 36 | 36 | 8 | 100.0% (36/36) | 1.97 | 2.0 | 2 | 2 | 100.0% | 100.0% |
| 36 Ranked (Learn rules) | `min_expected_survivors` | `candidate_only` | OFF | 36 | 36 | 8 | 100.0% (36/36) | 1.97 | 2.0 | 2 | 2 | 100.0% | 100.0% |
| 78 Learn (Learn rules) | `greedy_entropy` | `full_guess_pool` | ON (learn) | 78 | 78 | 8 | 100.0% (78/78) | 1.95 | 2.0 | 2 | 2 | 100.0% | 100.0% |
| 78 Learn (Learn rules) | `greedy_entropy` | `candidate_only` | ON (learn) | 78 | 78 | 8 | 100.0% (78/78) | 1.94 | 2.0 | 2 | 3 | 98.72% | 100.0% |
| 78 Learn (Learn rules) | `min_expected_survivors` | `full_guess_pool` | ON (learn) | 78 | 78 | 8 | 100.0% (78/78) | 1.95 | 2.0 | 2 | 2 | 100.0% | 100.0% |
| 78 Learn (Learn rules) | `min_expected_survivors` | `candidate_only` | ON (learn) | 78 | 78 | 8 | 100.0% (78/78) | 1.94 | 2.0 | 2 | 3 | 98.72% | 100.0% |
| 78 Learn (Learn rules) | `greedy_entropy` | `full_guess_pool` | OFF | 78 | 78 | 8 | 100.0% (78/78) | 2.13 | 2.0 | 3 | 3 | 85.9% | 100.0% |
| 78 Learn (Learn rules) | `greedy_entropy` | `candidate_only` | OFF | 78 | 78 | 8 | 100.0% (78/78) | 2.13 | 2.0 | 3 | 3 | 85.9% | 100.0% |
| 78 Learn (Learn rules) | `min_expected_survivors` | `full_guess_pool` | OFF | 78 | 78 | 8 | 100.0% (78/78) | 2.13 | 2.0 | 3 | 3 | 85.9% | 100.0% |
| 78 Learn (Learn rules) | `min_expected_survivors` | `candidate_only` | OFF | 78 | 78 | 8 | 100.0% (78/78) | 2.13 | 2.0 | 3 | 3 | 85.9% | 100.0% |
| 78 Learn (Ranked rules stress, cap=6) | `greedy_entropy` | `full_guess_pool` | OFF | 78 | 78 | 6 | 66.67% (52/78) | 5.35 | 5.0 | 7 | 7 | 3.85% | 66.67% |
| 78 Learn (Ranked rules stress, cap=6) | `greedy_entropy` | `candidate_only` | OFF | 78 | 78 | 6 | 41.03% (32/78) | 5.73 | 7.0 | 7 | 7 | 7.69% | 41.03% |
| 78 Learn (Ranked rules stress, cap=6) | `min_expected_survivors` | `full_guess_pool` | OFF | 78 | 78 | 6 | 66.67% (52/78) | 5.35 | 5.0 | 7 | 7 | 3.85% | 66.67% |
| 78 Learn (Ranked rules stress, cap=6) | `min_expected_survivors` | `candidate_only` | OFF | 78 | 78 | 6 | 41.03% (32/78) | 5.73 | 7.0 | 7 | 7 | 7.69% | 41.03% |
| 78 Learn (Ranked rules stress, cap=8) | `greedy_entropy` | `full_guess_pool` | ON (learn) | 78 | 78 | 8 | 100.0% (78/78) | 2.49 | 2.5 | 3 | 3 | 50.0% | 100.0% |
| 78 Learn (Ranked rules stress, cap=8) | `greedy_entropy` | `candidate_only` | ON (learn) | 78 | 78 | 8 | 100.0% (78/78) | 2.54 | 3.0 | 3 | 3 | 38.46% | 100.0% |
| 78 Learn (Ranked rules stress, cap=8) | `min_expected_survivors` | `full_guess_pool` | ON (learn) | 78 | 78 | 8 | 100.0% (78/78) | 2.49 | 2.5 | 3 | 3 | 50.0% | 100.0% |
| 78 Learn (Ranked rules stress, cap=8) | `min_expected_survivors` | `candidate_only` | ON (learn) | 78 | 78 | 8 | 100.0% (78/78) | 2.54 | 3.0 | 3 | 3 | 38.46% | 100.0% |

## 2. Learn Pilot Family & Category Pack Evaluation (`greedy_entropy`, Cap=8)

### 2.1 Pilot Biochemical Families (`42` New Learn Compounds)

| Pilot Family | Count | Full-Pool Mean (Max) | Family-Only Mean (Max) | Family Solved <=8 |
|---|---|---|---|---|
| Aromatic degradation | 6 | 1.83 (2) | 1.83 (2) | 100.0% |
| Microbial and related natural products | 6 | 1.33 (2) | 1.33 (2) | 100.0% |
| Neuroactive metabolites | 6 | 1.33 (2) | 1.33 (2) | 100.0% |
| Phenylpropanoids and plant metabolites | 6 | 1.33 (2) | 1.33 (2) | 100.0% |
| Primary metabolism | 6 | 1.5 (2) | 1.5 (2) | 100.0% |
| Purines and nucleobases | 6 | 1.67 (2) | 1.67 (2) | 100.0% |
| Therapeutic analogues | 6 | 1.67 (2) | 1.67 (2) | 100.0% |

### 2.2 UI Category Packs (Learn Mode, All `78` Compounds)

| Category Pack | Count | Exploratory (`k_min<3`) | Full-Pool Mean (Max) | Pack-Only Mean (Max) | Pack Solved <=8 |
|---|---|---|---|---|---|
| Alkaloid / Metabolite | 5 | `true` | 1.4 (2) | 1.4 (2) | 100.0% |
| Biogenic Amine | 3 | `false` | 1.67 (2) | 1.67 (2) | 100.0% |
| Cofactor / Vitamin | 2 | `true` | 1.0 (1) | 1.0 (1) | 100.0% |
| Hormone / Metabolite | 1 | `true` | 1.0 (1) | 1.0 (1) | 100.0% |
| Hormone / Neurotransmitter | 2 | `true` | 1.0 (1) | 1.0 (1) | 100.0% |
| Hydrocarbon Degradation | 13 | `true` | 1.77 (2) | 1.77 (2) | 100.0% |
| Microbial Natural Product | 10 | `true` | 1.6 (2) | 1.6 (2) | 100.0% |
| Neurotransmitter | 5 | `true` | 1.4 (2) | 1.4 (2) | 100.0% |
| Plant Natural Product | 11 | `true` | 1.55 (2) | 1.55 (2) | 100.0% |
| Primary Metabolite | 14 | `true` | 1.71 (2) | 1.71 (2) | 100.0% |
| Therapeutic Drug | 12 | `true` | 1.58 (2) | 1.58 (2) | 100.0% |

## 3. Structural Density & Collision Comparison (`36` Ranked vs `78` Learn)

| Metric | Ranked (`36`) | Learn (`78`) |
|---|---|---|
| Nearest-neighbor Tanimoto mean / median | 30.75% / 28.0% | 45.64% / 49.0% |
| Nearest-neighbor Tanimoto p95 / max | 60% / 60% | 70% / 74% |
| Distinct Murcko cyclic scaffolds | 24 | 34 |
| Shared Murcko scaffold collision groups | 2 | 8 |
| Exact 5-descriptor (`MW/cLogP/TPSA/HBD/HBA`) collisions | 1 | 4 |
| Identical feedback across all catalog probes (Ranked / Learn) | 0 / 0 | 0 / 0 |

- **Exact 5-property collisions in Learn (`78`):** `[['4_hydroxybenzoic_acid', 'salicylic_acid'], ['gentisic_acid', 'protocatechuic_acid'], ['m_xylene', 'o_xylene', 'p_xylene'], ['theobromine', 'theophylline']]` (separated cleanly by Morgan ECFP4 Tanimoto %).

- **Closest Learn Pairs (Tanimoto %):** `curcumin` vs `ferulic_acid` (74%); `l_phenylalanine` vs `l_tyrosine` (70%); `ferulic_acid` vs `sinapic_acid` (68%); `5_htp` vs `tryptophan` (68%); `ciprofloxacin` vs `norfloxacin` (67%); `l_dopa` vs `l_tyrosine` (65%)

## 4. Close-Series Solver Walkthroughs (Learn Mode, `greedy_entropy`, Full Pool)

### Methylxanthines (Caffeine / Theophylline / Theobromine)

- **`caffeine` (`2` guesses):** G1: `ciprofloxacin` (13→1 cand, H=3.7b) -> G2: `caffeine` (1→1 cand, H=0.0b)
- **`theobromine` (`2` guesses):** G1: `ciprofloxacin` (13→1 cand, H=3.7b) -> G2: `theobromine` (1→1 cand, H=0.0b)
- **`theophylline` (`2` guesses):** G1: `ciprofloxacin` (13→1 cand, H=3.7b) -> G2: `theophylline` (1→1 cand, H=0.0b)

### Hydroxycinnamates (Caffeic / Ferulic / p-Coumaric acid)

- **`caffeic_acid` (`2` guesses):** G1: `ciprofloxacin` (13→1 cand, H=3.7b) -> G2: `caffeic_acid` (1→1 cand, H=0.0b)
- **`ferulic_acid` (`2` guesses):** G1: `ciprofloxacin` (13→1 cand, H=3.7b) -> G2: `ferulic_acid` (1→1 cand, H=0.0b)
- **`p_coumaric_acid` (`2` guesses):** G1: `phenol` (19→1 cand, H=4.248b) -> G2: `p_coumaric_acid` (1→1 cand, H=0.0b)

### Catecholamines (Dopamine / Norepinephrine / Adrenaline)

- **`adrenaline` (`2` guesses):** G1: `ciprofloxacin` (13→1 cand, H=3.7b) -> G2: `adrenaline` (1→1 cand, H=0.0b)
- **`dopamine` (`2` guesses):** G1: `ferulic_acid` (16→1 cand, H=4.0b) -> G2: `dopamine` (1→1 cand, H=0.0b)
- **`norepinephrine` (`2` guesses):** G1: `ferulic_acid` (16→1 cand, H=4.0b) -> G2: `norepinephrine` (1→1 cand, H=0.0b)

## 5. Cold Mobile Chromium (`390×844`, `4×` CPU) Benchmark

- **Configuration:** 390x844 mobile viewport, `4x CPU + Fast-3G (1.6 Mbps / 750 kbps / 150ms RTT) HTTP; 4x CPU file://`, `5` cold HTTP Fast-3G runs + `5` cold offline `file://` standalone runs (uncompressed `puzzles_data.js` = `3027.99 KB`, `biowordle_standalone.html` = `3159.91 KB`).
- **External Network Requests:** `0` across all runs (`all_zero_external_requests=True`).

| Metric | HTTP Fast-3G Median (ms) | HTTP Fast-3G Max (ms) | `file://` Median (ms) | `file://` Max (ms) |
|---|---|---|---|---|
| Cold Boot (`#chemle-search-input` + well ready) | 16264.2 | 16487.0 | 1013.5 | 1163.0 |
| Toggle to Learn (`78 Mols` active) | 151.0 | 233.1 | 151.1 | 165.8 |
| Search Filter (`Xylene`) | 163.7 | 275.1 | 129.1 | 162.4 |
| Submit Guess (`m_xylene` card rendered) | 151.6 | 168.2 | 126.8 | 477.0 |

## 6. v0.3 Ranked Promotion Gate Evaluation

- **Decision:** `promote_to_ranked = False` — **Keep Ranked Chemle at `N=36`**, while shipping the curated `78`-compound expansion in **Learn Mode**.
- **Rationale:** Ranked Chemle stays at N=36 (schedule hash 9257bd168fa6) to preserve historical daily schedules and calibrated difficulty (N=36 Ranked mean=3.50, max=5, 100% solved<=6), while the 78-compound pool ships in Learn Mode (mean=1.95, 100% solved<=2 under Learn rules) ahead of empirical blind human playtesting and Ranked clue/threshold retuning (78 unconditioned Ranked stress yields mean=5.35, max=7, solved<=6=66.67%).

