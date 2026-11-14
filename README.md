# BioWordle v0.2 ("Restraint") — Chemle & Seqle (RDKit Edition)

**BioWordle** is a two-layer interactive cheminformatics and structural biology deduction suite comprising:
1. **🧪 Chemle**: Deduce a target metabolite, hydrocarbon degradation intermediate, natural product, or therapeutic drug (36 curated molecules) from **Chiral Morgan (`ECFP4`-like) Tanimoto similarity**, physicochemical directional deltas, functional group overlap, and progressive **Bemis–Murcko / Carbon-Skeleton 2D scaffolds**.
2. **🧬 Seqle**: Decode a 5-residue functional protein motif (`EKDEL`, `DEVDG`, `GRGDS`, `LPETG`, `GHSMG`, `HEIGH`, `KCVIM`, `DEADE`, `GVGKS`, `VDFGL`, `HRDLA`, `WSEWS`, `PPPVP`, `INFEK`, `YQRLN`, `KKMPQ`, `CAACH`, `NLTGY`, `DSGLD`, `KENIP`, `KKKRK`, `FGGFT`, `YGGFM`, `SYSGY`, `NPLYK`, `GHHPH`, plus Learn-only pedagogical abstractions `CCHHC` and `LDLAS`) using two-pass multiplicity-aware residue feedback and **Henderson–Hasselbalch** net charge at pH 7.4, Kyte–Doolittle GRAVY, and isoelectric point ($\text{pI}$).

---

## What's New in v0.2.1 (Correctness & Provenance Patch)

Built cleanly on top of v0.2.0 (`a30fb92`) without expanding the 36-molecule / 30-motif catalog or redesigning the UI:

### Gameplay, Runtime & Persistence Fixes
1. **Ranked Manual Hint Lockout**: In `ranked` mode, `#btn-hint` is disabled (`🔒 Auto-Clues (Ranked)`), `getUnlockedClueLevel()` ignores `hintLevelUnlocked` when not revealed, repeated pre-guess clicks cannot unlock clues early, and clue cards unlock strictly at `after_guesses` thresholds (`Level 1/2` after Guess 4, `Level 3` after Guess 5, scaffold preview after Guess 5).
2. **Learn-Mode Candidate Filtering Completeness**: `countCompatibleChemleCandidates()` and `countCompatibleSeqleCandidates()` incorporate **all player-visible feedback** (`SIM_PCT` exact Tanimoto, live shared `functional_groups`, `class` tiles in Learn Seqle, plus active Level-0 `k`-anonymous bands and unlocked clue/scaffold constraints).
3. **UTC Rollover Guard & Pre-Submit Stale Protection**: Each Ranked session tracks its `utcDate` and calls `checkAndHandleUtcRollover()` on the 1-second timer AND before accepting any guess, reveal, or mode/profile switch (testable via `window.__BIOWORDLE_TEST_UTC_DATE__`), preventing yesterday's guesses/wins from being credited to the next UTC day's streak.
4. **v1-to-v2 Stats Migration & Split Versioning**: `migrateV1Stats()` / `loadStatsV2()` preserves cumulative `chemWins` and `seqWins` from `biowordle_stats_v1` while resetting `streak` and `maxStreak` to `0` (since v1 lacked UTC win dates). Daily puzzle indexing and storage keys now hash `schedule_version_hash` (`9257bd168fa6`, derived solely from canonical catalog IDs and `game_eligibility`, with legacy `dataset_version_hash` read fallback) rather than `dataset_version_hash` (`8bb96ac156cc`), so future descriptor/feedback tuning does not scramble daily puzzles.

### Biological Reference Sequence & Provenance Audit (`chem_engine.py`, `build_puzzles.py`)
Every Seqle motif was audited against authoritative **UniProtKB/Swiss-Prot FASTA** sequences and **RCSB PDB** polymer entities and validated via `SEQLE_REFERENCE_REGISTRY` (`enforce_reference_registry=True`):
- **Human ITGB3 `NPLYR` -> `NPLYK` (`UniProt:P05106`)**: Removed the artificial `SEQLE_RAW_MOTIFS[28]["parent_sequence"]` override (`AKWDTANNPLYREATST`) and corrected the canonical ITGB3 cytoplasmic talin-binding turn (`770–774` in precursor / `744–748` mature, `PDB:1MK7`) inside `AKWDTANNPLYKEATSTFTNITYRGT` (`763–788`, offset `7`).
- **Corrected Misattributed / Inaccurate Native Exemplars (`7`)**: `GHSMG` -> `UniProt:G9BY57` (`163–167`, `PDB:4EB0`, replaces `A0A0K8P6T7` `GWSMG`); `HEIGH` -> `UniProt:P08254` (`MMP3_HUMAN`, `218–222`, `PDB:1B3D`); `VDFGL` -> `UniProt:P41279` (`MAP3K8`, `269–273`, `PDB:4Y83`); `HRDLA` -> `UniProt:P10721` (`KIT`, `790–794`, `PDB:1PKG`); `WSEWS` -> `UniProt:P32927` (`CSF2RB`, `425–429`, `PDB:1EGJ`, avoiding engineered `1EER` `A235E`); `CAACH` -> `UniProt:P07143` (`CY1_YEAST`, `101–105`, `PDB:3CX5`); `FGGFT` -> `UniProt:Q13519` cleaved mature Nociceptin (`130–146` `FGGFTGARKSARKLANQ`, `130–134`).
- **Canonical Precursor Coordinate Corrections (`4`)**: `GRGDS` (`P02751` `1614–1618`), `YIGSR` (`P07942` `950–954`), `SYSGY` (`P35637` `37–41`), `GHHPH` (`P04196` `383–387`).
- **Documented Consensus Exemplars (`4`) & Learn Abstractions (`2`)**: `KENIP`, `KKMPQ`, `NLTGY`, and `NPVYK` have `status: "consensus_exemplar"` (`parent_sequence == seq`, `window_offset == 0`) because their cited family proteins carry degenerate ELM/PROSITE variants rather than the literal 5-mer (we refuse to fabricate flanking sequences); `CCHHC` and `LDLAS` remain Learn-only `"abstraction"` entries.

### Clue-Conditioned Solver Calibration & Documented Bounds (`data/feedback_spec.json`, `tools/solver_bot.py`)
`run_chemle_greedy_solver` and `run_seqle_greedy_solver` now model Level-0 `k`-anonymous bands and progressive guess-threshold clue unlocks:

| Mode & Profile | Pre-Patch (`a30fb92`) | v0.2.1 (Clue-Conditioned) | Notes & Documented Limitations |
| :--- | :--- | :--- | :--- |
| **Chemle Ranked** | Cond: `mean=3.11`, `<=2=22.22%`, `dist={'1':1,'2':7,'3':18,'4':7,'5':3}` (Uncond: `3.89`, `dist={'1':1,'3':14,'4':11,'5':7,'6':3}`) | Cond: `mean=3.56`, `<=2=13.89%`, `max=6`, `er=13.61`, `dist={'2':5,'3':15,'4':9,'5':5,'6':2}` | Super-class moved to Guess 4, scaffold to Guess 5; widened Ranked `MW`/`cLogP`/Tanimoto bands (`[72, 88]`) |
| **Chemle Learn** | Cond: `mean=1.83`, `<=2=100%`, `dist={'1':6,'2':30}` (Uncond: `1.97`, `dist={'1':1,'2':35}`) | Cond: `mean=1.83`, `<=2=100%`, `max=2`, `er=1.0`, `dist={'1':6,'2':30}` | **Exploratory policy (Option A)**: Learn intentionally shows exact Tanimoto %, live shared FGs, & `TPSA`; gate relaxed rather than stripping educational telemetry |
| **Seqle Ranked** (`N=28`) | Uncond (`NPLYR`): `mean=1.96`, `max=2`, `YIGSR` splits `28/28`; Cond: `1.89`, `dist={'1':3,'2':25}` | Cond: `mean=2.00`, `<=2=92.86%`, `max=3`, `YIGSR` splits `19/28` (`max_probe=21`), `er=1.64`, `dist={'1':2,'2':24,'3':2}` (Uncond: `2.21`, `dist={'1':1,'2':20,'3':7}`) | Disabled positional `class` tiles in Ranked & widened non-directional strip tolerances. **Remaining $N=28$ bound**: 5-letter Wordle tiles + `charge_band` alone solve `28` targets in `mean=2.00` without strip cues (`N` catalog expansion deferred beyond v0.2.1) |
| **Seqle Learn** (`N=30`) | Cond: `mean=1.80`, `<=2=100%`, `dist={'1':6,'2':24}` | Cond: `mean=1.80`, `<=2=100%`, `max=2`, `YIGSR=28/30`, `dist={'1':6,'2':24}` (Uncond: `1.97`, `dist={'1':1,'2':29}`) | Full exploratory feedback (`class` tiles + tight `Q`/`GRAVY`/`pI`) |

---

## What's New in v0.2 ("Restraint")

- **Two Explicit Gameplay Profiles (`🏆 Ranked Daily` vs `🎓 Learn Mode`)** driven by a single authoritative shared specification (`data/feedback_spec.json`):
  - **Ranked Daily (`ranked`)**: 6 guesses, $k$-anonymous opening clues ($k_{\min} \ge 4$), 3-tier Tanimoto bands (`LO <61%`, `MD 61–84%`, `HI ≥85%`), calibrated descriptor tolerances (`MW`, `cLogP` in Chemle; `net_charge_ph74`, `gravy` in Seqle), and scaffold reveal gated until Guess 3.
  - **Learn Mode (`learn`)**: Pack-selectable practice mode, $k$-anonymous opening clues ($k_{\min} \ge 3$), exact Tanimoto %, live shared functional group chips, `TPSA` & `pI` strips, and free progressive hint/scaffold unlocks.
- **Build-Time $k$-Anonymity & Solver-Bot Calibration Gates (`tools/solver_bot.py`)**:
  - Level-0 opening clues are automatically selected and verified at build time so every equivalence class has size $\ge 4$ in `ranked` ($k_{\min}=10$ Chemle, $k_{\min}=8$ Seqle) and $\ge 3$ in `learn` ($k_{\min}=4$ Chemle, $k_{\min}=3$ Seqle). Pre-guess UI never exposes exact `MW`, exact `Category`/`family`, `difficulty_tier`, or `information_gain_bits`.
  - Greedy max-entropy solver (`tools/solver_bot.py`) enforces build-time difficulty gates:
    - **Ranked**: mean guesses $= 3.89 \in [3.5, 4.5]$, median $= 4.0$, $\le 2$-guess solves $= 2.78\% < 5\%$, $\le 6$-guess solves $= 100\% \ge 98\%$, minimum opening expected remaining $= 7.56 \ge 2.5$.
    - **Learn**: mean guesses $= 3.17 \le 3.5$, median $= 3.0$, $\le 2$-guess solves $= 19.44\% < 25\%$, $\le 6$-guess solves $= 100\%$.
- **Strict Separation of `guessPool` and `answerPool`**:
  - `guessPool` is **always** the full 36-molecule / 30-motif catalog so a puzzle is never locked out by a category filter.
  - In `ranked` Daily mode, `answerPool` is the full eligible catalog (`36` molecules in Chemle; `28` native/fragment motifs in Seqle) and category filtering is disabled with an explanatory tooltip.
- **Audited Seqle Motif Schema (`kind`, `consensus`, `parent_sequence`, `window_offset`, `terminal`, `exemplar_source`, `game_eligibility`)**:
  - Every motif is classified as `"exact"`, `"fragment"` (a 5-mer window of a longer biological motif), or `"abstraction"` (Learn-only compressed teaching representation) and validated at build time against its parent sequence, terminal placement, and PROSITE/ELM/UniProt/PDB citation.
- **Cross-Language Golden Vector Parity (`tests/golden_feedback.json` & `tests/verify_golden_and_dom.mjs`)**:
  - All **2,592 Chemle** (`36 × 36 × 2`) and **1,800 Seqle** (`30 × 30 × 2`) feedback signatures are verified byte-for-byte between Python (`chem_engine.py`) and browser JavaScript (`app.js`).
- **True ACS-1996 SVG Rendering & Non-Overlapping Carbonyl SMARTS**:
  - Calls `rdMolDraw2D.SetACS1996Mode(opts, Draw.MeanBondLength(mol))` with normalized 1-decimal SVG coordinates and per-record `svg_sha256`.
  - Splits generic carbonyl detection into non-overlapping SMARTS for `Aldehyde (-CHO)`, `Ketone (C-C(=O)-C)`, `Carboxylic Acid (-COOH)`, `Ester / Lactone`, and `Amide / Lactam`.

---

## Quick Start (Zero-Install Browser Play)

Open [`index.html`](index.html) or the self-contained single-file bundle [`dist/biowordle_standalone.html`](dist/biowordle_standalone.html) directly in any modern browser. No Python server or RDKit installation is required to play.

Optional URL query parameters:
- `?seed=CHEM101-WEEK3` — Classroom deterministic puzzle seed across all students.
- `?profile=learn` — Launch directly into Learn Mode.

---

## Rebuilding & Running the Verification Suite

```bash
# 1. Create and activate a Python 3.12 virtual environment
python -m venv .venv
# Windows PowerShell:
.\.venv\Scripts\Activate.ps1
# macOS / Linux:
# source .venv/bin/activate

# 2. Install pinned dependencies
pip install -r requirements.txt

# 3. Compile puzzles, SVGs, golden vectors, and standalone bundle
python build_puzzles.py

# 4. Run the standalone Solver-Bot calibration harness
python tools/solver_bot.py

# 5. Run Cross-Language Golden Vector & Runtime Invariant checks (Node.js)
node tests/verify_golden_and_dom.mjs

# 6. Run the full 11-test Python + Node audit suite
python -m unittest discover -v -s . -p "test_*.py"

# 7. Run headless Chromium / Edge End-to-End browser tests (uses Antigravity playwright-core)
node tests/test_playwright_e2e.mjs
```

---

## License

Released under the [MIT License](LICENSE). Pre-compiled chemical descriptors and 2D depictions are generated with [RDKit](https://www.rdkit.org/) (BSD 3-Clause License) and Peter Ertl's `sascorer` algorithm (BSD 3-Clause License).
