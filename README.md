# BioWordle v0.2 ("Restraint") — Chemle & Seqle (RDKit Edition)

**BioWordle** is a two-layer interactive cheminformatics and structural biology deduction suite comprising:
1. **🧪 Chemle**: Deduce a target metabolite, hydrocarbon degradation intermediate, natural product, or therapeutic drug (36 curated molecules) from **Chiral Morgan (`ECFP4`-like) Tanimoto similarity**, physicochemical directional deltas, functional group overlap, and progressive **Bemis–Murcko / Carbon-Skeleton 2D scaffolds**.
2. **🧬 Seqle**: Decode a 5-residue functional protein motif (`EKDEL`, `DEVDG`, `GRGDS`, `LPETG`, `HEIGH`, `KCVIM`, `DEADE`, `GVGKS`, `VDFGL`, `HRDLA`, `WSEWS`, `PPPVP`, `INFEK`, `YQRLN`, `KKMPQ`, `CAACH`, `NLTGY`, `DSGLD`, `KENIP`, `KKKRK`, `FGGFT`, `YGGFM`, `SYSGY`, `NPLYR`, `GHHPH`, plus Learn-only pedagogical abstractions `CCHHC` and `LDLAS`) using two-pass multiplicity-aware residue feedback and **Henderson–Hasselbalch** net charge at pH 7.4, Kyte–Doolittle GRAVY, and isoelectric point ($\text{pI}$).

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
```

---

## License

Released under the [MIT License](LICENSE). Pre-compiled chemical descriptors and 2D depictions are generated with [RDKit](https://www.rdkit.org/) (BSD 3-Clause License) and Peter Ertl's `sascorer` algorithm (BSD 3-Clause License).
