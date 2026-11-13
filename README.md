# BioWordle (`Chemle` + `Seqle`) — Audited RDKit Edition

An interactive cheminformatics and sequence-biophysics deduction suite where **RDKit** validates molecular identities, stereochemistry, chiral circular fingerprints, progressive scaffold depictions, peptide Henderson-Hasselbalch properties, and Shannon-entropy information gain at build time, compiling them into a zero-runtime-dependency browser application.

---

## Architecture & Scientific Specifications

1. **100% Independent Build Compiler (`chem_engine.py` + `build_puzzles.py`)**
   - Uses only public RDKit APIs and the Python standard library (zero external plugin/skill path dependencies).
   - Deterministic dataset version hashing (`SHA-256` over sorted canonical InChIKeys and motif sequences) paired with UTC calendar dates (`YYYY-MM-DD`) so daily challenges are identical across timezones.

2. **Molecular Identity vs. Chiral Morgan Fingerprint Similarity (`Chemle`)**
   - **Standardization Policy:** `rdMolStandardize.Cleanup` + `FragmentParent` + explicit `AssignStereochemistry` without destructive tautomer mutation, preserving curated keto/enol states and $(R)/(S)$ or $(E)/(Z)$ stereocenters. Both `curated_smiles` / `curated_inchikey` and `canonical_smiles` / `inchikey` are recorded in `data/puzzles.json`.
   - **Chiral ECFP4-like Fingerprint Spec:**
     ```python
     rdFingerprintGenerator.GetMorganGenerator(
         radius=2,
         fpSize=2048,
         includeChirality=True,
         useBondTypes=True,
     )
     ```
     Pairwise similarity uses the Tanimoto coefficient $T(A,B) = \frac{|A \cap B|}{|A \cup B|}$ over shared circular substructural features, while win detection strictly checks molecular identity (`inchikey`).
   - **Acyclic Scaffold Fallback:** Ring-containing molecules use their Bemis-Murcko framework (`MurckoScaffold.GetScaffoldForMol`) for intermediate visual clues; acyclic compounds (e.g., Citric Acid) fall back to an unfunctionalized carbon-chain skeleton (`scaffold_mode: "carbon_skeleton"`) so intermediate clues never leak the full functional structure.
   - **Post-Game Molecular Insights & Structure Explorer:** Advanced medicinal chemistry indicators (`Heuristic CNS-MPO Desirability` using Wager 2010 monotonic ramps and `PAINS / Brenk` screening catalogs) are presented in a dedicated post-game educational panel alongside positional-isomer notes (e.g., Protocatechuic Acid `3,4-DHB` vs. Gentisic Acid `2,5-DHB`) and a closest-probe functional group comparison.

3. **Build-Time Shannon Entropy Information Gain Evaluator**
   - For every molecule and peptide motif, `chem_engine.py` simulates all $N \times N$ feedback signatures $f(g, t)$ and computes the expected Shannon information gain in bits:
     $$IG(g) = H(T) - \mathbb{E}_f[H(T \mid f)]$$
   - Displays each puzzle's difficulty tier, nearest chemical-space neighbor, opening-probe $IG(g)$ ranking, and live candidate-pool reduction after every guess.

4. **Seqle Two-Pass Counting, Primary vs. Secondary AA Traits & Peptide $pK_a$ Model**
   - **Precedence & Multiplicity:** Strictly evaluates `exact` matches first, `present` (misplaced) matches second via a two-pass residue counter, and positional `class` matches third (`exact > present > class > absent`).
   - **Primary Class vs. Overlapping Secondary Traits:** Each of the 20 standard amino acids has one deterministic **Primary Class** (`Aliphatic`, `Aromatic`, `Polar`, `Positive`, `Negative`, `Special`) for gameplay rules plus explicit **Secondary Chemical Traits** (e.g., Histidine: `Positive (Imidazole, pKa 6.0) · Aromatic (Heteroaryl) · Metal Ligand`; Tyrosine: `Aromatic · Polar (Phenolic -OH)`).
   - **Motif Provenance & Context:** Every 5-residue puzzle stores its `pattern`, `context_requirement` (structural fold or terminal topology required for biological activity), and experimental `exemplar_source` (UniProt/PDB reference).
   - **Documented Peptide $pK_a$ Model:** Net charge at $\text{pH } 7.4$ and isoelectric point ($pI$) are computed via the Henderson-Hasselbalch equation for a free pentapeptide ($\text{N-terminus } pK_a = 8.0$, $\text{C-terminus } pK_a = 3.1$, $\text{Asp } 3.9$, $\text{Glu } 4.1$, $\text{His } 6.0$, $\text{Cys } 8.3$, $\text{Tyr } 10.1$, $\text{Lys } 10.5$, $\text{Arg } 12.5$).

---

## Rebuilding & Running the Audit Test Suite

```powershell
# 1. Compile the puzzle database and standalone HTML bundle
& ".\.venv\Scripts\python.exe" build_puzzles.py

# 2. Run the 9-test scientific & engineering verification suite
& ".\.venv\Scripts\python.exe" -m unittest discover -s . -p "test_*.py" -v
```
