# BioWordle (`Chemle` + `Seqle`) — RDKit Build-Time Edition

A zero-runtime-dependency dual-mode biochemistry and cheminformatics puzzle game where **RDKit** compiles and validates all chemical structures, Bemis-Murcko scaffolds, descriptors, Morgan fingerprint similarities, and pentapeptide 2D depictions at build time into a static browser bundle.

---

## Game Modes

### 1. 🧪 Chemle (Small Molecules, Metabolites & Hydrocarbon Intermediates)
- **36 Curated Compounds** across 4 categories:
  - *Hydrocarbon Degradation & Microbial Metabolites* (Catechol, Salicylic Acid, Protocatechuic Acid, Gentisic Acid, Naphthalene, Phenanthrene, Camphor, Pyocyanin, Prodigiosin, Violacein, Geosmin)
  - *Core Metabolites, Cofactors & Neurotransmitters* (Caffeine, Dopamine, Serotonin, Adrenaline, Melatonin, Histamine, L-Tryptophan, Ascorbic Acid, Biotin, Citric Acid, Shikimic Acid, Cholesterol)
  - *Plant Natural Products & Alkaloids* (Vanillin, Resveratrol, Curcumin, Capsaicin, Nicotine, Quinine, Artemisinin)
  - *Classic Therapeutics* (Aspirin, Ibuprofen, Acetaminophen, Penicillin G, Ciprofloxacin, Fluoxetine)
- **RDKit Build-Time Features**:
  - Canonical descriptors (`MW`, `cLogP`, `TPSA`, `HBD`, `HBA`, `RotB`, `Fsp3`, `HeavyAtoms`, `Rings`, `AromaticRings`, `QED`, `CNS_MPO`, `SA_Score`) + PAINS/Brenk alerts.
  - Full $36 \times 36$ pairwise **Morgan Fingerprint Tanimoto Similarity (%)** matrix.
  - SMARTS-verified functional group badges highlighting shared substructures with the hidden target.
  - Progressive 2D ACS-1996 SVG structure clues (Silhouette $\rightarrow$ Bemis-Murcko Scaffold Core $\rightarrow$ Full 2D Structure).

### 2. 🧬 Seqle (5-Residue Functional Protein Motifs)
- **30 Canonical 5-AA Motifs** (`GHSMG`, `KDELR`, `DEVDG`, `RGDVF`, `LPETG`, `HEIGH`, `CCHHC`, `YIGSR`, `NPVYK`, `CVIMK`, `DEADR`, `GAGTG`, `VDFGL`, `HRDLA`, `WSEWS`, `PPPVP`, `INFEK`, `YQRLN`, `KKMPQ`, `CAACH`, `NLTGY`, `DSGLS`, `KENKR`, `PKKKR`, `FGGFD`, `YGGFM`, `SYGYS`, `LDLAS`, `NPLYR`, `HPHPH`).
- **3-Tier Tile Feedback**:
  - 🟩 **Emerald (`Exact`)**: Right amino acid in the exact position.
  - 🟨 **Amber (`Present`)**: Amino acid is in the motif, different position.
  - 🟦 **Cyan (`Biochem Class Match`)**: Wrong exact amino acid, *but* the target residue at that position belongs to the same biochemical class (`Aliphatic`, `Aromatic`, `Polar`, `Positive`, `Negative`, `Special`).
- **Live Peptide Property Strip**: Compares Net Charge at pH 7.4, Kyte-Doolittle GRAVY hydropathy, and Isoelectric Point ($pI$) for every 5-mer guess against the target, plus an RDKit-rendered 2D pentapeptide chemical structure (`Chem.MolFromSequence`).

---

## Quick Start

### Play in Browser (Zero Runtime Dependencies)
Open `index.html` or `dist/biowordle_standalone.html` directly in any modern browser.

### Rebuild Puzzle Database After Adding New Molecules or Motifs
Edit `CHEMLE_RAW_COMPOUNDS` or `SEQLE_RAW_MOTIFS` in `build_puzzles.py`, then run:

```powershell
& "d:\W\Play W\biowordle\.venv\Scripts\python.exe" "d:\W\Play W\biowordle\build_puzzles.py"
```

