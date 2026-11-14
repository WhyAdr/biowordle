# Chemle v0.3 Molecular Audit

Offline and online structural provenance, stereochemical/protonation audit, close-pair similarity analysis, and pool-separation verification for BioWordle v0.3.

## 1. Scope & Pool Architecture

- **Canonical Manifest:** `data/chemle_compounds.json` (`schema_version`: `0.3.0`)
- **Total Catalog Compounds:** `78` (`71` `verified`, `7` `baseline_legacy`)
- **Ranked Pool (`ranked`):** `36` answers / `36` legal guesses (`chemle_schedule_hash`: `9257bd168fa6`, unchanged from v0.2.2)
- **Learn Pilot Pool (`learn`):** `78` answers / `78` legal guesses (`36` Ranked baseline + `42` Learn-only additions)
- **Seqle Baseline Parity:** `28` Ranked / `30` Learn motifs (`seqle_schedule_hash`: `623f02623c18`, unchanged)
- **Dataset Canonical Hash:** `488378e66071`
- **RDKit Compiler Version:** `2026.09.2` (`Morgan ECFP4 radius=2, 2048 bits`)

All `78` manifest compounds pass `chem_engine.validate_chemle_manifest()` offline (SMILES validity, single-component enforcement, parent uncharged InChIKey deduplication, formula consistency, defined stereo/E-Z verification, racemate check) and `tools/verify_chemle_manifest_online.py` online against PubChem PUG REST and EMBL-EBI OLS4 ChEBI (`78/78` resolved).

## 2. Baseline Ranked Pool Audit (`36` Compounds)

The `36` Ranked-eligible compounds are preserved byte-identical in canonical SMILES and descriptor/SVG hashes so Ranked schedules and difficulty calibrations remain completely isolated from Learn expansion.

- **`verified` (`29` compounds):** Match canonical neutral PubChem and ChEBI stereochemistry and connectivity. Includes two clinically marketed racemates (`fluoxetine`, CID `3386` / `CHEBI:5118`; `ibuprofen`, CID `3672` / `CHEBI:5855`) explicitly tagged `stereo_status="racemate"`.
- **`baseline_legacy` (`7` compounds):** Explicitly flagged in `data/chemle_compounds.json` and surfaced post-reveal in the UI provenance card rather than altering historical Ranked hashes:

| ID | Name | Formula | PubChem CID | ChEBI | Stereo / Protonation | InChIKey | Audit Note |
|---|---|---|---|---|---|---|---|
| `artemisinin` | Artemisinin | `C15H22O5` | `138059745` | `CHEBI:223316` | `defined_stereo` / `neutral_canonical` | `BLUAFEHZUWYNDE-GZEPNZDPSA-N` | Preserved v0.2.2 Ranked SMILES (matches PubChem CID 138059745 enantiomer/stereoisomer rather than canonical (+)-artemisinin CID 68827) to maintain schedule hash parity. |
| `cholesterol` | Cholesterol | `C27H46O` | `5997` | `CHEBI:16113` | `defined_stereo` / `neutral_canonical` | `HVYWMOMLDIMFJA-PZIGXULQSA-N` | Preserved v0.2.2 Ranked SMILES (inverted configuration at one chiral center vs PubChem CID 5997 HVYWMOMLDIMFJA-DPAQBDIFSA-N) to maintain schedule hash parity. |
| `geosmin` | Geosmin | `C12H22O` | `29746` | `CHEBI:46703` | `unspecified_legacy` / `neutral_canonical` | `SQFNUAVVZSZPMP-XMACVGRKSA-N` | Preserved v0.2.2 Ranked SMILES (legacy connectivity and 1 unspecified chiral center vs PubChem CID 29746 JLPUXFOGCDVKGO-TUAOUCFPSA-N) to maintain schedule hash parity. |
| `prodigiosin` | Prodigiosin | `C19H23N3O` | `172897611` | `CHEBI:82758` | `unspecified_legacy` / `neutral_canonical` | `LCNZNMPZJFZRAP-UHFFFAOYSA-N` | Preserved v0.2.2 Ranked SMILES (C4 sidechain C19H23N3O and unspecified methine double-bond stereo vs canonical C20H25N3O PubChem CID 172897611) to maintain schedule hash parity. |
| `pyocyanin` | Pyocyanin | `C13H10N2O` | `14360019` | `CHEBI:8653` | `achiral` / `neutral_major_tautomer` | `HAQCALIJITZJCC-UHFFFAOYSA-N` | Preserved v0.2.2 neutral keto-tautomer form (10-methylphenazin-2(10H)-one, PubChem CID 14360019) rather than zwitterionic betaine resonance contributor to maintain schedule hash parity. |
| `quinine` | Quinine | `C20H24N2O2` | `94175` | `CHEBI:15854` | `unspecified_legacy` / `neutral_canonical` | `LOUPRKONTZGTKE-AFHBHXEDSA-N` | Preserved v0.2.2 Ranked SMILES (matches PubChem CID 94175 epiquinidine / unspecified bridgehead nitrogen rather than canonical Quinine CID 3034034) to maintain schedule hash parity. |
| `violacein` | Violacein | `C28H18N4O3` | `11053` | `CHEBI:131914` | `unspecified_legacy` / `neutral_canonical` | `PJAOXHSAWDWRMO-UHFFFAOYSA-N` | Preserved v0.2.2 Ranked SMILES (extra indole substituent C28H18N4O3 and unspecified double bond vs canonical C20H13N3O3 PubChem CID 11053) to maintain schedule hash parity. |

## 3. Learn-Only Pilot Expansion (`42` Compounds Across 7 Biochemical Families)

All `42` additions carry `game_eligibility=["learn"]`, `guess_eligibility=["learn"]`, and `review_status="verified"`. Special curation notes verified during structural audit:
- **`sinapic_acid` (Trans `(E)`-sinapic acid):** Corrected to `(E)` double-bond stereoisomer (`CID:637775`, `CHEBI:15714`, InChIKey `PCMORTLOPMLEFB-ONEGZZNKSA-N`) matching its natural phenylpropanoid configuration.
- **`pyrrolnitrin` (Microbial dichlorophenylpyrrole-nitro):** Dipolar aromatic nitro group `[N+](=O)[O-]` verified neutral overall (`protonation_policy="neutral_canonical"`) via `_NITRO_GROUP_SMARTS` exception in `chem_engine.py`.
- **`ketoprofen` (Therapeutic NSAID):** Clinically marketed `(RS)`-racemate (`CID:3825`, `CHEBI:6128`, `stereo_status="racemate"`, `ChiralCenters=1`).

### 3.1 Microbial and related natural products (`6` compounds)

| ID | Name | Category | Formula | PubChem CID | ChEBI | Stereo | MW | cLogP | TPSA | HBD/HBA | HA | Rings (Ar) | InChIKey |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `1_hydroxyphenazine` | 1-Hydroxyphenazine | Microbial Natural Product | `C12H8N2O` | `135412648` | `CHEBI:62216` | `achiral` | 196.2 | 2.49 | 46.0 | 1/3 | 15 | 3 (3) | `SVRNCBGWUMMBQB-UHFFFAOYSA-N` |
| `2_4_diacetylphloroglucinol` | 2,4-Diacetylphloroglucinol | Microbial Natural Product | `C10H10O5` | `16547` | `CHEBI:78688` | `achiral` | 210.2 | 1.21 | 94.8 | 3/5 | 15 | 1 (1) | `PIFFQYJYNWXNGE-UHFFFAOYSA-N` |
| `indole_3_acetic_acid` | Indole-3-acetic acid | Microbial Natural Product | `C10H9NO2` | `802` | `CHEBI:16411` | `achiral` | 175.2 | 1.79 | 53.1 | 2/1 | 13 | 2 (2) | `SEOVTRFCIGRIMH-UHFFFAOYSA-N` |
| `phenazine` | Phenazine | Microbial Natural Product | `C12H8N2` | `4757` | `CHEBI:36674` | `achiral` | 180.2 | 2.78 | 25.8 | 0/2 | 14 | 3 (3) | `PCNDJXKNXGMECE-UHFFFAOYSA-N` |
| `phenazine_1_carboxylic_acid` | Phenazine-1-carboxylic acid | Microbial Natural Product | `C13H8N2O2` | `95069` | `CHEBI:62412` | `achiral` | 224.2 | 2.48 | 63.1 | 1/3 | 17 | 3 (3) | `JGCSKOVQDXEQHI-UHFFFAOYSA-N` |
| `pyrrolnitrin` | Pyrrolnitrin | Microbial Natural Product | `C10H6Cl2N2O2` | `13916` | `CHEBI:32079` | `achiral` | 257.1 | 3.9 | 58.9 | 1/2 | 16 | 2 (2) | `QJBZDBLBQWFTPZ-UHFFFAOYSA-N` |

### 3.2 Aromatic degradation (`6` compounds)

| ID | Name | Category | Formula | PubChem CID | ChEBI | Stereo | MW | cLogP | TPSA | HBD/HBA | HA | Rings (Ar) | InChIKey |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `4_hydroxybenzoic_acid` | 4-Hydroxybenzoic acid | Hydrocarbon Degradation | `C7H6O3` | `135` | `CHEBI:30763` | `achiral` | 138.1 | 1.09 | 57.5 | 2/2 | 10 | 1 (1) | `FJKROLUGYXJWQN-UHFFFAOYSA-N` |
| `benzoic_acid` | Benzoic acid | Hydrocarbon Degradation | `C7H6O2` | `243` | `CHEBI:30746` | `achiral` | 122.1 | 1.38 | 37.3 | 1/1 | 9 | 1 (1) | `WPYMKLBDIGXBTP-UHFFFAOYSA-N` |
| `m_xylene` | m-Xylene | Hydrocarbon Degradation | `C8H10` | `7929` | `CHEBI:28488` | `achiral` | 106.2 | 2.3 | 0.0 | 0/0 | 8 | 1 (1) | `IVSZLXZYQVIEFR-UHFFFAOYSA-N` |
| `o_xylene` | o-Xylene | Hydrocarbon Degradation | `C8H10` | `7237` | `CHEBI:28063` | `achiral` | 106.2 | 2.3 | 0.0 | 0/0 | 8 | 1 (1) | `CTQNGGLPUBDAKN-UHFFFAOYSA-N` |
| `p_xylene` | p-Xylene | Hydrocarbon Degradation | `C8H10` | `7809` | `CHEBI:27417` | `achiral` | 106.2 | 2.3 | 0.0 | 0/0 | 8 | 1 (1) | `URLKBWYHVLBVBO-UHFFFAOYSA-N` |
| `phenol` | Phenol | Hydrocarbon Degradation | `C6H6O` | `996` | `CHEBI:15882` | `achiral` | 94.1 | 1.39 | 20.2 | 1/1 | 7 | 1 (1) | `ISWSIDIOOBJBQZ-UHFFFAOYSA-N` |

### 3.3 Neuroactive metabolites (`6` compounds)

| ID | Name | Category | Formula | PubChem CID | ChEBI | Stereo | MW | cLogP | TPSA | HBD/HBA | HA | Rings (Ar) | InChIKey |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `5_htp` | 5-Hydroxy-L-tryptophan | Neurotransmitter | `C11H12N2O3` | `439280` | `CHEBI:17780` | `defined_stereo` | 220.2 | 0.83 | 99.3 | 4/3 | 16 | 2 (2) | `LDCYZAJDBXYCGN-VIFPVBQESA-N` |
| `gaba` | gamma-Aminobutyric acid | Neurotransmitter | `C4H9NO2` | `119` | `CHEBI:16865` | `achiral` | 103.1 | -0.19 | 63.3 | 2/2 | 7 | 0 (0) | `BTCSSZJGUNDROE-UHFFFAOYSA-N` |
| `l_dopa` | L-DOPA | Neurotransmitter | `C9H11NO4` | `6047` | `CHEBI:15765` | `defined_stereo` | 197.2 | 0.05 | 103.8 | 4/4 | 14 | 1 (1) | `WTDRDQBEARUVNC-LURJTMIESA-N` |
| `norepinephrine` | Norepinephrine | Hormone / Neurotransmitter | `C8H11NO3` | `439260` | `CHEBI:18357` | `defined_stereo` | 169.2 | 0.09 | 86.7 | 4/4 | 12 | 1 (1) | `SFLSHLFXELFNJZ-QMMMGPOBSA-N` |
| `tryptamine` | Tryptamine | Biogenic Amine | `C10H12N2` | `1150` | `CHEBI:16765` | `achiral` | 160.2 | 1.67 | 41.8 | 2/1 | 12 | 2 (2) | `APJYDQYYACXCRM-UHFFFAOYSA-N` |
| `tyramine` | Tyramine | Biogenic Amine | `C8H11NO` | `5610` | `CHEBI:15760` | `achiral` | 137.2 | 0.89 | 46.2 | 2/2 | 10 | 1 (1) | `DZGWFCGJZKJUFP-UHFFFAOYSA-N` |

### 3.4 Purines and nucleobases (`6` compounds)

| ID | Name | Category | Formula | PubChem CID | ChEBI | Stereo | MW | cLogP | TPSA | HBD/HBA | HA | Rings (Ar) | InChIKey |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `adenine` | Adenine | Primary Metabolite | `C5H5N5` | `190` | `CHEBI:16708` | `achiral` | 135.1 | -0.06 | 80.5 | 2/4 | 10 | 2 (2) | `GFFGJBXGBJISGV-UHFFFAOYSA-N` |
| `hypoxanthine` | Hypoxanthine | Primary Metabolite | `C5H4N4O` | `135398638` | `CHEBI:17368` | `achiral` | 136.1 | -0.35 | 74.4 | 2/3 | 10 | 2 (2) | `FDGQSTZJBFJUBT-UHFFFAOYSA-N` |
| `theobromine` | Theobromine | Alkaloid / Metabolite | `C7H8N4O2` | `5429` | `CHEBI:28946` | `achiral` | 180.2 | -1.04 | 72.7 | 1/3 | 13 | 2 (2) | `YAPQBXQYLJRXSA-UHFFFAOYSA-N` |
| `theophylline` | Theophylline | Alkaloid / Metabolite | `C7H8N4O2` | `2153` | `CHEBI:28177` | `achiral` | 180.2 | -1.04 | 72.7 | 1/3 | 13 | 2 (2) | `ZFXYFBGIUFBOJW-UHFFFAOYSA-N` |
| `uracil` | Uracil | Primary Metabolite | `C4H4N2O2` | `1174` | `CHEBI:17568` | `achiral` | 112.1 | -0.94 | 65.7 | 2/2 | 8 | 1 (1) | `ISAKRJDGNUQOIC-UHFFFAOYSA-N` |
| `xanthine` | Xanthine | Primary Metabolite | `C5H4N4O2` | `1188` | `CHEBI:15318` | `achiral` | 152.1 | -1.06 | 94.4 | 3/3 | 11 | 2 (2) | `LRFVTYWOQMYALW-UHFFFAOYSA-N` |

### 3.5 Therapeutic analogues (`6` compounds)

| ID | Name | Category | Formula | PubChem CID | ChEBI | Stereo | MW | cLogP | TPSA | HBD/HBA | HA | Rings (Ar) | InChIKey |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `amoxicillin` | Amoxicillin | Therapeutic Drug | `C16H19N3O5S` | `33613` | `CHEBI:2676` | `defined_stereo` | 365.4 | 0.02 | 133.0 | 4/6 | 25 | 3 (1) | `LSQZJLSUYDQPKJ-NJBDSQKTSA-N` |
| `diclofenac` | Diclofenac | Therapeutic Drug | `C14H11Cl2NO2` | `3033` | `CHEBI:47381` | `achiral` | 296.2 | 4.36 | 49.3 | 2/2 | 19 | 2 (2) | `DCOPUUMXTXDBNB-UHFFFAOYSA-N` |
| `ketoprofen` | Ketoprofen | Therapeutic Drug | `C16H14O3` | `3825` | `CHEBI:6128` | `racemate` | 254.3 | 3.11 | 54.4 | 1/2 | 19 | 2 (2) | `DKYWVDODHFEZIM-UHFFFAOYSA-N` |
| `mefenamic_acid` | Mefenamic acid | Therapeutic Drug | `C15H15NO2` | `4044` | `CHEBI:6717` | `achiral` | 241.3 | 3.75 | 49.3 | 2/2 | 18 | 2 (2) | `HYYBABOKPJLUIN-UHFFFAOYSA-N` |
| `naproxen` | Naproxen | Therapeutic Drug | `C14H14O3` | `156391` | `CHEBI:7476` | `defined_stereo` | 230.3 | 3.04 | 46.5 | 1/2 | 17 | 2 (2) | `CMWTZPSULFXXJA-VIFPVBQESA-N` |
| `norfloxacin` | Norfloxacin | Therapeutic Drug | `C16H18FN3O3` | `4539` | `CHEBI:100246` | `achiral` | 319.3 | 1.27 | 74.6 | 2/4 | 23 | 3 (2) | `OGJPXUAPXNRGGI-UHFFFAOYSA-N` |

### 3.6 Phenylpropanoids and plant metabolites (`6` compounds)

| ID | Name | Category | Formula | PubChem CID | ChEBI | Stereo | MW | cLogP | TPSA | HBD/HBA | HA | Rings (Ar) | InChIKey |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `caffeic_acid` | Caffeic acid | Plant Natural Product | `C9H8O4` | `689043` | `CHEBI:36281` | `defined_stereo` | 180.2 | 1.2 | 77.8 | 3/3 | 13 | 1 (1) | `QAIPRVGONGVQAS-DUXPYHPUSA-N` |
| `eugenol` | Eugenol | Plant Natural Product | `C10H12O2` | `3314` | `CHEBI:4917` | `achiral` | 164.2 | 2.13 | 29.5 | 1/2 | 12 | 1 (1) | `RRAFCDWBNXTKKO-UHFFFAOYSA-N` |
| `ferulic_acid` | Ferulic acid | Plant Natural Product | `C10H10O4` | `445858` | `CHEBI:193350` | `defined_stereo` | 194.2 | 1.5 | 66.8 | 2/3 | 14 | 1 (1) | `KSEBMYQBYZTDHS-HWKANZROSA-N` |
| `p_coumaric_acid` | p-Coumaric acid | Plant Natural Product | `C9H8O3` | `637542` | `CHEBI:32374` | `defined_stereo` | 164.2 | 1.49 | 57.5 | 2/2 | 12 | 1 (1) | `NGSWKAQJJWESNS-ZZXKWVIFSA-N` |
| `sinapic_acid` | Sinapic acid | Plant Natural Product | `C11H12O5` | `637775` | `CHEBI:15714` | `defined_stereo` | 224.2 | 1.51 | 76.0 | 2/4 | 16 | 1 (1) | `PCMORTLOPMLEFB-ONEGZZNKSA-N` |
| `vanillic_acid` | Vanillic acid | Plant Natural Product | `C8H8O4` | `8468` | `CHEBI:30816` | `achiral` | 168.1 | 1.1 | 66.8 | 2/3 | 12 | 1 (1) | `WKOLLVMJNQIZCI-UHFFFAOYSA-N` |

### 3.7 Primary metabolism (`6` compounds)

| ID | Name | Category | Formula | PubChem CID | ChEBI | Stereo | MW | cLogP | TPSA | HBD/HBA | HA | Rings (Ar) | InChIKey |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `fumaric_acid` | Fumaric acid | Primary Metabolite | `C4H4O4` | `444972` | `CHEBI:18012` | `defined_stereo` | 116.1 | -0.29 | 74.6 | 2/2 | 8 | 0 (0) | `VZCYOOQTPOCHFL-OWOJBTEDSA-N` |
| `l_malic_acid` | L-Malic acid | Primary Metabolite | `C4H6O5` | `222656` | `CHEBI:30797` | `defined_stereo` | 134.1 | -1.09 | 94.8 | 3/3 | 9 | 0 (0) | `BJEPYKJPYRNKOW-REOHCLBHSA-N` |
| `l_phenylalanine` | L-Phenylalanine | Primary Metabolite | `C9H11NO2` | `6140` | `CHEBI:17295` | `defined_stereo` | 165.2 | 0.64 | 63.3 | 2/2 | 12 | 1 (1) | `COLNVLDHVKWLRT-QMMMGPOBSA-N` |
| `l_tyrosine` | L-Tyrosine | Primary Metabolite | `C9H11NO3` | `6057` | `CHEBI:17895` | `defined_stereo` | 181.2 | 0.35 | 83.5 | 3/3 | 13 | 1 (1) | `OUYCCCASQSFEME-QMMMGPOBSA-N` |
| `pyruvic_acid` | Pyruvic acid | Primary Metabolite | `C3H4O3` | `1060` | `CHEBI:32816` | `achiral` | 88.1 | -0.34 | 54.4 | 1/2 | 6 | 0 (0) | `LCTONWCANYUPML-UHFFFAOYSA-N` |
| `succinic_acid` | Succinic acid | Primary Metabolite | `C4H6O4` | `1110` | `CHEBI:15741` | `achiral` | 118.1 | -0.06 | 74.6 | 2/2 | 8 | 0 (0) | `KDYFGRWQOYBRFD-UHFFFAOYSA-N` |

## 4. Close-Pair Similarity & Regioisomer/Congener Separation

Morgan ECFP4 (radius `2`, `2048` bits, rounded integer %) pairwise Tanimoto similarities and descriptor vectors demonstrate that Learn feedback separates even constitutional isomers and congeneric series:

### Methylxanthines

| Compound | Formula | MW | cLogP | TPSA | HBD/HBA | HeavyAtoms | RotB | ArRings | Functional Groups |
|---|---|---|---|---|---|---|---|---|---|
| Caffeine (`caffeine`) | `C8H10N4O2` | 194.2 | -1.03 | 61.8 | 0/3 | 14 | 0 | 2 | Heteroaromatic N |
| Theophylline (`theophylline`) | `C7H8N4O2` | 180.2 | -1.04 | 72.7 | 1/3 | 13 | 0 | 2 | Heteroaromatic N |
| Theobromine (`theobromine`) | `C7H8N4O2` | 180.2 | -1.04 | 72.7 | 1/3 | 13 | 0 | 2 | Heteroaromatic N |

- **Pairwise Morgan ECFP4 Tanimoto:** `caffeine` vs `theophylline` = `46%`; `caffeine` vs `theobromine` = `53%`; `theophylline` vs `theobromine` = `29%`

### Hydroxycinnamates

| Compound | Formula | MW | cLogP | TPSA | HBD/HBA | HeavyAtoms | RotB | ArRings | Functional Groups |
|---|---|---|---|---|---|---|---|---|---|
| p-Coumaric acid (`p_coumaric_acid`) | `C9H8O3` | 164.2 | 1.49 | 57.5 | 2/2 | 12 | 2 | 1 | Phenol / Enol (-OH), Carboxylic Acid (-COOH) |
| Caffeic acid (`caffeic_acid`) | `C9H8O4` | 180.2 | 1.2 | 77.8 | 3/3 | 13 | 2 | 1 | Phenol / Enol (-OH), Carboxylic Acid (-COOH) |
| Ferulic acid (`ferulic_acid`) | `C10H10O4` | 194.2 | 1.5 | 66.8 | 2/3 | 14 | 3 | 1 | Phenol / Enol (-OH), Carboxylic Acid (-COOH), Ether / Methoxy (C-O-C) |
| Sinapic acid (`sinapic_acid`) | `C11H12O5` | 224.2 | 1.51 | 76.0 | 2/4 | 16 | 4 | 1 | Phenol / Enol (-OH), Carboxylic Acid (-COOH), Ether / Methoxy (C-O-C) |

- **Pairwise Morgan ECFP4 Tanimoto:** `p_coumaric_acid` vs `caffeic_acid` = `61%`; `p_coumaric_acid` vs `ferulic_acid` = `50%`; `p_coumaric_acid` vs `sinapic_acid` = `45%`; `caffeic_acid` vs `ferulic_acid` = `64%`; `caffeic_acid` vs `sinapic_acid` = `46%`; `ferulic_acid` vs `sinapic_acid` = `68%`

### Catecholamines

| Compound | Formula | MW | cLogP | TPSA | HBD/HBA | HeavyAtoms | RotB | ArRings | Functional Groups |
|---|---|---|---|---|---|---|---|---|---|
| Adrenaline (`adrenaline`) | `C9H13NO3` | 183.2 | 0.35 | 72.7 | 4/4 | 13 | 3 | 1 | Aliphatic Alcohol (-OH), Phenol / Enol (-OH), Primary / Secondary Amine |
| Dopamine (`dopamine`) | `C8H11NO2` | 153.2 | 0.6 | 66.5 | 3/3 | 11 | 2 | 1 | Phenol / Enol (-OH), Primary / Secondary Amine |
| Norepinephrine (`norepinephrine`) | `C8H11NO3` | 169.2 | 0.09 | 86.7 | 4/4 | 12 | 2 | 1 | Aliphatic Alcohol (-OH), Phenol / Enol (-OH), Primary / Secondary Amine |
| L-DOPA (`l_dopa`) | `C9H11NO4` | 197.2 | 0.05 | 103.8 | 4/4 | 14 | 3 | 1 | Phenol / Enol (-OH), Carboxylic Acid (-COOH), Primary / Secondary Amine |

- **Pairwise Morgan ECFP4 Tanimoto:** `adrenaline` vs `dopamine` = `31%`; `adrenaline` vs `norepinephrine` = `58%`; `adrenaline` vs `l_dopa` = `29%`; `dopamine` vs `norepinephrine` = `42%`; `dopamine` vs `l_dopa` = `45%`; `norepinephrine` vs `l_dopa` = `35%`

### Xylenes

| Compound | Formula | MW | cLogP | TPSA | HBD/HBA | HeavyAtoms | RotB | ArRings | Functional Groups |
|---|---|---|---|---|---|---|---|---|---|
| o-Xylene (`o_xylene`) | `C8H10` | 106.2 | 2.3 | 0.0 | 0/0 | 8 | 0 | 1 | None |
| m-Xylene (`m_xylene`) | `C8H10` | 106.2 | 2.3 | 0.0 | 0/0 | 8 | 0 | 1 | None |
| p-Xylene (`p_xylene`) | `C8H10` | 106.2 | 2.3 | 0.0 | 0/0 | 8 | 0 | 1 | None |

- **Pairwise Morgan ECFP4 Tanimoto:** `o_xylene` vs `m_xylene` = `38%`; `o_xylene` vs `p_xylene` = `38%`; `m_xylene` vs `p_xylene` = `43%`

### Phenazines

| Compound | Formula | MW | cLogP | TPSA | HBD/HBA | HeavyAtoms | RotB | ArRings | Functional Groups |
|---|---|---|---|---|---|---|---|---|---|
| Pyocyanin (`pyocyanin`) | `C13H10N2O` | 210.2 | 2.04 | 34.9 | 0/2 | 16 | 0 | 1 | Heteroaromatic N |
| Phenazine-1-carboxylic acid (`phenazine_1_carboxylic_acid`) | `C13H8N2O2` | 224.2 | 2.48 | 63.1 | 1/3 | 17 | 1 | 3 | Phenol / Enol (-OH), Carboxylic Acid (-COOH), Heteroaromatic N |
| 1-Hydroxyphenazine (`1_hydroxyphenazine`) | `C12H8N2O` | 196.2 | 2.49 | 46.0 | 1/3 | 15 | 0 | 3 | Phenol / Enol (-OH), Heteroaromatic N |
| Phenazine (`phenazine`) | `C12H8N2` | 180.2 | 2.78 | 25.8 | 0/2 | 14 | 0 | 3 | Heteroaromatic N |

- **Pairwise Morgan ECFP4 Tanimoto:** `pyocyanin` vs `phenazine_1_carboxylic_acid` = `22%`; `pyocyanin` vs `1_hydroxyphenazine` = `21%`; `pyocyanin` vs `phenazine` = `27%`; `phenazine_1_carboxylic_acid` vs `1_hydroxyphenazine` = `53%`; `phenazine_1_carboxylic_acid` vs `phenazine` = `44%`; `1_hydroxyphenazine` vs `phenazine` = `52%`

### Salicylates & NSAIDs

| Compound | Formula | MW | cLogP | TPSA | HBD/HBA | HeavyAtoms | RotB | ArRings | Functional Groups |
|---|---|---|---|---|---|---|---|---|---|
| Aspirin (`aspirin`) | `C9H8O4` | 180.2 | 1.31 | 63.6 | 1/3 | 13 | 2 | 1 | Phenol / Enol (-OH), Carboxylic Acid (-COOH), Ester / Lactone |
| Salicylic Acid (`salicylic_acid`) | `C7H6O3` | 138.1 | 1.09 | 57.5 | 2/2 | 10 | 1 | 1 | Phenol / Enol (-OH), Carboxylic Acid (-COOH) |
| Ibuprofen (`ibuprofen`) | `C13H18O2` | 206.3 | 3.07 | 37.3 | 1/1 | 15 | 4 | 1 | Phenol / Enol (-OH), Carboxylic Acid (-COOH) |
| Naproxen (`naproxen`) | `C14H14O3` | 230.3 | 3.04 | 46.5 | 1/2 | 17 | 3 | 2 | Phenol / Enol (-OH), Carboxylic Acid (-COOH), Ether / Methoxy (C-O-C) |
| Ketoprofen (`ketoprofen`) | `C16H14O3` | 254.3 | 3.11 | 54.4 | 1/2 | 19 | 4 | 2 | Phenol / Enol (-OH), Carboxylic Acid (-COOH), Ketone (C-C(=O)-C) |
| Diclofenac (`diclofenac`) | `C14H11Cl2NO2` | 296.2 | 4.36 | 49.3 | 2/2 | 19 | 4 | 2 | Phenol / Enol (-OH), Carboxylic Acid (-COOH), Primary / Secondary Amine, Halogen (F/Cl/Br) |
| Mefenamic acid (`mefenamic_acid`) | `C15H15NO2` | 241.3 | 3.75 | 49.3 | 2/2 | 18 | 3 | 2 | Phenol / Enol (-OH), Carboxylic Acid (-COOH), Primary / Secondary Amine |

- **Pairwise Morgan ECFP4 Tanimoto:** `aspirin` vs `salicylic_acid` = `45%`; `aspirin` vs `ibuprofen` = `20%`; `aspirin` vs `naproxen` = `21%`; `aspirin` vs `ketoprofen` = `27%`; `aspirin` vs `diclofenac` = `20%`; `aspirin` vs `mefenamic_acid` = `36%`; `salicylic_acid` vs `ibuprofen` = `19%`; `salicylic_acid` vs `naproxen` = `18%`; `salicylic_acid` vs `ketoprofen` = `24%`; `salicylic_acid` vs `diclofenac` = `22%`; `salicylic_acid` vs `mefenamic_acid` = `38%`; `ibuprofen` vs `naproxen` = `32%`; `ibuprofen` vs `ketoprofen` = `39%`; `ibuprofen` vs `diclofenac` = `17%`; `ibuprofen` vs `mefenamic_acid` = `17%`; `naproxen` vs `ketoprofen` = `30%`; `naproxen` vs `diclofenac` = `13%`; `naproxen` vs `mefenamic_acid` = `19%`; `ketoprofen` vs `diclofenac` = `20%`; `ketoprofen` vs `mefenamic_acid` = `24%`; `diclofenac` vs `mefenamic_acid` = `30%`

## 5. k-Anonymity & Solver Calibration Audit

- **Chemle Ranked (`N=36`):** Level-0 features `['heavy_atom_tercile']`, `k_min = 10`, group sizes `[10, 13, 13]`.
- **Chemle Learn (`N=78`):** Level-0 features `['heavy_atom_tercile', 'polarity_band']`, `k_min = 6`, group sizes `[6, 11, 13, 13, 16, 19]`.
- **Seqle Ranked (`N=28`):** Level-0 features `['charge_band']`, `k_min = 7`, group sizes `[7, 8, 13]`.
- **Seqle Learn (`N=30`):** Level-0 features `['charge_band', 'gravy_band']`, `k_min = 3`, group sizes `[3, 4, 4, 5, 5, 9]`.

### Learn Category Pack k-Anonymity Breakdown

Category packs are filtered subsets for exploratory Learn study; small packs naturally fall below the whole-catalog `k_min >= 3` threshold and are explicitly marked `exploratory` in metadata:

| Category Pack | Count | `k_min` | Group Sizes | Exploratory (`k_min < 3`) |
|---|---|---|---|---|
| Alkaloid / Metabolite | 5 | 1 | `[1, 1, 3]` | `true` |
| Biogenic Amine | 3 | 3 | `[3]` | `false` |
| Cofactor / Vitamin | 2 | 1 | `[1, 1]` | `true` |
| Hormone / Metabolite | 1 | 1 | `[1]` | `true` |
| Hormone / Neurotransmitter | 2 | 1 | `[1, 1]` | `true` |
| Hydrocarbon Degradation | 13 | 1 | `[1, 2, 10]` | `true` |
| Microbial Natural Product | 10 | 1 | `[1, 2, 3, 4]` | `true` |
| Neurotransmitter | 5 | 1 | `[1, 2, 2]` | `true` |
| Plant Natural Product | 11 | 1 | `[1, 2, 2, 3, 3]` | `true` |
| Primary Metabolite | 14 | 1 | `[1, 1, 3, 9]` | `true` |
| Therapeutic Drug | 12 | 1 | `[1, 1, 1, 4, 5]` | `true` |

### Solver Gate Verification (Clue-Conditioned Greedy Entropy)

- **Ranked Chemle (`targets=36`, `probes=36`, `cap=6`):** `mean=3.5`, `median=3.0`, `max=5`, `<=2=13.89%`, `<=6=100.0%`, `min_opening_expected_remaining=13.61`, `dist={'2': 5, '3': 15, '4': 9, '5': 7}`.
- **Learn Chemle (`targets=78`, `probes=78`, `cap=8`):** `mean=1.95`, `median=2.0`, `max=2`, `<=2=100.0%`, `<=6=100.0%`, `min_opening_expected_remaining=1.28`, `dist={'1': 4, '2': 74}`.
