# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "pandas",
#     "numpy",
# ]
# ///
"""Build-time RDKit Puzzle Compiler for BioWordle (Chemle + Seqle).

100% self-contained within the repository via `chem_engine.py` (no external
machine-specific skill paths). Validates chemical identities, stereochemistry,
chiral ECFP4-like Morgan fingerprints, progressive scaffold/carbon-skeleton
SVGs, peptide Henderson-Hasselbalch biophysics, and Shannon entropy information
gain IG(g) = H(T) - E_f[H(T|f)] for every puzzle.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from typing import Any, Dict, List

from rdkit import Chem, DataStructs

from chem_engine import (
    FINGERPRINT_SPEC,
    PEPTIDE_PKA_SPEC,
    STANDARDIZATION_SPEC,
    build_scaffold_clue_svg,
    compute_chemle_information_metrics,
    compute_chiral_morgan_fp,
    compute_descriptors_and_insights,
    compute_peptide_biochem,
    compute_seqle_information_metrics,
    render_responsive_svg,
    standardize_and_audit,
)

# ---------------------------------------------------------------------------
# 1. Functional Group SMARTS Catalog
# ---------------------------------------------------------------------------
FUNCTIONAL_GROUP_SMARTS: List[tuple[str, str]] = [
    ("Hydroxyl (-OH)", "[OX2H]"),
    ("Phenol (Ar-OH)", "[OX2H][c]"),
    ("Carboxylic Acid (-COOH)", "[CX3](=O)[OX2H1]"),
    ("Ketone / Aldehyde (C=O)", "[CX3]=[OX1]"),
    ("Ester / Lactone", "[#6][CX3](=O)[OX2H0][#6]"),
    ("Ether / Methoxy (C-O-C)", "[OD2]([#6])[#6]"),
    ("Primary / Secondary Amine", "[NX3;H2,H1;!$(NC=O)]"),
    ("Tertiary Amine", "[NX3;H0;!$(NC=O);!$(N=O)]"),
    ("Amide / Lactam", "[NX3][CX3](=[OX1])[#6]"),
    ("Indole / Pyrrole Ring", "[nH1]1cccc1"),
    ("Heteroaromatic N", "[n]"),
    ("Halogen (F/Cl/Br)", "[F,Cl,Br,I]"),
    ("Sulfur / Thioether", "[#16]"),
]

COMPILED_FG_SMARTS = [
    (name, Chem.MolFromSmarts(smarts))
    for name, smarts in FUNCTIONAL_GROUP_SMARTS
]


def detect_functional_groups(mol: Chem.Mol) -> List[str]:
    """Detects functional groups present in `mol` using verified SMARTS."""
    matched: List[str] = []
    for name, pattern in COMPILED_FG_SMARTS:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            matched.append(name)
    return matched


# ---------------------------------------------------------------------------
# 2. Curated CHEMLE Compound Library (36 Molecules)
# ---------------------------------------------------------------------------
CHEMLE_RAW_COMPOUNDS: List[Dict[str, str]] = [
    # --- Pack 1: Core Metabolites, Neurotransmitters & Cofactors ---
    {
        "id": "caffeine",
        "name": "Caffeine",
        "Category": "Alkaloid / Metabolite",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "hint_role": "Adenosine A1/A2A receptor antagonist & methylxanthine CNS stimulant.",
        "hint_origin": "Biosynthesized from xanthosine in Coffea arabica and Camellia sinensis.",
        "structural_note": "1,3,7-Trimethylxanthine bicyclic purine-2,6-dione core with zero H-bond donors.",
    },
    {
        "id": "dopamine",
        "name": "Dopamine",
        "Category": "Neurotransmitter",
        "smiles": "NCCc1ccc(O)c(O)c1",
        "hint_role": "Catecholamine neurotransmitter governing reward, motivation, and motor control.",
        "hint_origin": "Synthesized from L-DOPA via aromatic L-amino acid decarboxylase (DOPA decarboxylase).",
        "structural_note": "Shares the 1,2-dihydroxybenzene ring of Catechol with a 2-aminoethyl side chain at C4.",
    },
    {
        "id": "serotonin",
        "name": "Serotonin",
        "Category": "Neurotransmitter",
        "smiles": "NCCc1c[nH]c2ccc(O)cc12",
        "hint_role": "Indolamine monoamine (5-HT) regulating mood, sleep, and GI motility.",
        "hint_origin": "Derived from L-tryptophan via tryptophan hydroxylase and decarboxylation.",
        "structural_note": "5-Hydroxytryptamine: differs from L-Tryptophan by C5 hydroxylation and loss of the alpha-carboxyl group.",
    },
    {
        "id": "adrenaline",
        "name": "Adrenaline",
        "Category": "Hormone / Neurotransmitter",
        "smiles": "CNC[C@H](O)c1ccc(O)c(O)c1",
        "hint_role": "Fight-or-flight catecholamine hormone acting on alpha- and beta-adrenergic receptors.",
        "hint_origin": "Produced by adrenal medulla chromaffin cells via PNMT methylation of noradrenaline.",
        "structural_note": "(R)-Enantiomer catecholamine bearing both a benzylic beta-hydroxyl stereocenter and an N-methyl secondary amine.",
    },
    {
        "id": "melatonin",
        "name": "Melatonin",
        "Category": "Hormone / Metabolite",
        "smiles": "CC(=O)NCCc1c[nH]c2ccc(OC)cc12",
        "hint_role": "Pineal gland methoxy-indole hormone synchronizing circadian rhythms.",
        "hint_origin": "Formed by N-acetylation and O-methylation of serotonin in the dark phase.",
        "structural_note": "Converts Serotonin's polar 5-OH and primary amine into a lipophilic 5-methoxy ether and N-acetyl amide.",
    },
    {
        "id": "histamine",
        "name": "Histamine",
        "Category": "Biogenic Amine",
        "smiles": "NCCc1c[nH]cn1",
        "hint_role": "Imidazole biogenic amine mediating allergic mast-cell degranulation and gastric acid secretion.",
        "hint_origin": "Produced by histidine decarboxylase removal of the carboxyl group from L-histidine.",
        "structural_note": "Smallest heteroaromatic biogenic amine in the set (111.1 Da), pairing a 5-membered imidazole ring with an ethylamine tail.",
    },
    {
        "id": "tryptophan",
        "name": "L-Tryptophan",
        "Category": "Primary Metabolite",
        "smiles": "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",
        "hint_role": "Largest standard proteinogenic amino acid and precursor to auxin, NAD+, and serotonin.",
        "hint_origin": "Encoded by the single codon UGG; contains an electron-rich indole side chain.",
        "structural_note": "(S)-Alpha-amino acid biosynthetic precursor to Serotonin, Melatonin, and the microbial bis-indole pigment Violacein.",
    },
    {
        "id": "ascorbic_acid",
        "name": "Ascorbic Acid",
        "Category": "Cofactor / Vitamin",
        "smiles": "O=C1O[C@H]([C@@H](O)CO)C(O)=C1O",
        "hint_role": "Redox-active enediol lactone cofactor (Vitamin C) required for prolyl/lysyl hydroxylases.",
        "hint_origin": "Essential dietary micronutrient in primates due to loss of L-gulonolactone oxidase (GULO).",
        "structural_note": "Five-membered gamma-lactone ring with a vinylogous carboxylic acid (enediol) system and two chiral centers.",
    },
    {
        "id": "biotin",
        "name": "Biotin",
        "Category": "Cofactor / Vitamin",
        "smiles": "O=C(O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12",
        "hint_role": "Covalent CO2-carrier cofactor (Vitamin B7) for acetyl-CoA and pyruvate carboxylases.",
        "hint_origin": "Fused ureido-tetrahydrothiophene bicycle with ultra-high affinity for avidin/streptavidin.",
        "structural_note": "Bicyclic cis-fused cyclic urea + tetrahydrothiophene core with three contiguous stereocenters and a valeric acid tail.",
    },
    {
        "id": "citric_acid",
        "name": "Citric Acid",
        "Category": "Primary Metabolite",
        "smiles": "O=C(O)CC(O)(CC(=O)O)C(=O)O",
        "hint_role": "Achiral prochiral tricarboxylic acid at the entry point of the Krebs / TCA cycle.",
        "hint_origin": "Formed by citrate synthase condensation of oxaloacetate and acetyl-CoA.",
        "structural_note": "Only completely acyclic molecule in the puzzle set (Rings = 0); symmetric C6 tertiary alcohol with three carboxyl groups.",
    },
    {
        "id": "shikimic_acid",
        "name": "Shikimic Acid",
        "Category": "Primary Metabolite",
        "smiles": "O=C(O)C1=C[C@@H](O)[C@@H](O)[C@H](O)C1",
        "hint_role": "Central cyclohexene intermediate in aromatic amino acid & phenylpropanoid biosynthesis.",
        "hint_origin": "Isolated from Illicium verum (star anise); industrial starting material for oseltamivir.",
        "structural_note": "Non-aromatic cyclohexene-1-carboxylic acid with three contiguous stereocenters (3R,4S,5R-trihydroxy) before aromatization to chorismate.",
    },
    {
        "id": "cholesterol",
        "name": "Cholesterol",
        "Category": "Primary Metabolite",
        "smiles": "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(C)CC[C@H]3[C@H]2CC=C4[C@@]3(C)CC[C@@H](O)C4",
        "hint_role": "Tetracyclic cyclopentanoperhydrophenanthrene sterol modulating eukaryotic membrane fluidity.",
        "hint_origin": "Formed by oxidosqualene cyclase chair-boat-chair-boat cyclization of 2,3-oxidosqualene.",
        "structural_note": "8 chiral centers, 4 aliphatic rings, 0 aromatic rings, and extreme lipophilicity (cLogP 7.39).",
    },
    # --- Pack 2: Hydrocarbon Degradation & Microbial Metabolites (CANT-HYD / NMZ focus) ---
    {
        "id": "catechol",
        "name": "Catechol",
        "Category": "Hydrocarbon Degradation",
        "smiles": "Oc1ccccc1O",
        "hint_role": "Central 1,2-dihydroxybenzene ring-cleavage hub in aerobic aromatic hydrocarbon degradation.",
        "hint_origin": "Cleaved via ortho-fission (CatA, intradiol) or meta-fission (DmpB/XylE, extradiol) dioxygenases.",
        "structural_note": "Unsubstituted ortho-diphenol formed by decarboxylative hydroxylation of Salicylic Acid (NahG).",
    },
    {
        "id": "salicylic_acid",
        "name": "Salicylic Acid",
        "Category": "Hydrocarbon Degradation",
        "smiles": "O=C(O)c1ccccc1O",
        "hint_role": "Key intermediate in upper naphthalene degradation (NahG salicylate hydroxylase) & plant defense hormone.",
        "hint_origin": "Converted to catechol by FAD-dependent salicylate 1-monooxygenase in Pseudomonas NAH7 plasmids.",
        "structural_note": "2-Hydroxybenzoic acid: links upper Naphthalene catabolism to Catechol, and is O-acetylated synthetically to make Aspirin.",
    },
    {
        "id": "protocatechuic_acid",
        "name": "Protocatechuic Acid",
        "Category": "Hydrocarbon Degradation",
        "smiles": "O=C(O)c1ccc(O)c(O)c1",
        "hint_role": "3,4-Dihydroxybenzoate funneling hub of the bacterial beta-ketoadipate pathway (PcaHG).",
        "hint_origin": "Central catabolic node for lignin-derived aromatics, phthalate, and p-hydroxybenzoate.",
        "structural_note": "Positional isomer of Gentisic Acid (both C7H6O4, MW 154.1 Da): bears a catechol-like 3,4-ortho-diol rather than a 2,5-para-diol.",
    },
    {
        "id": "gentisic_acid",
        "name": "Gentisic Acid",
        "Category": "Hydrocarbon Degradation",
        "smiles": "O=C(O)c1cc(O)ccc1O",
        "hint_role": "2,5-Dihydroxybenzoate intermediate cleaved to maleylpyruvate by gentisate 1,2-dioxygenase.",
        "hint_origin": "Alternative aerobic catabolic route for naphthalene, salicylate, and 3-hydroxybenzoate.",
        "structural_note": "Constitutional isomer of Protocatechuic Acid (C7H6O4) with identical MW (154.1), TPSA (77.8), and cLogP (0.80), distinguished by 2,5-hydroquinone substitution.",
    },
    {
        "id": "naphthalene",
        "name": "Naphthalene",
        "Category": "Hydrocarbon Degradation",
        "smiles": "c1ccc2ccccc2c1",
        "hint_role": "Simplest fused bicyclic polycyclic aromatic hydrocarbon (PAH) model substrate.",
        "hint_origin": "Oxidized to cis-1,2-dihydroxy-1,2-dihydronaphthalene by Rieske non-heme iron dioxygenase (NahAcAd).",
        "structural_note": "Pure aromatic hydrocarbon (C10H8, TPSA = 0, Fsp3 = 0) with 2 fused benzene rings.",
    },
    {
        "id": "phenanthrene",
        "name": "Phenanthrene",
        "Category": "Hydrocarbon Degradation",
        "smiles": "c1ccc2c(c1)ccc3ccccc23",
        "hint_role": "Angular tricyclic PAH featuring the reactive K-region (C9=C10) and bay-region architecture.",
        "hint_origin": "Signature CANT-HYD / marine oil-spill biodegradation benchmark substrate.",
        "structural_note": "Angular 3-ring homologue of Naphthalene (C14H10 vs. C10H8), increasing cLogP from 2.84 to 4.00.",
    },
    {
        "id": "camphor",
        "name": "Camphor",
        "Category": "Hydrocarbon Degradation",
        "smiles": "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2",
        "hint_role": "Bicyclic monoterpene ketone and canonical substrate of cytochrome P450cam (CYP101A1).",
        "hint_origin": "Regio- and stereoselectively hydroxylated at C5-exo by Pseudomonas putida P450cam.",
        "structural_note": "Bridged bicyclo[2.2.1]heptan-2-one (bornanone) monoterpene with Fsp3 = 0.90 and zero aromatic rings.",
    },
    {
        "id": "pyocyanin",
        "name": "Pyocyanin",
        "Category": "Microbial Natural Product",
        "smiles": "Cn1c2ccccc2nc2ccc(=O)cc21",
        "hint_role": "Blue redox-active tricyclic phenazine zwitterion & quorum-sensing virulence factor.",
        "hint_origin": "Secreted by Pseudomonas aeruginosa; generates reactive oxygen species via NAD(P)H oxidation.",
        "structural_note": "Planar N-methylphenazin-1-one tricyclic heteroaromatic pigment (Fsp3 = 0.08, 3 aromatic rings).",
    },
    {
        "id": "prodigiosin",
        "name": "Prodigiosin",
        "Category": "Microbial Natural Product",
        "smiles": "CCCCC1=CC(=C(N1)C=C2C=CC(=N2)C3=CC(=CN3)OC)C",
        "hint_role": "Blood-red tripyrrole bacterial secondary metabolite with H+/Cl- symporter & antimicrobial activity.",
        "hint_origin": "Biosynthesized via hybrid NRPS-PKS pig/red gene cluster in Serratia marcescens and marine Vibrio.",
        "structural_note": "Contains three conjugated pyrrole rings (methoxybipyrrole + alkylpyrrole) linked by a methene bridge.",
    },
    {
        "id": "violacein",
        "name": "Violacein",
        "Category": "Microbial Natural Product",
        "smiles": "O=C1NC(=C(c2c[nH]c3ccccc23)C1=C1Nc2ccccc2C1=O)c1c[nH]c2ccc(O)cc12",
        "hint_role": "Deep-violet bis-indole pigment assembled by the 5-enzyme VioABCDE pathway from two L-tryptophans.",
        "hint_origin": "Signature metabolite of Chromobacterium violaceum and marine Pseudoalteromonas.",
        "structural_note": "Condenses two L-Tryptophan molecules into a 5-ring system (4 aromatic + 1 aliphatic pyrrolidone, Fsp3 = 0.00).",
    },
    {
        "id": "geosmin",
        "name": "Geosmin",
        "Category": "Microbial Natural Product",
        "smiles": "C[C@H]1CCCC2(C)[C@H]1CC[C@@]2(C)O",
        "hint_role": "Trans-decalin sesquiterpene alcohol responsible for the scent of wet earth (petrichor).",
        "hint_origin": "Formed from farnesyl diphosphate by bifunctional geosmin synthase (Sco6073) in Streptomyces.",
        "structural_note": "Fully saturated bicyclic tertiary alcohol (Fsp3 = 1.00, 2 aliphatic rings, 3 stereocenters).",
    },
    # --- Pack 3: Plant Natural Products & Secondary Metabolites ---
    {
        "id": "vanillin",
        "name": "Vanillin",
        "Category": "Plant Natural Product",
        "smiles": "COc1cc(C=O)ccc1O",
        "hint_role": "Phenolic benzaldehyde derivative and primary flavor molecule of cured vanilla beans.",
        "hint_origin": "Also formed during microbial ferulic acid and lignin depolymerization.",
        "structural_note": "4-Hydroxy-3-methoxybenzaldehyde: precedes Protocatechuic Acid in microbial lignin/vanillate O-demethylation.",
    },
    {
        "id": "resveratrol",
        "name": "Resveratrol",
        "Category": "Plant Natural Product",
        "smiles": "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1",
        "hint_role": "Trans-stilbene phytoalexin polyphenol synthesised by stilbene synthase in grapevines.",
        "hint_origin": "Combines one p-coumaroyl-CoA and three malonyl-CoA extender units via type III PKS.",
        "structural_note": "(E)-Stilbene triol with two phenolic rings linked by a trans-alkene bridge (Fsp3 = 0.00, HBD = 3).",
    },
    {
        "id": "curcumin",
        "name": "Curcumin",
        "Category": "Plant Natural Product",
        "smiles": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
        "hint_role": "Symmetric diarylheptanoid beta-diketone pigment giving turmeric its bright yellow fluorescence.",
        "hint_origin": "Classic pan-assay interference (PAINS) / Michael acceptor benchmark in medicinal chemistry.",
        "structural_note": "Two vanillin-like feruloyl halves joined at a central methylene active beta-diketone.",
    },
    {
        "id": "capsaicin",
        "name": "Capsaicin",
        "Category": "Plant Natural Product",
        "smiles": "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O",
        "hint_role": "Vanilloid amide agonist of the heat-activated TRPV1 ion channel nociceptor.",
        "hint_origin": "Synthesized in Capsicum placenta from vanillylamine and a branched-chain fatty acyl-CoA.",
        "structural_note": "Shares the 4-hydroxy-3-methoxybenzyl headgroup of Vanillin, coupled via an amide bond to a 9-carbon hydrophobic tail (RotB = 9).",
    },
    {
        "id": "nicotine",
        "name": "Nicotine",
        "Category": "Alkaloid / Metabolite",
        "smiles": "CN1CCC[C@H]1c2cccnc2",
        "hint_role": "Pyridine-pyrrolidine alkaloid agonist of nicotinic acetylcholine receptors (nAChRs).",
        "hint_origin": "Biosynthesized in Nicotiana roots by coupling ornithine-derived N-methylpyrrolinium with nicotinic acid.",
        "structural_note": "(S)-Enantiomer biaryl-like system linking 1 aromatic pyridine ring with 1 aliphatic N-methylpyrrolidine ring.",
    },
    {
        "id": "quinine",
        "name": "Quinine",
        "Category": "Alkaloid / Metabolite",
        "smiles": "COc1ccc2nccc([C@@H](O)[C@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1",
        "hint_role": "Cinchona quinoline-quinuclidine alkaloid that inhibits Plasmodium hemozoin biocrystallization.",
        "hint_origin": "Derived from strictosidine via tryptophan-terpenoid indole-to-quinoline rearrangement.",
        "structural_note": "Contains 4 rings (2 aromatic quinoline + 2 bridged aliphatic quinuclidine) and 4 stereocenters.",
    },
    {
        "id": "artemisinin",
        "name": "Artemisinin",
        "Category": "Plant Natural Product",
        "smiles": "C[C@@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@]23OO4",
        "hint_role": "Sesquiterpene lactone containing an essential endoperoxide bridge cleaved by heme iron(II).",
        "hint_origin": "Isolated from Artemisia annua by Tu Youyou (2015 Nobel Prize in Physiology or Medicine).",
        "structural_note": "Tetracyclic cage sesquiterpene with 7 stereocenters, Fsp3 = 0.93, and zero aromatic rings.",
    },
    # --- Pack 4: Classic Therapeutics & Bioactive Probes ---
    {
        "id": "aspirin",
        "name": "Aspirin",
        "Category": "Therapeutic Drug",
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "hint_role": "Irreversible covalent inhibitor acetylating Ser530/Ser516 of COX-1 and COX-2.",
        "hint_origin": "O-Acetyl derivative of salicylic acid introduced clinically in 1899.",
        "structural_note": "Acetylsalicylic acid: masks the phenolic -OH of Salicylic Acid as a labile phenyl acetate ester.",
    },
    {
        "id": "ibuprofen",
        "name": "Ibuprofen",
        "Category": "Therapeutic Drug",
        "smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
        "hint_role": "2-Arylpropionic acid reversible non-selective cyclooxygenase (COX) NSAID inhibitor.",
        "hint_origin": "Features an isobutyl-substituted phenyl ring and an alpha-methyl propanoate group.",
        "structural_note": "Lipophilic NSAID (cLogP 3.07) with a single benzene ring, carboxylic acid, and zero nitrogen atoms.",
    },
    {
        "id": "acetaminophen",
        "name": "Acetaminophen",
        "Category": "Therapeutic Drug",
        "smiles": "CC(=O)Nc1ccc(O)cc1",
        "hint_role": "Para-aminophenol analgesic/antipyretic (Paracetamol); oxidized by CYP2E1 to reactive NAPQI.",
        "hint_origin": "Minimalist 151 Da acylated phenol.",
        "structural_note": "N-(4-Hydroxyphenyl)acetamide: pairs a para-phenol with an acetamide group.",
    },
    {
        "id": "penicillin_g",
        "name": "Penicillin G",
        "Category": "Therapeutic Drug",
        "smiles": "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
        "hint_role": "Fused beta-lactam / thiazolidine antibiotic mimicking D-Ala-D-Ala to acylate bacterial PBPs.",
        "hint_origin": "Assembled by ACV synthetase (NRPS) and isopenicillin N synthase in Penicillium chrysogenum.",
        "structural_note": "Contains 1 aromatic phenylacetyl ring + 2 fused aliphatic rings (strained 4-membered beta-lactam + 5-membered thiazolidine).",
    },
    {
        "id": "ciprofloxacin",
        "name": "Ciprofloxacin",
        "Category": "Therapeutic Drug",
        "smiles": "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
        "hint_role": "Second-generation fluoroquinolone trapping bacterial DNA gyrase and topoisomerase IV.",
        "hint_origin": "Features a cyclopropyl N1 substituent, C6 fluorine, and C7 piperazine ring.",
        "structural_note": "4-Ring system (2 aromatic quinolone + cyclopropyl + piperazine) containing a single fluorine atom.",
    },
    {
        "id": "fluoxetine",
        "name": "Fluoxetine",
        "Category": "Therapeutic Drug",
        "smiles": "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",
        "hint_role": "Archetypal selective serotonin reuptake inhibitor (SSRI) blocking the SERT transporter.",
        "hint_origin": "Para-trifluoromethylphenoxy phenylpropylamine with high blood-brain barrier permeability.",
        "structural_note": "Contains two benzene rings, a secondary methylamine, an aryl ether, and a -CF3 group.",
    },
]


# ---------------------------------------------------------------------------
# 3. Curated SEQLE 5-Residue Motif Library with Primary Class & Secondary Traits
# ---------------------------------------------------------------------------
AA_PROPERTIES: Dict[str, Dict[str, Any]] = {
    "A": {"name": "Alanine", "three": "Ala", "group": "Aliphatic", "traits": ["Aliphatic", "Nonpolar", "Small"], "hydropathy": 1.8, "charge": 0.0, "mw": 89.09},
    "V": {"name": "Valine", "three": "Val", "group": "Aliphatic", "traits": ["Aliphatic", "Beta-branched", "Hydrophobic"], "hydropathy": 4.2, "charge": 0.0, "mw": 117.15},
    "I": {"name": "Isoleucine", "three": "Ile", "group": "Aliphatic", "traits": ["Aliphatic", "Beta-branched", "Hydrophobic"], "hydropathy": 4.5, "charge": 0.0, "mw": 131.17},
    "L": {"name": "Leucine", "three": "Leu", "group": "Aliphatic", "traits": ["Aliphatic", "Gamma-branched", "Hydrophobic"], "hydropathy": 3.8, "charge": 0.0, "mw": 131.17},
    "M": {"name": "Methionine", "three": "Met", "group": "Aliphatic", "traits": ["Aliphatic", "Thioether Sulfur", "Hydrophobic"], "hydropathy": 1.9, "charge": 0.0, "mw": 149.21},
    "F": {"name": "Phenylalanine", "three": "Phe", "group": "Aromatic", "traits": ["Aromatic", "Hydrophobic"], "hydropathy": 2.8, "charge": 0.0, "mw": 165.19},
    "Y": {"name": "Tyrosine", "three": "Tyr", "group": "Aromatic", "traits": ["Aromatic", "Polar (Phenolic -OH)", "Phosphorylation Site"], "hydropathy": -1.3, "charge": 0.0, "mw": 181.19},
    "W": {"name": "Tryptophan", "three": "Trp", "group": "Aromatic", "traits": ["Aromatic (Indole)", "Polar (Indole N-H)", "Bulky"], "hydropathy": -0.9, "charge": 0.0, "mw": 204.23},
    "S": {"name": "Serine", "three": "Ser", "group": "Polar", "traits": ["Polar (-OH)", "Nucleophilic", "Small"], "hydropathy": -0.8, "charge": 0.0, "mw": 105.09},
    "T": {"name": "Threonine", "three": "Thr", "group": "Polar", "traits": ["Polar (-OH)", "Beta-branched"], "hydropathy": -0.7, "charge": 0.0, "mw": 119.12},
    "N": {"name": "Asparagine", "three": "Asn", "group": "Polar", "traits": ["Polar (Carboxamide)", "H-Bond Donor/Acceptor"], "hydropathy": -3.5, "charge": 0.0, "mw": 132.12},
    "Q": {"name": "Glutamine", "three": "Gln", "group": "Polar", "traits": ["Polar (Carboxamide)", "H-Bond Donor/Acceptor"], "hydropathy": -3.5, "charge": 0.0, "mw": 146.15},
    "K": {"name": "Lysine", "three": "Lys", "group": "Positive", "traits": ["Positive (Primary Amine, pKa 10.5)", "Polar"], "hydropathy": -3.9, "charge": 1.0, "mw": 146.19},
    "R": {"name": "Arginine", "three": "Arg", "group": "Positive", "traits": ["Positive (Guanidinium, pKa 12.5)", "Bidentate H-Bonding"], "hydropathy": -4.5, "charge": 1.0, "mw": 174.20},
    "H": {"name": "Histidine", "three": "His", "group": "Positive", "traits": ["Positive (Imidazole, pKa 6.0)", "Aromatic (Heteroaryl)", "Metal Ligand"], "hydropathy": -3.2, "charge": 0.1, "mw": 155.16},
    "D": {"name": "Aspartate", "three": "Asp", "group": "Negative", "traits": ["Negative (Carboxylate, pKa 3.9)", "Polar"], "hydropathy": -3.5, "charge": -1.0, "mw": 133.10},
    "E": {"name": "Glutamate", "three": "Glu", "group": "Negative", "traits": ["Negative (Carboxylate, pKa 4.1)", "Polar"], "hydropathy": -3.5, "charge": -1.0, "mw": 147.13},
    "C": {"name": "Cysteine", "three": "Cys", "group": "Special", "traits": ["Special (Thiol -SH, pKa 8.3)", "Polar / Redox Disulfide", "Metal Ligand"], "hydropathy": 2.5, "charge": 0.0, "mw": 121.16},
    "G": {"name": "Glycine", "three": "Gly", "group": "Special", "traits": ["Special (Achiral, R=H)", "Conformational Flexibility"], "hydropathy": -0.4, "charge": 0.0, "mw": 75.07},
    "P": {"name": "Proline", "three": "Pro", "group": "Special", "traits": ["Special (Cyclic Secondary Amide)", "Helix Breaker / Turn Inducer"], "hydropathy": -1.6, "charge": 0.0, "mw": 115.13},
}

SEQLE_RAW_MOTIFS: List[Dict[str, str]] = [
    {
        "seq": "GHSMG",
        "title": "Serine Hydrolase / Cutinase GXSXG Elbow",
        "pattern": "G-X-S-X-G",
        "family": "Catalytic Active Site",
        "hint_role": "Nucleophilic elbow of alpha/beta-hydrolases, lipases, PETases, and esterases.",
        "hint_detail": "Flanking glycines prevent steric clashes at the sharp beta-turn carrying the catalytic Ser.",
        "context_requirement": "Embedded at a tight beta-strand-to-alpha-helix nucleophilic elbow coupled to a His-Asp catalytic dyad.",
        "exemplar_source": "Ideonella sakaiensis PETase / Fusarium solani Cutinase (PDB: 5XJH, 1CUS)",
    },
    {
        "seq": "KDELR",
        "title": "ER Luminal Retention & Retrieval Signal",
        "pattern": "K-D-E-L-[R/K]",
        "family": "Subcellular Targeting",
        "hint_role": "Recognized in the cis-Golgi by KDEL receptors to retrograde-transport escaped ER chaperones.",
        "hint_detail": "Combines a basic N-terminal Lys with an acidic Asp-Glu pair and hydrophobic Leu.",
        "context_requirement": "Must reside at the extreme C-terminus of soluble ER luminal proteins (BiP/GRP78, PDI, calreticulin).",
        "exemplar_source": "Human BiP / HSPA5 & KDELR2 Complex (PDB: 6I6H)",
    },
    {
        "seq": "DEVDG",
        "title": "Caspase-3 Apoptotic Protease Cleavage Site",
        "pattern": "D-E-V-D-G",
        "family": "Protease Cleavage Motif",
        "hint_role": "Canonical P4-P1' recognition sequence cleaved after the second Asp in PARP-1 during apoptosis.",
        "hint_detail": "Strictly requires Asp at both P4 and P1 positions followed by a small P1' Gly.",
        "context_requirement": "Requires an exposed solvent-accessible loop between the N-terminal DNA-binding and C-terminal catalytic domains.",
        "exemplar_source": "Human Poly [ADP-ribose] polymerase 1 (PARP1, UniProt: P09874, residues 211-215)",
    },
    {
        "seq": "RGDVF",
        "title": "Fibronectin Integrin-Binding Motif",
        "pattern": "R-G-D-V-F",
        "family": "Cell Adhesion",
        "hint_role": "Extracellular matrix tripeptide core recognized by alphaV-beta3 and alpha5-beta1 integrins.",
        "hint_detail": "Starts with a cationic Arg, flexible Gly spacer, and anionic Asp; basis of cilengitide.",
        "context_requirement": "Presented at the apex of a flexible type III fibronectin beta-hairpin loop; Asp coordinates the MIDAS Mg2+ ion.",
        "exemplar_source": "Integrin alphaV-beta3 in complex with cyclic RGDfV peptide (PDB: 1L5G)",
    },
    {
        "seq": "LPETG",
        "title": "Sortase A Cell-Wall Anchoring Motif",
        "pattern": "L-P-X-T-G",
        "family": "Transpeptidation Tag",
        "hint_role": "Cleaved between Thr and Gly by Staphylococcus aureus Sortase A and ligated to pentaglycine.",
        "hint_detail": "Widely used for site-specific enzymatic protein bioconjugation ('sortagging').",
        "context_requirement": "Located upstream of a hydrophobic transmembrane segment and positively charged C-terminal tail in Gram-positive surface proteins.",
        "exemplar_source": "Staphylococcus aureus Sortase A substrate complex (PDB: 1T2W)",
    },
    {
        "seq": "HEIGH",
        "title": "Zincin Metalloprotease Catalytic Core",
        "pattern": "H-E-X-X-H",
        "family": "Metal-Binding Catalytic Site",
        "hint_role": "Two histidines coordinate the catalytic Zn2+ ion while glutamate acts as general base.",
        "hint_detail": "Found in thermolysin, matrix metalloproteinases (MMPs), and botulinum neurotoxin.",
        "context_requirement": "Must fold as an active-site alpha-helix so i and i+4 His imidazole rings project onto the same face toward Zn2+.",
        "exemplar_source": "Human Matrix Metalloproteinase-1 / Collagenase (MMP1, UniProt: P03956, PDB: 1HFC)",
    },
    {
        "seq": "CCHHC",
        "title": "Retrotransposon / Gag Zinc Knuckle Core",
        "pattern": "C-X2-C-X4-H-X4-C",
        "family": "Nucleic Acid Binding",
        "hint_role": "Tetrahedral Zn2+-coordinating Cys/His cluster that binds single-stranded viral RNA.",
        "hint_detail": "Compressed 5-residue coordination signature representing the Cys/His ligand sphere.",
        "context_requirement": "Requires deprotonated Cys thiolates and neutral His imidazole spaced by short beta-hairpin knuckles around Zn2+.",
        "exemplar_source": "Retroviral Nucleocapsid CCHC Zinc Knuckle (PDB: 1MFS)",
    },
    {
        "seq": "YIGSR",
        "title": "Laminin B1 Chain Receptor-Binding Pentapeptide",
        "pattern": "Y-I-G-S-R",
        "family": "Cell Adhesion",
        "hint_role": "Basement membrane laminin-1 motif promoting epithelial cell attachment and inhibiting metastasis.",
        "hint_detail": "Begins with an aromatic Tyr and hydrophobic Ile, ending with a basic Arg.",
        "context_requirement": "Located within cysteine-rich EGF-like repeats of the laminin beta-1 short arm.",
        "exemplar_source": "Human Laminin Subunit Beta-1 (LAMB1, UniProt: P07942, residues 929-933)",
    },
    {
        "seq": "NPVYK",
        "title": "Clathrin Endocytosis NPXY Internalization Signal",
        "pattern": "N-P-X-Y-[K/R]",
        "family": "Membrane Trafficking",
        "hint_role": "Cytoplasmic tail motif of LDL, insulin, and integrin receptors recruiting PTB-domain adaptors.",
        "hint_detail": "Asn-Pro forms a tight beta-turn positioning the aromatic Tyr into the adaptor pocket.",
        "context_requirement": "Located in the unstructured cytosolic juxtamembrane tail of Type I transmembrane receptors.",
        "exemplar_source": "Human LDL Receptor & Autosomal Recessive Hypercholesterolemia (ARH) PTB Complex (PDB: 3SO6)",
    },
    {
        "seq": "CVIMK",
        "title": "K-Ras4B C-Terminal CAAX Prenylation Box",
        "pattern": "C-A-A-X",
        "family": "Post-Translational Lipidation",
        "hint_role": "Farnesyltransferase attaches a C15 isoprenoid lipid to the conserved Cys thiol.",
        "hint_detail": "Terminal Met directs high-affinity recognition by farnesyltransferase over GGTase-I.",
        "context_requirement": "Strictly requires C-terminal positioning preceded by a polybasic Lys stretch; -VIM is cleaved off by Rce1 after farnesylation.",
        "exemplar_source": "Human GTPase KRas Isoform 2B (KRAS, UniProt: P01116, PDB: 1D8D)",
    },
    {
        "seq": "DEADR",
        "title": "DEAD-Box RNA Helicase Walker B Motif",
        "pattern": "D-E-A-D-[R/K]",
        "family": "NTPase / Helicase Core",
        "hint_role": "Coordinates Mg2+ and polarizes the attacking water molecule for ATP hydrolysis in RNA unwinding.",
        "hint_detail": "Iconic sequence motif that gave the entire DEAD-box helicase superfamily its name.",
        "context_requirement": "Positioned at the C-terminal tip of RecA-like domain 1 strand beta-3 adjacent to the Walker A P-loop.",
        "exemplar_source": "Human ATP-Dependent RNA Helicase DDX4 / eIF4A1 (UniProt: P60842, PDB: 2ZU6)",
    },
    {
        "seq": "GAGTG",
        "title": "Walker A / P-Loop Nucleotide-Binding Core",
        "pattern": "G-X4-G-K-[T/S]",
        "family": "NTPase / Helicase Core",
        "hint_role": "Glycine-rich phosphate-binding loop wrapping around beta/gamma phosphates of ATP or GTP.",
        "hint_detail": "Alternating glycines allow backbone amides to form a tight anion-binding nest.",
        "context_requirement": "Connects strand beta-1 and helix alpha-1 in Rossmann-fold NTPases (Ras, myosin, ABC transporters, kinases).",
        "exemplar_source": "F1-ATPase / RecA P-Loop NTPase Core (PROSITE: PS00017)",
    },
    {
        "seq": "VDFGL",
        "title": "Protein Kinase Activation Loop DFG Motif",
        "pattern": "[V/I]-D-F-G-[L/M]",
        "family": "Kinase Switch",
        "hint_role": "Conformational switch (DFG-in vs. DFG-out) governing Mg2+-ATP coordination in all eukaryotic kinases.",
        "hint_detail": "Target of Type II kinase inhibitors like imatinib; contains Asp-Phe-Gly at positions 2-4.",
        "context_requirement": "Marks the N-terminal anchor of the kinase activation segment; Asp chelates Mg2+ while Phe packs into the regulatory R-spine.",
        "exemplar_source": "Human Proto-Oncogene Tyrosine-Protein Kinase ABL1 / B-Raf (UniProt: P00519, PDB: 1IEP)",
    },
    {
        "seq": "HRDLA",
        "title": "Protein Kinase Catalytic Loop HRD Motif",
        "pattern": "H-R-D-L-A",
        "family": "Catalytic Active Site",
        "hint_role": "Catalytic Asp acts as proton acceptor for the substrate hydroxyl group during phosphoryl transfer.",
        "hint_detail": "Arg neutralizes the phosphorylated activation loop in RD-class tyrosine and Ser/Thr kinases.",
        "context_requirement": "Located in the C-lobe catalytic loop between beta-6 and beta-7 strands directly facing the ATP gamma-phosphate.",
        "exemplar_source": "Human Src / EGFR / ABL1 Tyrosine Kinase Catalytic Loop (UniProt: P12931, PDB: 2SRC)",
    },
    {
        "seq": "WSEWS",
        "title": "Type I Cytokine Receptor WSXWS Box",
        "pattern": "W-S-X-W-S",
        "family": "Structural Folding Motif",
        "hint_role": "Forms a Tryptophan-Arginine zipper ('Trp-stack') essential for folding EPO, IL-2, and GH receptors.",
        "hint_detail": "Symmetric pair of Trp-Ser dipeptides separated by a central acidic Glu.",
        "context_requirement": "Located in the membrane-proximal fibronectin type III domain of Class I hematopoietic cytokine receptors; often C-mannosylated.",
        "exemplar_source": "Human Erythropoietin Receptor (EPOR, UniProt: P19235, PDB: 1EER)",
    },
    {
        "seq": "PPPVP",
        "title": "Polyproline Type II Helix SH3 Ligand",
        "pattern": "P-X-X-P",
        "family": "Protein-Protein Interaction",
        "hint_role": "Left-handed PPII helix docking into aromatic grooves of Src Homology 3 (SH3) signaling domains.",
        "hint_detail": "Four rigid cyclic Pro residues flanking a single branched aliphatic Val.",
        "context_requirement": "Adopts an extended 3.0-residues-per-turn left-handed PPII helix recognized by conserved Trp/Tyr grooves on SH3 domains.",
        "exemplar_source": "SOS1 / GRB2 SH3 Proline-Rich Peptide Complex (PDB: 1GBR)",
    },
    {
        "seq": "INFEK",
        "title": "Ovalbumin SIINFEKL MHC Class I Epitope Core",
        "pattern": "S-I-I-N-F-E-K-L",
        "family": "Immunology Epitope",
        "hint_role": "Central TCR-contact and H-2Kb anchor core of the benchmark CD8+ T-cell octapeptide SIINFEKL.",
        "hint_detail": "Aromatic Phe5 serves as the deep central MHC pocket anchor while Asn4 and Glu6/Lys7 contact the TCR.",
        "context_requirement": "Processed by the immunoproteasome and TAP transporter to bind the groove of murine MHC Class I H-2Kb.",
        "exemplar_source": "Gallus gallus Ovalbumin 257-264 bound to MHC H-2Kb (PDB: 1VAC)",
    },
    {
        "seq": "YQRLN",
        "title": "Trans-Golgi Network Tyrosine Sorting Motif",
        "pattern": "Y-X-X-[L/I/V]",
        "family": "Membrane Trafficking",
        "hint_role": "Binds the mu-2 subunit of AP-2 clathrin adaptor complexes for rapid endocytic sorting.",
        "hint_detail": "Starts with an essential unphosphorylated Tyr and places hydrophobic Leu at Y+3.",
        "context_requirement": "Must be in an extended conformation within a cytosolic receptor tail; phosphorylation of Tyr blocks AP-2 binding.",
        "exemplar_source": "TGN38 Sorting Peptide bound to AP-2 Mu2 Subunit (PDB: 1BXX)",
    },
    {
        "seq": "KKMPQ",
        "title": "COPI Coatomer Dilysine ER Retrieval Signal",
        "pattern": "K-K-X-X",
        "family": "Subcellular Targeting",
        "hint_role": "C-terminal cytosolic tail motif of ER membrane proteins recognized by alpha/beta'-COP.",
        "hint_detail": "Starts with a twin positively charged Lys-Lys dyad.",
        "context_requirement": "Requires positioning at the -3/-4 or -4/-5 positions from the cytosolic C-terminus of Type I ER membrane proteins.",
        "exemplar_source": "WD40 Beta'-COP Domain bound to Dilysine Retrieval Peptide (PDB: 4J77)",
    },
    {
        "seq": "CAACH",
        "title": "Cytochrome c CxxCH Heme-Binding Motif",
        "pattern": "C-X-X-C-H",
        "family": "Cofactor Covalent Attachment",
        "hint_role": "Two cysteines form thioether bonds to heme vinyl groups while His serves as axial Fe ligand.",
        "hint_detail": "Signature motif of c-type cytochromes in respiratory and marine electron transport chains.",
        "context_requirement": "Matured by cytochrome c heme lyase / Ccm system in the periplasm or mitochondrial intermembrane space.",
        "exemplar_source": "Paracoccus / Shewanella Multiheme c-Type Cytochromes (PROSITE: PS51007)",
    },
    {
        "seq": "NLTGY",
        "title": "Eukaryotic N-Linked Glycosylation Sequon",
        "pattern": "N-X-[S/T]",
        "family": "Post-Translational Modification",
        "hint_role": "Oligosaccharyltransferase transfers Glc3Man9GlcNAc2 onto the Asn amide nitrogen in the ER.",
        "hint_detail": "Requires Asn at position 1, any non-Pro at position 2, and Thr/Ser at position 3.",
        "context_requirement": "Must translocate into the ER lumen and remain unfolded/accessible for STT3 catalytic subunit recognition.",
        "exemplar_source": "Oligosaccharyltransferase (OST) Acceptor Peptide Complex (PDB: 3WAK)",
    },
    {
        "seq": "DSGLS",
        "title": "Beta-TrCP SCF Ubiquitin Ligase Phosphodegron",
        "pattern": "D-S-G-X-X-S",
        "family": "Proteasomal Degradation",
        "hint_role": "Dual phosphorylation of the Ser residues triggers ubiquitination of I-kappa-B and beta-catenin.",
        "hint_detail": "Contains an acidic Asp followed by Ser-Gly-Leu-Ser.",
        "context_requirement": "Requires phosphorylation of both Ser residues by IKK or GSK3-beta to dock into the WD40 propeller of Beta-TrCP.",
        "exemplar_source": "Human Beta-Catenin / I-kappa-B-alpha bound to SCF-Beta-TrCP1 (PDB: 1P22)",
    },
    {
        "seq": "KENKR",
        "title": "APC/C Cdh1 Mitotic Exit KEN-Box Degron",
        "pattern": "K-E-N-X-R",
        "family": "Proteasomal Degradation",
        "hint_role": "Targets cell-cycle kinases (Aurora B, Nek2, Cdc6) for destruction during late mitosis and G1.",
        "hint_detail": "Named directly after its N-terminal Lys-Glu-Asn tripeptide.",
        "context_requirement": "Folds into an underwound conformation on the top surface of the Cdh1 WD40 beta-propeller.",
        "exemplar_source": "Anaphase-Promoting Complex Coactivator Cdh1-KEN Box Structure (PDB: 4BH6)",
    },
    {
        "seq": "PKKKR",
        "title": "SV40 Large T-Antigen Nuclear Localization Signal",
        "pattern": "P-K-K-K-R-K-V",
        "family": "Subcellular Targeting",
        "hint_role": "Prototype monopartite classical NLS recognized by Importin-alpha for nuclear pore import.",
        "hint_detail": "Initiated by a helix-breaking Pro followed by four consecutive cationic residues.",
        "context_requirement": "Binds in an extended conformation across the major ARM-repeat groove of Importin-alpha.",
        "exemplar_source": "SV40 Large T-Antigen NLS bound to Karyopherin / Importin-Alpha (PDB: 1BK6)",
    },
    {
        "seq": "FGGFD",
        "title": "Nociceptin / Orphanin FQ Neuropeptide Core",
        "pattern": "F-G-G-F-[T/D]",
        "family": "Neuropeptide GPCR Message",
        "hint_role": "N-terminal 'message' domain activating the NOP opioid-related GPCR (vs. YGGFM in enkephalin).",
        "hint_detail": "Two Phe rings separated by a flexible Gly-Gly spacer.",
        "context_requirement": "Requires a free N-terminal alpha-amino group on Phe1 to form a salt bridge with conserved Asp3.32 in the GPCR orthosteric pocket.",
        "exemplar_source": "Human Nociceptin / Orphanin FQ Receptor Complex (OPRL1, UniProt: P41146, PDB: 4EA3)",
    },
    {
        "seq": "YGGFM",
        "title": "Met-Enkephalin Endogenous Opioid Pentapeptide",
        "pattern": "Y-G-G-F-M",
        "family": "Neuropeptide GPCR Message",
        "hint_role": "Endogenous delta/mu-opioid receptor agonist pentapeptide cleaved from proenkephalin.",
        "hint_detail": "N-terminal Tyr1 phenol and amine mimic the tyramine moiety of morphine.",
        "context_requirement": "Functions as a free endogenous pentapeptide neurotransmitter released by proteolytic cleavage of Proenkephalin-A.",
        "exemplar_source": "Human Proenkephalin-A / Met-Enkephalin (PENK, UniProt: P01210, PDB: 1PLW)",
    },
    {
        "seq": "SYGYS",
        "title": "FUS / FET Low-Complexity Phase-Separation Motif",
        "pattern": "[S/G]-Y-[G/S]-Y-S",
        "family": "Biomolecular Condensate",
        "hint_role": "Pi-pi and cation-pi sticker-and-spacer repeat driving liquid-liquid phase separation of stress granules.",
        "hint_detail": "Palindrome of polar Ser, flexible Gly, and aromatic Tyr stickers.",
        "context_requirement": "Repeated in tandem (~27 times) within the intrinsically disordered N-terminal prion-like domain of FUS/EWSR1/TAF15.",
        "exemplar_source": "Human RNA-Binding Protein FUS Low-Complexity Fibril Core (UniProt: P35637, PDB: 5W3N)",
    },
    {
        "seq": "LDLAS",
        "title": "CRM1 / Exportin-1 Leucine-Rich Nuclear Export Core",
        "pattern": "L-X(2,3)-L-X(2,3)-L-X-L",
        "family": "Subcellular Targeting",
        "hint_role": "Hydrophobic ridge docking into the NES cleft of Exportin-1 (XPO1) in the presence of RanGTP.",
        "hint_detail": "Two Leu residues separated by an acidic Asp.",
        "context_requirement": "Forms an amphipathic alpha-helix-plus-loop that inserts hydrophobic side chains into five pockets of the CRM1 HEAT-repeat groove.",
        "exemplar_source": "Protein Kinase Inhibitor (PKI) Nuclear Export Signal bound to CRM1-RanGTP (PDB: 3NBY)",
    },
    {
        "seq": "NPLYR",
        "title": "Integrin Beta-Tail Talin-Binding PTB Motif",
        "pattern": "N-P-L-Y-R",
        "family": "Cell Adhesion / Inside-Out Signaling",
        "hint_role": "Binds the FERM domain of talin to switch integrins from bent-inactive to extended-active states.",
        "hint_detail": "Classic Asn-Pro-Leu-Tyr beta-turn capped by a C-terminal Arg.",
        "context_requirement": "Located in the membrane-proximal cytoplasmic tail of Integrin Beta-3/Beta-1; Tyr phosphorylation switches binding from Talin to Dok1.",
        "exemplar_source": "Human Integrin Beta-3 Cytoplasmic Tail bound to Talin-1 F3 FERM Domain (PDB: 1MK7)",
    },
    {
        "seq": "HPHPH",
        "title": "Histidine-Proline Cu2+/Zn2+ Antimicrobial Repeat",
        "pattern": "H-P-H-P-H",
        "family": "Metal Chelation / Host Defense",
        "hint_role": "Conformationally restricted metal-chelating motif found in histidine-rich glycoproteins.",
        "hint_detail": "Strictly alternating His and Pro residues.",
        "context_requirement": "Tandemly repeated in the central histidine-proline-rich (HPR) domain of plasma Histidine-Rich Glycoprotein (HRG) for Zn2+/Cu2+ and heme sequestration.",
        "exemplar_source": "Human Histidine-Rich Glycoprotein (HRG, UniProt: P04196)",
    },
]


def build_all() -> Dict[str, Any]:
    """Compiles, validates, and audits the complete Chemle + Seqle puzzle database."""
    chemle_records: List[Dict[str, Any]] = []
    chemle_fps: List[Any] = []
    seen_ids: set[str] = set()
    seen_inchikeys: set[str] = set()
    seen_smiles: set[str] = set()

    print(
        f"[BioWordle Build] Auditing & compiling {len(CHEMLE_RAW_COMPOUNDS)} Chemle molecules "
        f"with RDKit {Chem.rdBase.rdkitVersion}..."
    )
    for item in CHEMLE_RAW_COMPOUNDS:
        cid = item["id"]
        if cid in seen_ids:
            raise ValueError(f"Duplicate Chemle ID: {cid}")
        seen_ids.add(cid)

        mol, audit = standardize_and_audit(item["smiles"])
        if audit["inchikey"] in seen_inchikeys:
            raise ValueError(f"Duplicate InChIKey after standardization for {cid}: {audit['inchikey']}")
        if audit["canonical_smiles"] in seen_smiles:
            raise ValueError(f"Duplicate canonical SMILES for {cid}: {audit['canonical_smiles']}")
        seen_inchikeys.add(audit["inchikey"])
        seen_smiles.add(audit["canonical_smiles"])

        desc = compute_descriptors_and_insights(mol)
        fg_list = detect_functional_groups(mol)

        svg_full = render_responsive_svg(mol, width=320, height=220, dark_mode=False)
        svg_scaffold, scaffold_smiles, scaffold_mode = build_scaffold_clue_svg(mol)

        fp = compute_chiral_morgan_fp(mol)
        chemle_fps.append(fp)

        record = {
            "id": cid,
            "name": item["name"],
            "Category": item["Category"],
            **audit,
            **desc,
            "functional_groups": fg_list,
            "hint_role": item["hint_role"],
            "hint_origin": item["hint_origin"],
            "structural_note": item["structural_note"],
            "scaffold_smiles": scaffold_smiles,
            "scaffold_mode": scaffold_mode,
            "has_scaffold": scaffold_mode in ("bemis_murcko", "ring_core", "carbon_skeleton"),
            "svg_scaffold": svg_scaffold,
            "svg_full": svg_full,
        }
        chemle_records.append(record)

    # Compute N x N pairwise Chiral Morgan Tanimoto similarity matrix (0 - 100%)
    n_chem = len(chemle_records)
    similarity_matrix: Dict[str, Dict[str, int]] = {}
    for i in range(n_chem):
        id_i = chemle_records[i]["id"]
        similarity_matrix[id_i] = {}
        for j in range(n_chem):
            id_j = chemle_records[j]["id"]
            if i == j:
                similarity_matrix[id_i][id_j] = 100
            else:
                sim = DataStructs.TanimotoSimilarity(chemle_fps[i], chemle_fps[j])
                similarity_matrix[id_i][id_j] = int(round(sim * 100))

    # Attach Shannon Entropy Information Gain & Difficulty metrics to Chemle
    chem_info = compute_chemle_information_metrics(chemle_records, similarity_matrix)
    for rec in chemle_records:
        rec["info_metrics"] = chem_info[rec["id"]]

    print(
        f"[BioWordle Build] Auditing & compiling {len(SEQLE_RAW_MOTIFS)} Seqle 5-AA motifs..."
    )
    seqle_records: List[Dict[str, Any]] = []
    seen_seqs: set[str] = set()
    for item in SEQLE_RAW_MOTIFS:
        seq = item["seq"].strip().upper()
        if len(seq) != 5 or any(ch not in AA_PROPERTIES for ch in seq):
            raise ValueError(f"Invalid 5-AA sequence: {seq}")
        if seq in seen_seqs:
            raise ValueError(f"Duplicate Seqle sequence: {seq}")
        seen_seqs.add(seq)

        pep_mol = Chem.MolFromSequence(seq)
        if pep_mol is None:
            raise ValueError(f"RDKit MolFromSequence failed for: {seq}")

        pep_desc = compute_descriptors_and_insights(pep_mol)
        biochem = compute_peptide_biochem(seq, AA_PROPERTIES)
        svg_peptide = render_responsive_svg(pep_mol, width=380, height=210, dark_mode=False)

        seqle_records.append(
            {
                "seq": seq,
                "title": item["title"],
                "pattern": item["pattern"],
                "family": item["family"],
                "hint_role": item["hint_role"],
                "hint_detail": item["hint_detail"],
                "context_requirement": item["context_requirement"],
                "exemplar_source": item["exemplar_source"],
                "formula": pep_desc["formula"],
                "MW": pep_desc["MW"],
                "cLogP": pep_desc["cLogP"],
                "TPSA": pep_desc["TPSA"],
                "HBD": pep_desc["HBD"],
                "HBA": pep_desc["HBA"],
                "HeavyAtoms": pep_desc["HeavyAtoms"],
                **biochem,
                "svg_peptide": svg_peptide,
            }
        )

    seq_info = compute_seqle_information_metrics(seqle_records, AA_PROPERTIES)
    for rec in seqle_records:
        rec["info_metrics"] = seq_info[rec["seq"]]

    # Compute deterministic dataset version hash
    hash_input = "|".join(
        sorted(r["inchikey"] for r in chemle_records)
        + sorted(r["seq"] for r in seqle_records)
    )
    dataset_hash = hashlib.sha256(hash_input.encode("utf-8")).hexdigest()[:12]

    return {
        "metadata": {
            "generator": "BioWordle Independent RDKit Build Compiler (chem_engine.py)",
            "rdkit_version": Chem.rdBase.rdkitVersion,
            "dataset_version_hash": dataset_hash,
            "daily_seed_policy": "UTC calendar date (YYYY-MM-DD) combined with dataset_version_hash",
            "chemle_count": len(chemle_records),
            "seqle_count": len(seqle_records),
            "fingerprint_spec": FINGERPRINT_SPEC,
            "standardization_spec": STANDARDIZATION_SPEC,
            "peptide_pka_spec": PEPTIDE_PKA_SPEC,
        },
        "aa_properties": AA_PROPERTIES,
        "chemle": chemle_records,
        "chemle_similarity": similarity_matrix,
        "seqle": seqle_records,
    }


def main() -> None:
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, "data")
    os.makedirs(data_dir, exist_ok=True)

    payload = build_all()

    json_path = os.path.join(data_dir, "puzzles.json")
    tmp_json = json_path + ".tmp"
    with open(tmp_json, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    os.replace(tmp_json, json_path)

    js_path = os.path.join(data_dir, "puzzles_data.js")
    tmp_js = js_path + ".tmp"
    with open(tmp_js, "w", encoding="utf-8") as f:
        f.write("// Auto-generated by build_puzzles.py using RDKit — do not edit manually.\n")
        f.write("window.BIOWORDLE_DATA = ")
        json.dump(payload, f, separators=(",", ":"))
        f.write(";\n")
    os.replace(tmp_js, js_path)

    print(
        f"[BioWordle Build] Done! Dataset hash: {payload['metadata']['dataset_version_hash']}\n"
        f"  - {json_path}\n"
        f"  - {js_path}"
    )

    index_path = os.path.join(base_dir, "index.html")
    css_path = os.path.join(base_dir, "styles.css")
    app_js_path = os.path.join(base_dir, "app.js")
    if all(os.path.isfile(p) for p in (index_path, css_path, app_js_path)):
        with open(index_path, "r", encoding="utf-8") as f:
            html_src = f.read()
        with open(css_path, "r", encoding="utf-8") as f:
            css_src = f.read()
        with open(js_path, "r", encoding="utf-8") as f:
            data_js_src = f.read()
        with open(app_js_path, "r", encoding="utf-8") as f:
            app_js_src = f.read()

        bundled = html_src.replace(
            '<link rel="stylesheet" href="styles.css" />',
            f"<style>\n{css_src}\n</style>",
        )
        bundled = bundled.replace(
            '<script src="data/puzzles_data.js"></script>\n  <script src="app.js"></script>',
            f"<script>\n{data_js_src}\n</script>\n<script>\n{app_js_src}\n</script>",
        )

        dist_dir = os.path.join(base_dir, "dist")
        os.makedirs(dist_dir, exist_ok=True)
        standalone_path = os.path.join(dist_dir, "biowordle_standalone.html")
        tmp_standalone = standalone_path + ".tmp"
        with open(tmp_standalone, "w", encoding="utf-8") as f:
            f.write(bundled)
        os.replace(tmp_standalone, standalone_path)
        print(f"  - {standalone_path} (Single-file standalone bundle)")

        if len(sys.argv) > 1 and sys.argv[1] == "--artifact-output" and len(sys.argv) > 2:
            art_path = sys.argv[2]
            tmp_art = art_path + ".tmp"
            with open(tmp_art, "w", encoding="utf-8") as f:
                f.write(bundled)
            os.replace(tmp_art, art_path)
            print(f"  - {art_path} (Artifact preview updated)")


if __name__ == "__main__":
    main()
