# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "pandas",
#     "numpy",
# ]
# ///
"""Build-time RDKit Puzzle Compiler for BioWordle v0.2 ("Restraint" Release).

Implements all P0 and P1 patch requirements from BIOWORDLE_PATCH.md:
- P0-1: k-anonymous Level-0 opening clues (k_min >= 4 Ranked, >= 3 Learn) + post-game probe quality
- P0-2: Solver-bot build gate verifying mean guesses, <=2-guess rate, and min expected remaining
- P0-3: Separation of answerPool and guessPool metadata
- P0-4: Seqle motif schema (kind, consensus, parent_sequence, window_offset, terminal, exemplar_source, game_eligibility) + validator
- P1-1: Shared `data/feedback_spec.json` + cross-language `tests/golden_feedback.json`
- P1-3: Order-independent canonical dataset hash over all gameplay fields
- P1-5: Authentic RDKit SetACS1996Mode rendering, split carbonyl SMARTS, normalized SVG floats
- P2-3: Automated hint-text non-leakage audit
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from rdkit import Chem, DataStructs

from chem_engine import (
    FINGERPRINT_SPEC,
    PEPTIDE_PKA_SPEC,
    PROVENANCE_VERIFIED_DATE,
    SEQLE_REFERENCE_REGISTRY,
    STANDARDIZATION_SPEC,
    UNIPROT_VERIFIED_RELEASE,
    assert_seqle_solver_gate,
    assert_solver_gate,
    audit_hint_text_leakage,
    build_scaffold_clue_svg,
    compute_chiral_morgan_fp,
    compute_descriptors_and_insights,
    compute_peptide_biochem,
    compute_seqle_information_metrics,
    detect_functional_groups,
    evaluate_chemle_feedback_from_spec,
    evaluate_seqle_feedback_from_spec,
    render_responsive_svg,
    run_chemle_greedy_solver,
    run_seqle_greedy_solver,
    select_k_anonymous_level0_features,
    standardize_and_audit,
    validate_chemle_manifest,
    validate_seqle_motif_record,
    validate_seqle_registry_completeness,
)

SCHEMA_VERSION = "0.2.2"

# ---------------------------------------------------------------------------
# 1. Curated CHEMLE Compound Library (36 Molecules, P2-3 Leakage-Audited)
# ---------------------------------------------------------------------------
SUPER_CLASS_MAP: Dict[str, str] = {
    "Alkaloid / Metabolite": "Bioactive Nitrogen Heterocycle / Amine",
    "Neurotransmitter": "Bioactive Nitrogen Heterocycle / Amine",
    "Hormone / Neurotransmitter": "Bioactive Nitrogen Heterocycle / Amine",
    "Hormone / Metabolite": "Bioactive Nitrogen Heterocycle / Amine",
    "Biogenic Amine": "Bioactive Nitrogen Heterocycle / Amine",
    "Primary Metabolite": "Core Cellular Metabolite / Cofactor",
    "Cofactor / Vitamin": "Core Cellular Metabolite / Cofactor",
    "Hydrocarbon Degradation": "Microbial & Aromatic Catabolite",
    "Microbial Natural Product": "Microbial & Aromatic Catabolite",
    "Plant Natural Product": "Botanical Secondary Metabolite / Drug",
    "Therapeutic Drug": "Botanical Secondary Metabolite / Drug",
}

_DEFAULT_MANIFEST_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "data", "chemle_compounds.json"
)


def load_chemle_manifest(manifest_path: Optional[str] = None) -> List[Dict[str, Any]]:
    """Loads and offline-validates the curated Chemle compound manifest."""
    path = manifest_path or _DEFAULT_MANIFEST_PATH
    with open(path, "r", encoding="utf-8") as fh:
        payload = json.load(fh)
    compounds = payload.get("compounds") if isinstance(payload, dict) else payload
    validate_chemle_manifest(compounds)
    return list(compounds)


CHEMLE_RAW_COMPOUNDS: List[Dict[str, Any]] = load_chemle_manifest()


# ---------------------------------------------------------------------------
# 2. Curated SEQLE 5-Residue Motif Library (30 Motifs, P0-4 & P2-3 Audited)
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

SEQLE_RAW_MOTIFS: List[Dict[str, Any]] = [
    {
        "seq": "GHSMG",
        "title": "Serine Hydrolase / Cutinase GXSXG Elbow",
        "kind": "exact",
        "consensus": "G-X-S-X-G",
        "parent_sequence": "AVAGHSMGGGG",
        "window_offset": 3,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Catalytic Active Site",
        "hint_role": "Nucleophilic beta-turn of alpha/beta-hydrolases, lipases, PETases, and esterases.",
        "hint_detail": "Flanking achiral residues prevent steric clashes at the sharp turn carrying the nucleophilic alcohol side chain.",
        "context_requirement": "Embedded at a tight beta-strand-to-alpha-helix nucleophilic elbow coupled to a His-Asp catalytic dyad.",
        "exemplar_source": "UniProt:G9BY57 (Leaf-branch compost cutinase LCC 163-167); PDB:4EB0",
    },
    {
        "seq": "EKDEL",
        "title": "ER Luminal C-Terminal Retention Signal",
        "kind": "exact",
        "consensus": "[DE]-[KRH]-[DENQ]-E-L",
        "parent_sequence": "TAEKDEL",
        "window_offset": 2,
        "terminal": "C",
        "game_eligibility": ["learn", "ranked"],
        "family": "Subcellular Targeting",
        "hint_role": "Recognized in the cis-Golgi by seven-transmembrane retrieval receptors to return escaped ER chaperones.",
        "hint_detail": "Extreme C-terminal tetrapeptide preceded by an acidic glutamate in human ER Hsp70 chaperone BiP.",
        "context_requirement": "Must reside at the extreme C-terminus (-COOH free) of soluble ER luminal proteins (BiP/HSPA5, PDI).",
        "exemplar_source": "UniProt:P11021 (Human HSPA5/BiP C-terminus 650-654); PDB:6I6H",
    },
    {
        "seq": "DEVDG",
        "title": "Apoptotic Executioner Protease Cleavage Site",
        "kind": "exact",
        "consensus": "D-E-V-D-G",
        "parent_sequence": "GDEVDGVDEV",
        "window_offset": 1,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Protease Cleavage Motif",
        "hint_role": "Canonical P4-P1' recognition sequence cleaved after the second acidic residue in PARP-1 during apoptosis.",
        "hint_detail": "Strictly requires aspartate at both P4 and P1 positions followed by a small unbranched P1' residue.",
        "context_requirement": "Requires an exposed solvent-accessible loop between the N-terminal DNA-binding and C-terminal catalytic domains.",
        "exemplar_source": "UniProt:P09874 (Human PARP1 211-215); ELM:CLV_Caspase3-7",
    },
    {
        "seq": "GRGDS",
        "title": "Fibronectin III-10 Integrin-Binding Loop",
        "kind": "exact",
        "consensus": "G-R-G-D-S",
        "parent_sequence": "VTGRGDSPASS",
        "window_offset": 2,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Cell Adhesion",
        "hint_role": "Extracellular matrix cell-attachment loop recognized by alpha5-beta1 and alphaV-beta3 integrins.",
        "hint_detail": "Flanks the central cationic-glycine-anionic tripeptide with an N-terminal Gly and C-terminal polar Ser.",
        "context_requirement": "Presented at the apex of the 10th fibronectin type III beta-hairpin loop; Asp coordinates the MIDAS Mg2+ ion.",
        "exemplar_source": "UniProt:P02751 (Human Fibronectin FN1 1614-1618; mature 1492-1496); PDB:1FNF",
    },
    {
        "seq": "LPETG",
        "title": "Sortase A Gram-Positive Anchor Motif",
        "kind": "exact",
        "consensus": "L-P-X-T-G",
        "parent_sequence": "QALPETGEEN",
        "window_offset": 2,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Transpeptidation Tag",
        "hint_role": "Cleaved between the fourth and fifth residues by a bacterial cysteine transpeptidase and ligated to peptidoglycan cross-bridges.",
        "hint_detail": "Widely engineered for site-specific enzymatic protein bioconjugation via a thioacyl intermediate.",
        "context_requirement": "Located upstream of a hydrophobic transmembrane helix and positively charged C-terminal tail in surface proteins.",
        "exemplar_source": "PDB:1T2W; UniProt:P02976 (Staphylococcus aureus Protein A sorting signal)",
    },
    {
        "seq": "HEIGH",
        "title": "Zincin Metalloprotease Catalytic Helix",
        "kind": "exact",
        "consensus": "H-E-X-X-H",
        "parent_sequence": "VAAHEIGHSLGL",
        "window_offset": 3,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Metal-Binding Catalytic Site",
        "hint_role": "Two imidazole side chains coordinate the active-site Zn2+ ion while glutamate polarizes the attacking water.",
        "hint_detail": "Conserved active-site signature across matrixins (MMPs), thermolysin, and neurotoxins.",
        "context_requirement": "Must fold as an active-site alpha-helix so i and i+4 His imidazole rings project onto the same face toward Zn2+.",
        "exemplar_source": "UniProt:P08254 (Human MMP3 Stromelysin-1 218-222); PDB:1B3D; PROSITE:PS00142",
    },
    {
        "seq": "CCHHC",
        "title": "Retroviral Zinc Knuckle Ligand Abstraction",
        "kind": "abstraction",
        "consensus": "C-X2-C-X4-H-X4-C",
        "parent_sequence": "CCHHC",
        "window_offset": 0,
        "terminal": "none",
        "game_eligibility": ["learn"],
        "family": "Nucleic Acid Binding",
        "hint_role": "Compressed pedagogical representation of Cys/His tetrahedral Zn2+ coordination in viral nucleocapsid domains.",
        "hint_detail": "Illustrates the sulfur- and imidazole-rich metal-chelation chemistry of viral gag CCHC motifs (Learn mode only).",
        "context_requirement": "In native proteins, ligands are spaced across 14 residues (C-X2-C-X4-H-X4-C) rather than a contiguous pentapeptide.",
        "exemplar_source": "PDB:1MFS; PROSITE:PS50158 (Retroviral Gag CCHC Zinc Knuckle)",
    },
    {
        "seq": "YIGSR",
        "title": "Laminin Beta-1 Epithelial Adhesion Pentapeptide",
        "kind": "exact",
        "consensus": "Y-I-G-S-R",
        "parent_sequence": "CDPGYIGSRCD",
        "window_offset": 4,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Cell Adhesion",
        "hint_role": "Basement membrane glycoprotein sequence promoting cell attachment via the 67-kDa high-affinity receptor.",
        "hint_detail": "Begins with a phenolic aromatic residue and beta-branched isoleucine, ending with a guanidinium side chain.",
        "context_requirement": "Located within cysteine-rich EGF-like repeat 9 of the laminin subunit beta-1 short arm.",
        "exemplar_source": "UniProt:P07942 (Human Laminin Subunit Beta-1 950-954; mature 929-933)",
    },
    {
        "seq": "NPVYK",
        "title": "Clathrin Endocytosis NPXY Internalization Turn",
        "kind": "exact",
        "consensus": "N-P-X-Y-X",
        "parent_sequence": "NPVYK",
        "window_offset": 0,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Membrane Trafficking",
        "hint_role": "Cytoplasmic tail sequence of lipoprotein and growth-factor receptors recruiting PTB-domain coat adaptors.",
        "hint_detail": "Conserved Asn-Pro-X-Tyr turn docks the unphosphorylated phenolic ring into DAB2/ARH PTB adaptors (teaching example suffix X=Lys).",
        "context_requirement": "Synthetic consensus example for the F-X-N-P-X-Y internalization signal (human LDLR UniProt:P01130 825-829 carries NPVYQ in PDB:3SO6).",
        "exemplar_source": "ELM:TRG_ENDOCYTIC_2; PDB:3SO6 (motif-family context only: human LDLR UniProt:P01130 carries NPVYQ)",
    },
    {
        "seq": "KCVIM",
        "title": "KRAS4B C-Terminal CAAX Prenylation Tail",
        "kind": "fragment",
        "consensus": "[KR]-C-[AVILM]-[AVILM]-[MCSQ]",
        "parent_sequence": "KKKKKKSKTKCVIM",
        "window_offset": 9,
        "terminal": "C",
        "game_eligibility": ["learn", "ranked"],
        "family": "Post-Translational Lipidation",
        "hint_role": "Farnesyltransferase attaches a C15 isoprenoid lipid to the fourth-from-last thiol before tripeptide proteolysis.",
        "hint_detail": "Preceded by a polybasic lysine track; final methionine residue directs high-affinity farnesylation.",
        "context_requirement": "Strictly requires extreme C-terminal positioning; -VIM is cleaved by Rce1 and alpha-carboxyl-methylated by ICMT.",
        "exemplar_source": "UniProt:P01116-2 (Human KRAS Isoform 2B C-terminus 184-188); PDB:1D8D",
    },
    {
        "seq": "DEADE",
        "title": "RNA Helicase Walker B Mg2+-Binding Box",
        "kind": "fragment",
        "consensus": "V-L-D-E-A-D-E-M-L",
        "parent_sequence": "VLDEADEMLSR",
        "window_offset": 2,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "NTPase / Helicase Core",
        "hint_role": "Coordinates divalent metal via aspartate and polarizes the hydrolytic water molecule in ATP-driven duplex unwinding.",
        "hint_detail": "Five-residue window capturing the signature acidic-acidic-aliphatic-acidic-acidic box in human eIF4A1.",
        "context_requirement": "Positioned at the C-terminal tip of RecA-like domain 1 strand beta-3 adjacent to the phosphate-binding loop.",
        "exemplar_source": "UniProt:P60842 (Human eIF4A1 / DDX2A 182-186); PDB:2ZU6",
    },
    {
        "seq": "GVGKS",
        "title": "Walker A / P-Loop Phosphate-Binding Tail",
        "kind": "fragment",
        "consensus": "G-X4-G-K-[TS]",
        "parent_sequence": "GAGGVGKS",
        "window_offset": 3,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "NTPase / Helicase Core",
        "hint_role": "Distal half of the glycine-rich nucleotide loop that wraps around beta/gamma-phosphates of GTP and ATP.",
        "hint_detail": "Displays the conserved C-terminal triad where an invariant primary amine neutralizes nucleotide polyphosphate negative charge.",
        "context_requirement": "Connects strand beta-1 and helix alpha-1 in Rossmann-fold GTPases and ATPases.",
        "exemplar_source": "UniProt:P01116 (Human KRAS P-loop 13-17); PROSITE:PS00017; PDB:4OBE",
    },
    {
        "seq": "VDFGL",
        "title": "Kinase Activation Segment Mg2+ Switch",
        "kind": "exact",
        "consensus": "[VI]-D-F-G-[LM]",
        "parent_sequence": "KAVLVDFGLSVQM",
        "window_offset": 4,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Kinase Switch",
        "hint_role": "Conformational flip (in vs. out) at the start of the regulatory loop governing ATP-metal catalysis.",
        "hint_detail": "Central aspartate chelates divalent magnesium while the aromatic ring packs into the regulatory R-spine; targeted by type-II kinase inhibitors.",
        "context_requirement": "Marks the N-terminal anchor of the eukaryotic protein kinase activation loop between beta-8 and beta-9.",
        "exemplar_source": "UniProt:P41279 (Human MAP3K8/COT Kinase 269-273); PDB:4Y83",
    },
    {
        "seq": "HRDLA",
        "title": "Tyrosine Kinase Catalytic Loop Base",
        "kind": "exact",
        "consensus": "H-R-D-L-A",
        "parent_sequence": "IHRDLAARN",
        "window_offset": 1,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Catalytic Active Site",
        "hint_role": "Central aspartate acts as general proton acceptor for the substrate hydroxyl group during gamma-phosphoryl transfer.",
        "hint_detail": "Cationic arginine neutralizes the phosphorylated regulatory loop to align the active site.",
        "context_requirement": "Located in the C-lobe catalytic loop between beta-6 and beta-7 strands directly facing the ATP gamma-phosphate.",
        "exemplar_source": "UniProt:P10721 (Human KIT Receptor Tyrosine Kinase 790-794); PDB:1PKG; PROSITE:PS00109",
    },
    {
        "seq": "WSEWS",
        "title": "Type I Cytokine Receptor Folding Zipper",
        "kind": "exact",
        "consensus": "W-S-X-W-S",
        "parent_sequence": "YNGIWSEWSEARS",
        "window_offset": 4,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Structural Folding Motif",
        "hint_role": "Forms a Tryptophan-Arginine pi-cation stack essential for ER maturation and export of hematopoietin family proteins.",
        "hint_detail": "Symmetric pair of indole-serine dipeptides separated by a central acidic glutamate.",
        "context_requirement": "Located in the membrane-proximal fibronectin type III domain of Class I cytokine receptors; often C-mannosylated.",
        "exemplar_source": "UniProt:P32927 (Human Cytokine Receptor Common Beta Chain CSF2RB 425-429); PDB:1EGJ; PROSITE:PS01355",
    },
    {
        "seq": "PPPVP",
        "title": "Polyproline Type II Helix SH3 Docking Tract",
        "kind": "exact",
        "consensus": "P-P-X-[VIL]-P",
        "parent_sequence": "VPPPVPPRRR",
        "window_offset": 1,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Protein-Protein Interaction",
        "hint_role": "Left-handed 3.0-residues-per-turn secondary structure binding into aromatic Trp/Tyr grooves of Src-homology-3 domains.",
        "hint_detail": "Four cyclic imino acid residues flanking a single beta-branched aliphatic valine.",
        "context_requirement": "Adopts an extended left-handed PPII helical conformation recognized by Class II SH3 adaptor domains.",
        "exemplar_source": "UniProt:Q07889 (Human SOS1 Proline-Rich Domain 1150-1154); ELM:LIG_SH3_2; PDB:1GBR (motif-family context only: murine Sos1 Q02384)",
    },
    {
        "seq": "INFEK",
        "title": "Ovalbumin H-2Kb MHC Class I Central Core",
        "kind": "fragment",
        "consensus": "S-I-I-N-F-E-K-L",
        "parent_sequence": "SIINFEKL",
        "window_offset": 2,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Immunology Epitope",
        "hint_role": "Internal P3-P7 segment of the benchmark CD8+ T-cell octapeptide presented by murine major histocompatibility complexes.",
        "hint_detail": "Aromatic Phe at P5 buries into the deep C-pocket while P4 Asn and P6 Glu/P7 Lys face the T-cell receptor.",
        "context_requirement": "Internal 5-mer window (residues 3-7) of the 8-residue SIINFEKL peptide bound in the MHC-I groove.",
        "exemplar_source": "UniProt:P01012 (Gallus gallus Ovalbumin 260-264 within 258-265); PDB:1VAC",
    },
    {
        "seq": "YQRLN",
        "title": "Clathrin AP-2 Mu2 Tyrosine Sorting Signal",
        "kind": "exact",
        "consensus": "Y-X-X-[LIVF]-N",
        "parent_sequence": "SDYQRLNLKL",
        "window_offset": 2,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Membrane Trafficking",
        "hint_role": "Binds the medium subunit of heterotetrameric vesicle adaptor complexes for rapid endocytosis from the plasma membrane.",
        "hint_detail": "Requires an unphosphorylated N-terminal phenolic residue and a bulky hydrophobic side chain three positions downstream.",
        "context_requirement": "Must be in an extended conformation within a cytosolic receptor tail; Tyr phosphorylation blocks AP-2 binding.",
        "exemplar_source": "UniProt:P19814 (Rat TGN38 Cytoplasmic Tail 350-354); PDB:1BXX; ELM:TRG_ENDOCYTIC_1",
    },
    {
        "seq": "KEKSD",
        "title": "COPI Coatomer Dilysine ER Retrieval Signal",
        "kind": "exact",
        "consensus": "K-X-K-X-X",
        "parent_sequence": "LGLFIFSIVFLHMKEKEKSD",
        "window_offset": 15,
        "terminal": "C",
        "game_eligibility": ["learn", "ranked"],
        "family": "Subcellular Targeting",
        "hint_role": "Cytosolic extreme C-terminal retrieval signal on ER-resident membrane proteins recognized by alpha/beta'-COP WD40 propellers.",
        "hint_detail": "Alternating cationic primary amines at positions -5 and -3 from the C-terminus pair with an acidic side chain, polar hydroxyl, and terminal carboxylate.",
        "context_requirement": "Canonical extreme C-terminal K-X-K-X-X retrieval signal (lysines at positions -5 and -3; co-crystallized AKEKSD peptide in PDB:4J77 entity 1).",
        "exemplar_source": "UniProt:P39656 (Human OST48/DDOST 452-456); PDB:4J77 (AKEKSD entity 1); ELM:TRG_ER_diLys_1",
    },
    {
        "seq": "CAACH",
        "title": "Cytochrome c CxxCH Heme Attachment Site",
        "kind": "exact",
        "consensus": "C-X-X-C-H",
        "parent_sequence": "YREVCAACHSLDR",
        "window_offset": 4,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Cofactor Covalent Attachment",
        "hint_role": "Two cysteine thiols form covalent thioether bonds to porphyrin vinyl groups while histidine serves as axial Fe ligand.",
        "hint_detail": "Defining bioinorganic signature of electron-transport hemoproteins in mitochondrial and bacterial respiratory chains.",
        "context_requirement": "Matured by cytochrome c heme lyase / Ccm system in the periplasm or mitochondrial intermembrane space.",
        "exemplar_source": "UniProt:P07143 (Yeast Mitochondrial Cytochrome c1 101-105); PDB:3CX5; PROSITE:PS51007",
    },
    {
        "seq": "NLTGY",
        "title": "ER Luminal N-Linked Glycosylation Sequon",
        "kind": "exact",
        "consensus": "N-[^P]-[ST]-X-X",
        "parent_sequence": "NLTGY",
        "window_offset": 0,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Post-Translational Modification",
        "hint_role": "Oligosaccharyltransferase transfers a 14-sugar Glc3Man9GlcNAc2 precursor onto the primary carboxamide nitrogen.",
        "hint_detail": "Requires Asn at position 1, any residue except Pro at position 2, and Thr or Ser at position 3 (synthetic G-Y padding completes the 5-mer).",
        "context_requirement": "Synthetic N-[^P]-[ST] eukaryotic acceptor sequon example (no universal Tyr5 requirement; PDB:3WAK is apo archaeal AglB UniProt:O29867).",
        "exemplar_source": "ELM:MOD_N-GLC_1; PDB:3WAK (motif-family context only: apo Archaeoglobus fulgidus AglB UniProt:O29867)",
    },
    {
        "seq": "DSGLD",
        "title": "SCF-Beta-TrCP Phosphodegron N-Core",
        "kind": "fragment",
        "consensus": "D-S-G-X-X-S",
        "parent_sequence": "DDRHDSGLDSMKD",
        "window_offset": 4,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Proteasomal Degradation",
        "hint_role": "Phosphorylation by IKK kinase recruits the cullin-RING E3 ubiquitin ligase to trigger NF-kappa-B activation.",
        "hint_detail": "First five residues of the human I-kappa-B-alpha destruction box bound to the SCF-beta-TrCP WD40 receptor in PDB:9T9W.",
        "context_requirement": "Requires dual serine phosphorylation within the disordered N-terminal tail (31-35) of I-kappa-B-alpha.",
        "exemplar_source": "UniProt:P25963 (Human NFKBIA / I-kappa-B-alpha 31-35); PDB:9T9W (entity 2); ELM:DEG_SCF_TRCP1_1",
    },
    {
        "seq": "KENIP",
        "title": "APC/C-Cdh1 Mitotic Exit Destruction Box",
        "kind": "fragment",
        "consensus": "K-E-N-X-X",
        "parent_sequence": "KENIP",
        "window_offset": 0,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Proteasomal Degradation",
        "hint_role": "Recruits the Anaphase-Promoting Complex coactivator to ubiquitinate cell-cycle kinases during late mitosis and G1.",
        "hint_detail": "Begins with an invariant cationic-anionic-carboxamide core followed by an aliphatic-imino teaching suffix (downstream positions vary across native substrates).",
        "context_requirement": "Synthetic K-E-N degron example without a universal Pro requirement (PDB:4BH6 Acm1 UniProt:Q08981 carries KENLS; Hsl1 UniProt:P34244 carries NKENE).",
        "exemplar_source": "ELM:DEG_APCC_KENBOX_2; PDB:4BH6 (motif-family context only: yeast Acm1 UniProt:Q08981 carries KENLS)",
    },
    {
        "seq": "KKKRK",
        "title": "SV40 Large T-Antigen Classical NLS Core",
        "kind": "fragment",
        "consensus": "P-K-K-K-R-K-V",
        "parent_sequence": "PKKKRKV",
        "window_offset": 1,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Subcellular Targeting",
        "hint_role": "Polycationic five-residue segment of the prototype monopartite nuclear import signal recognized by Importin-alpha.",
        "hint_detail": "Five consecutive positively charged basic residues (four primary amines and one guanidinium) docking into ARM pockets P1-P5.",
        "context_requirement": "Binds in an extended conformation across the major groove of Karyopherin/Importin-alpha (residues 127-131 of PKKKRKV).",
        "exemplar_source": "UniProt:P03070 (SV40 Large T-Antigen 127-131); PDB:1BK6; ELM:TRG_NLS_MonoCore_2",
    },
    {
        "seq": "FGGFT",
        "title": "Orphanin FQ / Nociceptin Receptor Agonist",
        "kind": "fragment",
        "consensus": "F-G-G-F-T-G-A-R-K-S",
        "parent_sequence": "FGGFTGARKSARKLANQ",
        "window_offset": 0,
        "terminal": "N",
        "game_eligibility": ["learn", "ranked"],
        "family": "Neuropeptide GPCR Message",
        "hint_role": "Extreme N-terminal 'message' pentapeptide activating the NOP opioid-related GPCR while excluding classical mu/delta subtypes.",
        "hint_detail": "Replaces the N-terminal phenolic Tyr1 of enkephalins with an unhydroxylated Phe1, a flexible achiral dipeptide spacer, Phe4, and Thr5.",
        "context_requirement": "Requires a free N-terminal alpha-amino group on Phe1 to form a salt bridge with Asp130 in the NOP receptor pocket.",
        "exemplar_source": "UniProt:Q13519 (Human Prepronociceptin PNOC 130-134; mature Nociceptin 130-146)",
    },
    {
        "seq": "YGGFM",
        "title": "Met-Enkephalin Endogenous Opioid Peptide",
        "kind": "exact",
        "consensus": "Y-G-G-F-M",
        "parent_sequence": "YGGFM",
        "window_offset": 0,
        "terminal": "N",
        "game_eligibility": ["learn", "ranked"],
        "family": "Neuropeptide GPCR Message",
        "hint_role": "Complete native delta- and mu-receptor agonist pentapeptide cleaved from proenkephalin-A.",
        "hint_detail": "N-terminal Tyr1 phenol and free alpha-amine mimic the tyramine pharmacophore of morphine; ends with methionine.",
        "context_requirement": "Functions as an intact free pentapeptide neurotransmitter released by prohormone convertases.",
        "exemplar_source": "UniProt:P01210 (Human Proenkephalin-A PENK 100-104); PDB:1PLW",
    },
    {
        "seq": "SYSGY",
        "title": "FUS Low-Complexity Prion-Like Condensate Repeat",
        "kind": "fragment",
        "consensus": "[SG]-Y-[SG]-G-Y",
        "parent_sequence": "GQQSYSGYSQ",
        "window_offset": 3,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Biomolecular Condensate",
        "hint_role": "Pi-pi and sp2/pi sticker-and-spacer motif driving liquid-liquid phase separation of ribonucleoprotein granules.",
        "hint_detail": "Combines polar serine and flexible glycine spacers with two phenolic tyrosine 'stickers' in nuclear RNA-binding proteins.",
        "context_requirement": "Tandemly repeated across the intrinsically disordered N-terminal low-complexity domain (residues 1-214) of FUS.",
        "exemplar_source": "UniProt:P35637 (Human RNA-Binding Protein FUS 37-41); PDB:5W3N",
    },
    {
        "seq": "LDLAS",
        "title": "NES Amphipathic Export Ridge Abstraction",
        "kind": "abstraction",
        "consensus": "L-X(2,3)-L-X(2,3)-L-X-L",
        "parent_sequence": "LDLAS",
        "window_offset": 0,
        "terminal": "none",
        "game_eligibility": ["learn"],
        "family": "Subcellular Targeting",
        "hint_role": "Simplified 5-residue teaching representation of a CRM1-binding nuclear translocation signal (Learn mode only).",
        "hint_detail": "Shows two hydrophobic branched side chains flanking an acidic aspartate followed by a small methyl and hydroxyl tail.",
        "context_requirement": "Native CRM1/XPO1 nuclear export signals span 10-12 residues (four spaced Phi positions) rather than a single pentapeptide.",
        "exemplar_source": "PDB:3NBY; ELM:TRG_NES_CRM1_1 (CRM1-RanGTP Nuclear Export Complex)",
    },
    {
        "seq": "NPLYK",
        "title": "Integrin Beta-3 Talin-Binding Cytoplasmic Turn",
        "kind": "exact",
        "consensus": "N-P-L-Y-K",
        "parent_sequence": "AKWDTANNPLYKEATSTFTNITYRGT",
        "window_offset": 7,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Cell Adhesion / Inside-Out Signaling",
        "hint_role": "Binds the F3 FERM domain of cytoskeletal adaptors to switch adhesion receptors from bent-inactive to extended-active states.",
        "hint_detail": "Carboxamide-proline-leucine-phenol beta-turn capped by a basic primary-amine side chain in membrane-proximal receptor tails.",
        "context_requirement": "Located in the membrane-proximal cytoplasmic tail of Integrin Beta-3; Tyr phosphorylation switches binding to Dok1.",
        "exemplar_source": "UniProt:P05106 (Human Integrin Beta-3 ITGB3 770-774; mature 744-748); PDB:1MK7",
    },
    {
        "seq": "GHHPH",
        "title": "Plasma Glycoprotein Zn2+/Cu2+ Chelation Repeat",
        "kind": "exact",
        "consensus": "G-H-H-P-H",
        "parent_sequence": "GHHPHGHHPH",
        "window_offset": 0,
        "terminal": "none",
        "game_eligibility": ["learn", "ranked"],
        "family": "Metal Chelation / Host Defense",
        "hint_role": "Conformationally restricted divalent-metal and heme-binding pentapeptide repeated 12 times in serum HRG.",
        "hint_detail": "Pairs an N-terminal glycine with three imidazole rings and a turn-inducing cyclic proline.",
        "context_requirement": "Tandemly repeated in the central histidine-proline-rich (HPR) domain of human plasma Histidine-Rich Glycoprotein.",
        "exemplar_source": "UniProt:P04196 (Human Histidine-Rich Glycoprotein HRG 383-387)",
    },
]


def _build_chemle_candidate_extractors(records: Sequence[Dict[str, Any]]) -> List[Tuple[str, Callable[[Dict[str, Any]], str]]]:
    ha_vals = sorted(r["HeavyAtoms"] for r in records)
    q1, q2 = ha_vals[len(ha_vals) // 3], ha_vals[(2 * len(ha_vals)) // 3]

    def ha_tercile(r: Dict[str, Any]) -> str:
        if r["HeavyAtoms"] <= q1:
            return f"Size Band: Small (≤{q1} heavy atoms)"
        if r["HeavyAtoms"] <= q2:
            return f"Size Band: Medium ({q1+1}–{q2} heavy atoms)"
        return f"Size Band: Large (≥{q2+1} heavy atoms)"

    def ring_band(r: Dict[str, Any]) -> str:
        return "Topology: 0–1 Rings" if r["Rings"] <= 1 else "Topology: Polycyclic (≥2 Rings)"

    def polarity_band(r: Dict[str, Any]) -> str:
        return "Polarity: Low-Moderate TPSA (<60 Å²)" if r["TPSA"] < 60.0 else "Polarity: High TPSA (≥60 Å²)"

    return [
        ("heavy_atom_tercile", ha_tercile),
        ("ring_topology_band", ring_band),
        ("polarity_band", polarity_band),
    ]


def _build_seqle_candidate_extractors() -> List[Tuple[str, Callable[[Dict[str, Any]], str]]]:
    def charge_band(r: Dict[str, Any]) -> str:
        q = r["net_charge_ph74"]
        if q < -0.5:
            return "Charge Band (pH 7.4): Acidic (Q < -0.5)"
        if q > 0.5:
            return "Charge Band (pH 7.4): Basic (Q > +0.5)"
        return "Charge Band (pH 7.4): Near-Neutral (-0.5 ≤ Q ≤ +0.5)"

    def gravy_band(r: Dict[str, Any]) -> str:
        return (
            "Hydropathy Band: Hydrophilic (GRAVY < -1.0)"
            if r["gravy"] < -1.0
            else "Hydropathy Band: Moderate / Hydrophobic (GRAVY ≥ -1.0)"
        )

    def arom_band(r: Dict[str, Any]) -> str:
        return (
            "Aromatic Content: Contains Aromatic Residue(s)"
            if r["aromatic_count"] > 0
            else "Aromatic Content: No Aromatic Residues"
        )

    return [
        ("charge_band", charge_band),
        ("gravy_band", gravy_band),
        ("arom_band", arom_band),
    ]


def compute_canonical_dataset_hash(
    chemle_records: Sequence[Dict[str, Any]],
    seqle_records: Sequence[Dict[str, Any]],
    feedback_spec: Dict[str, Any],
) -> str:
    """Computes an order-independent SHA-256 hash over all gameplay-relevant data (P1-3)."""
    canonical_chemle = [
        {
            "id": r["id"],
            "inchikey": r["inchikey"],
            "canonical_smiles": r["canonical_smiles"],
            "Category": r["Category"],
            "super_class": r["super_class"],
            "MW": r["MW"],
            "cLogP": r["cLogP"],
            "TPSA": r["TPSA"],
            "HBD": r["HBD"],
            "HBA": r["HBA"],
            "Rings": r["Rings"],
            "AromaticRings": r["AromaticRings"],
            "HeavyAtoms": r["HeavyAtoms"],
            "functional_groups": r["functional_groups"],
            "game_eligibility": r.get("game_eligibility", ["learn", "ranked"]),
            "guess_eligibility": r.get("guess_eligibility", r.get("game_eligibility", ["learn", "ranked"])),
            "clue_ladder": r["clue_ladder"],
        }
        for r in sorted(chemle_records, key=lambda x: x["id"])
    ]
    canonical_seqle = [
        {
            "seq": r["seq"],
            "kind": r["kind"],
            "consensus": r["consensus"],
            "parent_sequence": r["parent_sequence"],
            "window_offset": r["window_offset"],
            "terminal": r["terminal"],
            "game_eligibility": r["game_eligibility"],
            "family": r["family"],
            "functional_rule": r.get("functional_rule", ""),
            "provenance_status": r.get("provenance_status", ""),
            "net_charge_ph74": r["net_charge_ph74"],
            "gravy": r["gravy"],
            "isoelectric_point": r["isoelectric_point"],
            "residue_groups": r["residue_groups"],
            "clue_ladder": r["clue_ladder"],
        }
        for r in sorted(seqle_records, key=lambda x: x["seq"])
    ]
    payload = {
        "schema_version": SCHEMA_VERSION,
        "fingerprint_spec": FINGERPRINT_SPEC,
        "standardization_spec": STANDARDIZATION_SPEC,
        "peptide_pka_spec": PEPTIDE_PKA_SPEC,
        "feedback_spec": feedback_spec,
        "chemle": canonical_chemle,
        "seqle": canonical_seqle,
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:12]


def compute_chemle_rules_hash(chemle_records, similarity_matrix, feedback_spec):
    """Identity of Chemle gameplay inputs; increment engine revision for algorithm changes.

    This is independent of release and roster schedule identities. It includes the
    full playable catalogue because any of its compounds can be a submitted guess.
    """
    fields = ("id", "inchikey", "canonical_smiles", "MW", "cLogP", "TPSA",
              "HBD", "HBA", "Rings", "AromaticRings", "HeavyAtoms",
              "functional_groups", "Category", "super_class", "clue_ladder",
              "game_eligibility", "guess_eligibility")
    payload = {"engine_revision": 1, "fingerprint_spec": FINGERPRINT_SPEC,
               "standardization_spec": STANDARDIZATION_SPEC,
               "feedback_spec": feedback_spec,
               "records": [{k: r.get(k) for k in fields}
                           for r in sorted(chemle_records, key=lambda x: x["id"])],
               "similarity": similarity_matrix}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()[:12]


LEGACY_SCHEDULE_HASHES = ["9257bd168fa6", "8bb96ac156cc", "919ad8afb190"]
V021_CHEMLE_ROSTER_DIGEST = "87689dabf71408a8"


def compute_chemle_schedule_hash(chemle_records: Sequence[Dict[str, Any]]) -> str:
    """Computes per-game daily schedule hash strictly from Ranked-eligible Chemle answer IDs."""
    ranked_ids = sorted(
        r["id"]
        for r in chemle_records
        if "ranked" in r.get("game_eligibility", ["learn", "ranked"])
    )
    roster_digest = hashlib.sha256(
        json.dumps(ranked_ids, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:16]
    # Preserve the v0.2.1 Ranked Chemle daily & classroom schedule when the answer roster is unchanged.
    if roster_digest == V021_CHEMLE_ROSTER_DIGEST:
        return "9257bd168fa6"
    payload = {"game": "chemle", "chooser": "fnv1a_v1", "ranked_ids": ranked_ids}
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:12]


def compute_seqle_schedule_hash(seqle_records: Sequence[Dict[str, Any]]) -> str:
    """Computes per-game daily schedule hash strictly from Ranked-eligible Seqle answer sequences."""
    ranked_seqs = sorted(
        r["seq"]
        for r in seqle_records
        if "ranked" in r.get("game_eligibility", ["learn"])
    )
    payload = {"game": "seqle", "chooser": "fnv1a_v1", "ranked_seqs": ranked_seqs}
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:12]


def compute_schedule_version_hash(
    chemle_records: Sequence[Dict[str, Any]],
    seqle_records: Sequence[Dict[str, Any]],
) -> str:
    """Combines per-game Ranked schedule hashes for backwards-compatible callers."""
    chem_h = compute_chemle_schedule_hash(chemle_records)
    seq_h = compute_seqle_schedule_hash(seqle_records)
    blob = f"chemle:{chem_h}|seqle:{seq_h}".encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:12]


def build_all(
    feedback_spec: Dict[str, Any],
    chemle_raw_compounds: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Compiles, validates, and audits the complete Chemle + Seqle puzzle database."""
    # Sort raw catalogs canonically by ID / sequence first so output order is invariant (P1-3)
    raw_chem_source = load_chemle_manifest() if chemle_raw_compounds is None else chemle_raw_compounds
    if raw_chem_source and all("input_isomeric_smiles" in r for r in raw_chem_source):
        validate_chemle_manifest(raw_chem_source)
    sorted_raw_chem = sorted(raw_chem_source, key=lambda x: x["id"])
    sorted_raw_seq = sorted(SEQLE_RAW_MOTIFS, key=lambda x: x["seq"])

    chemle_records: List[Dict[str, Any]] = []
    chemle_fps: List[Any] = []
    seen_ids: set[str] = set()
    seen_inchikeys: set[str] = set()
    seen_smiles: set[str] = set()

    print(
        f"[BioWordle Build v{SCHEMA_VERSION}] Auditing & compiling {len(sorted_raw_chem)} Chemle molecules "
        f"with RDKit {Chem.rdBase.rdkitVersion}..."
    )
    for item in sorted_raw_chem:
        cid = item["id"]
        if cid in seen_ids:
            raise ValueError(f"Duplicate Chemle ID: {cid}")
        seen_ids.add(cid)

        game_elig = list(item.get("game_eligibility", ["learn", "ranked"]))
        guess_elig = list(item.get("guess_eligibility", game_elig))
        if not set(game_elig).issubset({"learn", "ranked"}) or not set(guess_elig).issubset({"learn", "ranked"}):
            raise ValueError(f"Invalid eligibility for {cid}: game={game_elig}, guess={guess_elig}")
        if not set(game_elig).issubset(set(guess_elig)):
            raise ValueError(f"Eligible answer profile(s) missing from guess_eligibility for {cid}: game={game_elig}, guess={guess_elig}")
        if not game_elig and not guess_elig:
            # Draft/rejected non-playable manifest entries are excluded from compiled gameplay pools
            continue

        raw_smiles = str(item.get("input_isomeric_smiles") or item.get("smiles") or "")
        mol, audit = standardize_and_audit(raw_smiles)
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
        svg_sha = hashlib.sha256(svg_full.encode("utf-8")).hexdigest()[:12]

        fp = compute_chiral_morgan_fp(mol)
        chemle_fps.append(fp)

        record: Dict[str, Any] = {
            "id": cid,
            "name": item["name"],
            "aliases": list(item.get("aliases", [])),
            "Category": item["Category"],
            "super_class": SUPER_CLASS_MAP[item["Category"]],
            "primary_class": item.get("primary_class") or SUPER_CLASS_MAP[item["Category"]],
            "role_tags": list(item.get("role_tags", [])),
            "stereo_status": item.get("stereo_status", "unspecified_legacy"),
            "protonation_policy": item.get("protonation_policy", "neutral_canonical"),
            "structure_refs": dict(item.get("structure_refs", {})),
            "evidence_refs": list(item.get("evidence_refs", [])),
            "review_status": item.get("review_status", "baseline_legacy"),
            "review_notes": item.get("review_notes", ""),
            "game_eligibility": game_elig,
            "guess_eligibility": guess_elig,
            **audit,
            **desc,
            "functional_groups": fg_list,
            "hint_role": item["hint_role"],
            "hint_origin": item["hint_origin"],
            "structural_note": item["structural_note"],
            "scaffold_smiles": scaffold_smiles,
            "scaffold_mode": scaffold_mode,
            "has_scaffold": scaffold_mode in ("bemis_murcko", "ring_core", "carbon_skeleton"),
            "svg_sha256": svg_sha,
            "svg_scaffold": svg_scaffold,
            "svg_full": svg_full,
        }
        # P2-3 Hint Leakage Audit
        audit_hint_text_leakage(record, ("hint_role", "hint_origin"), "id")
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

    # P0-1: Separate profile answer pools & legal guess pools; fit extractors on each profile's answer universe
    chem_pools = {
        prof: {
            "answers": [r for r in chemle_records if prof in r["game_eligibility"]],
            "guesses": [r for r in chemle_records if prof in r["guess_eligibility"]],
        }
        for prof in ("ranked", "learn")
    }
    for prof in ("ranked", "learn"):
        ans_ids = {r["id"] for r in chem_pools[prof]["answers"]}
        guess_ids = {r["id"] for r in chem_pools[prof]["guesses"]}
        if not ans_ids or not ans_ids.issubset(guess_ids):
            raise ValueError(f"Invalid Chemle {prof} pools: answers={len(ans_ids)} guesses={len(guess_ids)}")

    r_min_k = int(feedback_spec["profiles"]["ranked"]["min_k_level0"])
    l_min_k = int(feedback_spec["profiles"]["learn"]["min_k_level0"])

    r_chem_answers = chem_pools["ranked"]["answers"]
    l_chem_answers = chem_pools["learn"]["answers"]

    r_chem_extractors = _build_chemle_candidate_extractors(r_chem_answers)
    r_chem_feats, r_chem_lvl0, r_chem_kmin, r_chem_groups = select_k_anonymous_level0_features(
        r_chem_answers, r_chem_extractors, min_k=r_min_k, max_features=1
    )
    l_chem_extractors = _build_chemle_candidate_extractors(l_chem_answers)
    l_chem_feats, l_chem_lvl0, l_chem_kmin, l_chem_groups = select_k_anonymous_level0_features(
        l_chem_answers, l_chem_extractors, min_k=l_min_k, max_features=2
    )
    r_lvl0_by_id = {rec["id"]: r_chem_lvl0[idx] for idx, rec in enumerate(r_chem_answers)}
    l_lvl0_by_id = {rec["id"]: l_chem_lvl0[idx] for idx, rec in enumerate(l_chem_answers)}

    learn_pack_audit: Dict[str, Any] = {}
    for cat in sorted({r["Category"] for r in l_chem_answers}):
        cat_recs = [r for r in l_chem_answers if r["Category"] == cat]
        pack_counts: Dict[str, int] = {}
        for r in cat_recs:
            txt = l_lvl0_by_id[r["id"]]
            pack_counts[txt] = pack_counts.get(txt, 0) + 1
        sizes = sorted(pack_counts.values())
        k_min_pack = min(sizes) if sizes else 0
        learn_pack_audit[cat] = {
            "count": len(cat_recs),
            "k_min": k_min_pack,
            "group_sizes": sizes,
            "exploratory": k_min_pack < l_min_k,
        }

    print(
        f"[P0-1 Chemle k-Anonymity] Ranked Level-0 ({r_chem_feats}, N={len(r_chem_answers)}): k_min={r_chem_kmin}, groups={r_chem_groups} | "
        f"Learn Level-0 ({l_chem_feats}, N={len(l_chem_answers)}): k_min={l_chem_kmin}, groups={l_chem_groups}"
    )

    ranked_super_unlock = int(feedback_spec["profiles"]["ranked"].get("super_class_unlock_after_guesses", 4))
    learn_super_unlock = int(feedback_spec["profiles"]["learn"].get("super_class_unlock_after_guesses", 1))
    for rec in chemle_records:
        cid = rec["id"]
        ladders: Dict[str, List[Dict[str, Any]]] = {}
        if cid in r_lvl0_by_id:
            ladders["ranked"] = [
                {"level": 0, "after_guesses": 0, "title": "Clue 1 · Molecular Size Band (k-Anonymous)", "text": r_lvl0_by_id[cid]},
                {"level": 1, "after_guesses": ranked_super_unlock, "title": "Clue 2 · Chemical Super-Class", "text": f"Broad Super-Class: {rec['super_class']}"},
                {"level": 2, "after_guesses": 4, "title": "Clue 3 · Biochemical Role & Mechanism", "text": rec["hint_role"]},
                {"level": 3, "after_guesses": 5, "title": "Clue 4 · Biosynthesis & Pathway Origin", "text": rec["hint_origin"]},
            ]
        if cid in l_lvl0_by_id:
            ladders["learn"] = [
                {"level": 0, "after_guesses": 0, "title": "Clue 1 · Chemical Space Band (k-Anonymous)", "text": l_lvl0_by_id[cid]},
                {"level": 1, "after_guesses": learn_super_unlock, "title": "Clue 2 · Chemical Super-Class", "text": f"Broad Super-Class: {rec['super_class']}"},
                {"level": 2, "after_guesses": 2, "title": "Clue 3 · Biochemical Role & Mechanism", "text": rec["hint_role"]},
                {"level": 3, "after_guesses": 3, "title": "Clue 4 · Biosynthesis & Pathway Origin", "text": rec["hint_origin"]},
            ]
        rec["clue_ladder"] = ladders

    # P0-2: Run clue-conditioned greedy max-entropy solver over each profile's answers & legal probes
    solver_reports: Dict[str, Any] = {}
    for prof_id, prof_spec in feedback_spec["profiles"].items():
        ans_ids = [r["id"] for r in chem_pools[prof_id]["answers"]]
        probe_ids = [r["id"] for r in chem_pools[prof_id]["guesses"]]
        stats = run_chemle_greedy_solver(
            chemle_records,
            similarity_matrix,
            prof_spec,
            profile_id=prof_id,
            condition_on_clues=True,
            answer_ids=ans_ids,
            probe_ids=probe_ids,
        )
        assert_solver_gate(stats, prof_spec["solver_gate"], prof_id)
        solver_reports[prof_id] = {
            k: v
            for k, v in stats.items()
            if k not in ("probe_metrics", "per_target_guesses", "per_target_results")
        }
        print(
            f"[P0-2 Solver Gate PASS · {prof_id}] targets={stats['eligible_count']} probes={stats['probe_count']} "
            f"mean={stats['mean']} median={stats['median']} "
            f"max={stats['max']} <=2:{stats['pct_solved_le2']}% <=6:{stats['pct_solved_le6']}% "
            f"min_exp_rem={stats['min_opening_expected_remaining']} dist={stats['distribution']}"
        )
        for rec in chemle_records:
            if rec["id"] in stats["probe_metrics"]:
                rec.setdefault("opening_probe_quality", {})[prof_id] = stats["probe_metrics"][rec["id"]]

    # Attach post-game nearest-neighbor info (explicitly labeled for post-game only)
    for rec in chemle_records:
        cid = rec["id"]
        other_sims = sorted(
            (similarity_matrix[cid][other["id"]], other["name"])
            for other in chemle_records
            if other["id"] != cid
        )
        max_sim, nearest_name = other_sims[-1] if other_sims else (0, "None")
        r_qual = rec.get("opening_probe_quality", {}).get("ranked", {})
        l_qual = rec.get("opening_probe_quality", {}).get("learn", {})
        rec["postgame_insights"] = {
            "nearest_neighbor": nearest_name,
            "nearest_neighbor_tanimoto": max_sim,
            "opening_ig_ranked_bits": r_qual.get("information_gain_bits"),
            "opening_ig_learn_bits": l_qual.get("information_gain_bits"),
            "expected_remaining_ranked": r_qual.get("expected_remaining_candidates"),
        }

    print(
        f"[BioWordle Build v{SCHEMA_VERSION}] Auditing & compiling {len(sorted_raw_seq)} Seqle 5-AA motifs..."
    )
    seqle_records: List[Dict[str, Any]] = []
    seen_seqs: set[str] = set()
    for item in sorted_raw_seq:
        seq = item["seq"].strip().upper()
        if seq in seen_seqs:
            raise ValueError(f"Duplicate Seqle sequence: {seq}")
        seen_seqs.add(seq)

        # P0-4 Schema & Provenance Validator
        validate_seqle_motif_record(item, AA_PROPERTIES)

        pep_mol = Chem.MolFromSequence(seq)
        if pep_mol is None:
            raise ValueError(f"RDKit MolFromSequence failed for: {seq}")

        pep_desc = compute_descriptors_and_insights(pep_mol)
        biochem = compute_peptide_biochem(seq, AA_PROPERTIES)
        svg_peptide = render_responsive_svg(pep_mol, width=380, height=210, dark_mode=False)
        svg_sha = hashlib.sha256(svg_peptide.encode("utf-8")).hexdigest()[:12]

        ref_prov = SEQLE_REFERENCE_REGISTRY.get(seq) or {}
        seq_rec: Dict[str, Any] = {
            "seq": seq,
            "title": item["title"],
            "kind": item["kind"],
            "consensus": item["consensus"],
            "pattern": item["consensus"],
            "functional_rule": ref_prov.get("functional_rule", item["consensus"]),
            "parent_sequence": item["parent_sequence"],
            "window_offset": item["window_offset"],
            "terminal": item["terminal"],
            "game_eligibility": item["game_eligibility"],
            "family": item["family"],
            "hint_role": item["hint_role"],
            "hint_detail": item["hint_detail"],
            "context_requirement": item["context_requirement"],
            "exemplar_source": item["exemplar_source"],
            "provenance_status": ref_prov.get("status", "unknown"),
            "provenance_verified_date": ref_prov.get("provenance_verified_date", PROVENANCE_VERIFIED_DATE),
            "uniprot_release": ref_prov.get("uniprot_release", UNIPROT_VERIFIED_RELEASE),
            "citations": ref_prov.get("citations", []),
            "reference_provenance": ref_prov,
            "three_letter": "-".join(AA_PROPERTIES[aa]["three"] for aa in seq),
            "formula": pep_desc["formula"],
            "Formula": pep_desc["formula"],
            "MW": pep_desc["MW"],
            "cLogP": pep_desc["cLogP"],
            "TPSA": pep_desc["TPSA"],
            "HBD": pep_desc["HBD"],
            "HBA": pep_desc["HBA"],
            "HeavyAtoms": pep_desc["HeavyAtoms"],
            **biochem,
            "svg_sha256": svg_sha,
            "svg_peptide": svg_peptide,
            "svg_full": svg_peptide,
        }
        validate_seqle_motif_record(seq_rec, AA_PROPERTIES, enforce_reference_registry=True)
        audit_hint_text_leakage(seq_rec, ("hint_role", "hint_detail"), "seq")
        seqle_records.append(seq_rec)

    validate_seqle_registry_completeness(seqle_records)

    # P0-1: Select k-anonymous Level-0 opening clues for Seqle
    seq_extractors = _build_seqle_candidate_extractors()
    ranked_seq_pool = [r for r in seqle_records if "ranked" in r["game_eligibility"]]
    r_seq_feats, r_seq_lvl0_map, r_seq_kmin, r_seq_groups = select_k_anonymous_level0_features(
        ranked_seq_pool, seq_extractors, min_k=r_min_k, max_features=1
    )
    l_seq_feats, l_seq_lvl0, l_seq_kmin, l_seq_groups = select_k_anonymous_level0_features(
        seqle_records, seq_extractors, min_k=l_min_k, max_features=2
    )
    print(
        f"[P0-1 Seqle k-Anonymity] Ranked Level-0 ({r_seq_feats}): k_min={r_seq_kmin}, groups={r_seq_groups} | "
        f"Learn Level-0 ({l_seq_feats}): k_min={l_seq_kmin}, groups={l_seq_groups}"
    )

    # Map Ranked Level-0 extractor onto all records
    r_seq_fn = dict(seq_extractors)[r_seq_feats[0]]
    for idx, rec in enumerate(seqle_records):
        kind_badge = (
            f"5-mer window of {rec['parent_sequence']} (offset {rec['window_offset']})"
            if rec["kind"] == "fragment" and rec["parent_sequence"] != rec["seq"]
            else (
                f"Consensus teaching example ({rec['functional_rule']})"
                if rec.get("provenance_status") == "consensus_exemplar"
                else ("Pedagogical coordination abstraction (Learn mode)" if rec["kind"] == "abstraction" else "Canonical 5-residue motif")
            )
        )
        rec["clue_ladder"] = {
            "ranked": [
                {"level": 0, "after_guesses": 0, "title": "Clue 1 · Peptide Charge Band (k-Anonymous)", "text": r_seq_fn(rec)},
                {"level": 1, "after_guesses": 2, "title": "Clue 2 · Motif Architecture & Terminal Context", "text": f"Kind: {rec['kind'].upper()} · Terminal: {rec['terminal']}"},
                {"level": 2, "after_guesses": 4, "title": "Clue 3 · Biological Function", "text": rec["hint_role"]},
                {"level": 3, "after_guesses": 5, "title": "Clue 4 · Biochemical Mechanism", "text": rec["hint_detail"]},
            ],
            "learn": [
                {"level": 0, "after_guesses": 0, "title": "Clue 1 · Biophysics Band (k-Anonymous)", "text": l_seq_lvl0[idx]},
                {"level": 1, "after_guesses": 1, "title": "Clue 2 · Functional Family & Window Context", "text": f"{rec['family']} · {kind_badge}"},
                {"level": 2, "after_guesses": 2, "title": "Clue 3 · Biological Function", "text": rec["hint_role"]},
                {"level": 3, "after_guesses": 3, "title": "Clue 4 · Biochemical Mechanism", "text": rec["hint_detail"]},
            ],
        }

    seq_info = compute_seqle_information_metrics(seqle_records, AA_PROPERTIES)
    for rec in seqle_records:
        rec["postgame_insights"] = seq_info[rec["seq"]]

    seqle_solver_reports: Dict[str, Any] = {}
    for prof_id, prof_spec in feedback_spec["profiles"].items():
        seq_stats = run_seqle_greedy_solver(
            seqle_records,
            AA_PROPERTIES,
            prof_spec,
            profile_id=prof_id,
            condition_on_clues=True,
        )
        assert_seqle_solver_gate(seq_stats, prof_spec["seqle_solver_gate"], prof_id)
        seqle_solver_reports[prof_id] = {
            k: v
            for k, v in seq_stats.items()
            if k not in ("per_target_guesses", "per_target_results")
        }

    dataset_hash = compute_canonical_dataset_hash(chemle_records, seqle_records, feedback_spec)
    chemle_schedule_hash = compute_chemle_schedule_hash(chemle_records)
    seqle_schedule_hash = compute_seqle_schedule_hash(seqle_records)
    schedule_hash = compute_schedule_version_hash(chemle_records, seqle_records)

    return {
        "metadata": {
            "schema_version": SCHEMA_VERSION,
            "release_version": "0.3.1",
            "chemle_rules_hash": compute_chemle_rules_hash(chemle_records, similarity_matrix, feedback_spec),
            "generator": "BioWordle Independent RDKit Build Compiler (chem_engine.py)",
            "rdkit_version": Chem.rdBase.rdkitVersion,
            "dataset_version_hash": dataset_hash,
            "schedule_version_hash": schedule_hash,
            "chemle_schedule_hash": chemle_schedule_hash,
            "seqle_schedule_hash": seqle_schedule_hash,
            "legacy_schedule_hashes": LEGACY_SCHEDULE_HASHES,
            "retired_seqle_answers": ["KKMPQ"],
            "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
            "uniprot_release": UNIPROT_VERIFIED_RELEASE,
            "daily_seed_policy": "FNV-1a over '{game}:{profile}:{utcDate}:{gameScheduleHash}' indexing canonically sorted Ranked answer pool",
            "chemle_count": len(chemle_records),
            "chemle_ranked_count": len(chem_pools["ranked"]["answers"]),
            "chemle_ranked_guess_count": len(chem_pools["ranked"]["guesses"]),
            "chemle_learn_count": len(chem_pools["learn"]["answers"]),
            "chemle_learn_guess_count": len(chem_pools["learn"]["guesses"]),
            "seqle_count": len(seqle_records),
            "seqle_ranked_count": len(ranked_seq_pool),
            "fingerprint_spec": FINGERPRINT_SPEC,
            "standardization_spec": STANDARDIZATION_SPEC,
            "peptide_pka_spec": PEPTIDE_PKA_SPEC,
            "k_anonymity_audit": {
                "chemle_ranked": {"features": r_chem_feats, "k_min": r_chem_kmin, "group_sizes": r_chem_groups},
                "chemle_learn": {"features": l_chem_feats, "k_min": l_chem_kmin, "group_sizes": l_chem_groups},
                "chemle_learn_packs": learn_pack_audit,
                "seqle_ranked": {"features": r_seq_feats, "k_min": r_seq_kmin, "group_sizes": r_seq_groups},
                "seqle_learn": {"features": l_seq_feats, "k_min": l_seq_kmin, "group_sizes": l_seq_groups},
            },
            "solver_gate_audit": solver_reports,
            "seqle_solver_gate_audit": seqle_solver_reports,
        },
        "feedback_spec": feedback_spec,
        "aa_properties": AA_PROPERTIES,
        "chemle": chemle_records,
        "chemle_similarity": similarity_matrix,
        "seqle": seqle_records,
    }


def write_golden_feedback_vectors(payload: Dict[str, Any], base_dir: str) -> str:
    """Writes `tests/golden_feedback.json` for cross-language Python <-> JS verification (P1-1)."""
    tests_dir = os.path.join(base_dir, "tests")
    os.makedirs(tests_dir, exist_ok=True)

    chemle = payload["chemle"]
    sim = payload["chemle_similarity"]
    seqle = payload["seqle"]
    aa_props = payload["aa_properties"]
    f_spec = payload["feedback_spec"]

    chemle_vectors: Dict[str, Dict[str, str]] = {"ranked": {}, "learn": {}}
    for prof_id, prof_spec in f_spec["profiles"].items():
        for g in chemle:
            for t in chemle:
                key = f"{g['id']}->{t['id']}"
                res = evaluate_chemle_feedback_from_spec(
                    g, t, sim[g["id"]][t["id"]], prof_spec
                )
                chemle_vectors[prof_id][key] = res["signature"]

    seqle_vectors: Dict[str, Dict[str, str]] = {"ranked": {}, "learn": {}}
    for prof_id, prof_spec in f_spec["profiles"].items():
        for g in seqle:
            for t in seqle:
                key = f"{g['seq']}->{t['seq']}"
                res = evaluate_seqle_feedback_from_spec(
                    g, t, g["seq"], t["seq"], aa_props, prof_spec
                )
                seqle_vectors[prof_id][key] = res["signature"]

    golden_path = os.path.join(tests_dir, "golden_feedback.json")
    tmp_golden = golden_path + ".tmp"
    with open(tmp_golden, "w", encoding="utf-8", newline="\n") as f:
        json.dump(
            {
                "dataset_version_hash": payload["metadata"]["dataset_version_hash"],
                "chemle_signatures": chemle_vectors,
                "seqle_signatures": seqle_vectors,
            },
            f,
            indent=2,
            sort_keys=True,
        )
    os.replace(tmp_golden, golden_path)
    return golden_path


def main() -> None:
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, "data")
    os.makedirs(data_dir, exist_ok=True)

    spec_path = os.path.join(data_dir, "feedback_spec.json")
    with open(spec_path, "r", encoding="utf-8") as f:
        feedback_spec = json.load(f)

    payload = build_all(feedback_spec)

    json_path = os.path.join(data_dir, "puzzles.json")
    tmp_json = json_path + ".tmp"
    with open(tmp_json, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=2)
    os.replace(tmp_json, json_path)

    js_path = os.path.join(data_dir, "puzzles_data.js")
    tmp_js = js_path + ".tmp"
    with open(tmp_js, "w", encoding="utf-8", newline="\n") as f:
        f.write("// Auto-generated by build_puzzles.py using RDKit — do not edit manually.\n")
        f.write("window.BIOWORDLE_DATA = ")
        json.dump(payload, f, separators=(",", ":"))
        f.write(";\n")
    os.replace(tmp_js, js_path)

    golden_path = write_golden_feedback_vectors(payload, base_dir)

    print(
        f"[BioWordle Build] Done! Dataset hash: {payload['metadata']['dataset_version_hash']}\n"
        f"  - {json_path}\n"
        f"  - {js_path}\n"
        f"  - {golden_path}"
    )

    index_path = os.path.join(base_dir, "index.html")
    css_path = os.path.join(base_dir, "styles.css")
    app_js_path = os.path.join(base_dir, "app.js")
    if all(os.path.isfile(p) for p in (index_path, css_path, app_js_path)):
        with open(index_path, "r", encoding="utf-8") as f:
            html_src = f.read().replace("\r\n", "\n")
        with open(css_path, "r", encoding="utf-8") as f:
            css_src = f.read().replace("\r\n", "\n")
        with open(js_path, "r", encoding="utf-8") as f:
            data_js_src = f.read().replace("\r\n", "\n")
        with open(app_js_path, "r", encoding="utf-8") as f:
            app_js_src = f.read().replace("\r\n", "\n")

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
        with open(tmp_standalone, "w", encoding="utf-8", newline="\n") as f:
            f.write(bundled)
        os.replace(tmp_standalone, standalone_path)
        print(f"  - {standalone_path} (Single-file standalone bundle)")

        if len(sys.argv) > 1 and sys.argv[1] == "--artifact-output" and len(sys.argv) > 2:
            art_path = sys.argv[2]
            tmp_art = art_path + ".tmp"
            with open(tmp_art, "w", encoding="utf-8", newline="\n") as f:
                f.write(bundled)
            os.replace(tmp_art, art_path)
            print(f"  - {art_path} (Artifact preview updated)")


if __name__ == "__main__":
    main()
