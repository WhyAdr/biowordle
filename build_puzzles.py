# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "rdkit",
#     "pandas",
#     "numpy",
# ]
# ///
"""Build-time RDKit Puzzle Compiler for BioWordle (Seqle + Chemle).

Uses the verified `rdkit-general` and `cheminformatics-core` skill modules to:
1. Standardize and validate 36 iconic metabolites, hydrocarbon-degradation
   intermediates, natural products, and drugs for CHEMLE.
2. Compute canonical 2D descriptors (MW, cLogP, TPSA, HBD, HBA, RotB, Fsp3,
   HeavyAtoms, Rings, AromaticRings, AliphaticRings, QED, SA_Score, CNS_MPO),
   Lipinski Rule-of-5 / PAINS / Brenk alerts, SMARTS functional groups, and
   the full N x N Morgan fingerprint Tanimoto similarity matrix.
3. Render responsive ACS-1996 2D SVGs for every molecule and its Bemis-Murcko
   scaffold core.
4. Build and validate 30 real 5-amino-acid functional protein motifs for SEQLE
   using `Chem.MolFromSequence`, computing both peptide-level RDKit descriptors
   and residue-level physicochemical properties + 2D pentapeptide SVGs.
5. Emit `data/puzzles.json` and `data/puzzles_data.js` atomically.
"""

from __future__ import annotations

import json
import os
import re
import sys
from typing import Any, Dict, List, Optional

# ---------------------------------------------------------------------------
# 1. Discover and attach the Science Skills (`rdkit_general` & `cheminformatics_core`)
# ---------------------------------------------------------------------------
_DEFAULT_SKILL_ROOT = os.path.expanduser(
    r"~\.gemini\config\plugins\science\skills"
)
_RDKIT_GENERAL_SCRIPTS = os.path.join(
    _DEFAULT_SKILL_ROOT, "rdkit_general", "scripts"
)
if os.path.isdir(_RDKIT_GENERAL_SCRIPTS) and _RDKIT_GENERAL_SCRIPTS not in sys.path:
    sys.path.insert(0, _RDKIT_GENERAL_SCRIPTS)

try:
    from skill_paths import add_skill_scripts

    add_skill_scripts("rdkit_general")
    add_skill_scripts("cheminformatics_core")
except Exception:
    _core_scripts = os.path.join(
        _DEFAULT_SKILL_ROOT, "cheminformatics_core", "scripts"
    )
    if os.path.isdir(_core_scripts) and _core_scripts not in sys.path:
        sys.path.insert(0, _core_scripts)

from rdkit import Chem, DataStructs, RDConfig  # noqa: E402
from rdkit.Chem import (  # noqa: E402
    Crippen,
    Descriptors,
    Lipinski,
    QED,
    rdDepictor,
    rdFingerprintGenerator,
    rdMolDescriptors,
)
from rdkit.Chem.Draw import rdMolDraw2D  # noqa: E402
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams  # noqa: E402
from rdkit.Chem.MolStandardize import rdMolStandardize  # noqa: E402
from rdkit.Chem.Scaffolds import MurckoScaffold  # noqa: E402

try:
    from compat import (  # noqa: E402
        get_morgan_fingerprint,
        safe_substruct_match,
        standardize_molecule,
    )
    from descriptors import (  # noqa: E402
        calculate_cns_mpo,
        compute_standard_descriptors,
    )
    from filter_catalogs import (  # noqa: E402
        calculate_rule_of_5,
        screen_structural_alerts,
    )
except ImportError:
    # Built-in portable RDKit stand-ins when running outside the Antigravity Science skill tree
    _contrib_sa = os.path.join(RDConfig.RDContribDir, "SA_Score")
    if os.path.isdir(_contrib_sa) and _contrib_sa not in sys.path:
        sys.path.insert(0, _contrib_sa)
    try:
        import sascorer as _sascorer  # type: ignore
    except Exception:
        _sascorer = None

    def standardize_molecule(mol: Chem.Mol, uncharge: bool = False) -> Chem.Mol:
        if mol is None:
            return mol
        try:
            clean = rdMolStandardize.Cleanup(mol)
            parent = rdMolStandardize.FragmentParent(clean)
            if uncharge:
                parent = rdMolStandardize.Uncharger().uncharge(parent)
            te = rdMolStandardize.TautomerEnumerator()
            return te.Canonicalize(parent)
        except Exception:
            return mol

    def get_morgan_fingerprint(mol: Chem.Mol, radius: int = 2, n_bits: int = 2048) -> Any:
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
        return gen.GetFingerprint(mol)

    def safe_substruct_match(mol: Chem.Mol, pattern: Chem.Mol) -> bool:
        if mol is None or pattern is None:
            return False
        try:
            return bool(mol.HasSubstructMatch(pattern))
        except Exception:
            return False

    def calculate_cns_mpo(mol: Chem.Mol) -> float:
        mw = float(Descriptors.MolWt(mol))
        clogp = float(Crippen.MolLogP(mol))
        tpsa = float(rdMolDescriptors.CalcTPSA(mol))
        hbd = int(Lipinski.NumHDonors(mol))
        basic_amine = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=[O,S]);!$(NS=O);!$(Nc);!$(N#*)]")
        has_basic_amine = mol.HasSubstructMatch(basic_amine) if basic_amine else False
        pka = 9.5 if has_basic_amine else 5.0
        clogd = clogp - (pka - 7.4) if pka > 7.4 else clogp

        def _rescale(x: float, min_x: float, max_x: float) -> float:
            return min(1.0, max(0.0, (max_x - x) / (max_x - min_x)))

        s_clogp = _rescale(clogp, 3.0, 5.0)
        s_clogd = _rescale(clogd, 2.0, 4.0)
        s_mw = _rescale(mw, 360.0, 500.0)
        if 40.0 <= tpsa <= 90.0:
            s_tpsa = 1.0
        elif 20.0 < tpsa < 40.0:
            s_tpsa = (tpsa - 20.0) / 20.0
        elif 90.0 < tpsa < 120.0:
            s_tpsa = _rescale(tpsa, 90.0, 120.0)
        else:
            s_tpsa = 0.0
        s_hbd = _rescale(float(hbd), 0.5, 3.5)
        s_pka = _rescale(pka, 8.0, 10.0)
        return round(float(s_clogp + s_clogd + s_mw + s_tpsa + s_hbd + s_pka), 2)

    def compute_standard_descriptors(
        mol: Chem.Mol, style: str = "short", include_sa_score: bool = True
    ) -> Dict[str, Any]:
        del style
        sa_val = 2.5
        if include_sa_score and _sascorer is not None:
            try:
                sa_val = float(_sascorer.calculateScore(mol))
            except Exception:
                sa_val = 2.5
        return {
            "MW": float(Descriptors.MolWt(mol)),
            "cLogP": float(Crippen.MolLogP(mol)),
            "TPSA": float(rdMolDescriptors.CalcTPSA(mol)),
            "HBD": int(Lipinski.NumHDonors(mol)),
            "HBA": int(Lipinski.NumHAcceptors(mol)),
            "RotB": int(Lipinski.NumRotatableBonds(mol)),
            "Fsp3": float(rdMolDescriptors.CalcFractionCSP3(mol)),
            "HeavyAtoms": int(mol.GetNumHeavyAtoms()),
            "Rings": int(rdMolDescriptors.CalcNumRings(mol)),
            "AromaticRings": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
            "AliphaticRings": int(rdMolDescriptors.CalcNumAliphaticRings(mol)),
            "QED": float(QED.qed(mol)),
            "SA_Score": sa_val,
        }

    def calculate_rule_of_5(mol: Chem.Mol) -> tuple[int, List[str]]:
        violations: List[str] = []
        mw = Descriptors.MolWt(mol)
        clogp = Crippen.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        if mw > 500.0:
            violations.append(f"MW > 500 ({mw:.1f})")
        if clogp > 5.0:
            violations.append(f"cLogP > 5.0 ({clogp:.2f})")
        if hbd > 5:
            violations.append(f"HBD > 5 ({hbd})")
        if hba > 10:
            violations.append(f"HBA > 10 ({hba})")
        return len(violations), violations

    def screen_structural_alerts(
        mol: Chem.Mol, pains: bool = True, brenk: bool = True, nih: bool = False
    ) -> Dict[str, List[str]]:
        results: Dict[str, List[str]] = {"pains": [], "brenk": [], "nih": []}
        if pains:
            p_params = FilterCatalogParams()
            p_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
            p_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
            p_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
            cat = FilterCatalog(p_params)
            results["pains"] = [m.GetDescription() for m in cat.GetMatches(mol)]
        if brenk:
            b_params = FilterCatalogParams()
            b_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
            cat = FilterCatalog(b_params)
            results["brenk"] = [m.GetDescription() for m in cat.GetMatches(mol)]
        if nih:
            n_params = FilterCatalogParams()
            n_params.AddCatalog(FilterCatalogParams.FilterCatalogs.NIH)
            cat = FilterCatalog(n_params)
            results["nih"] = [m.GetDescription() for m in cat.GetMatches(mol)]
        return results


# ---------------------------------------------------------------------------
# 2. Responsive ACS-1996 SVG Renderer (following rdkit-general guidelines)
# ---------------------------------------------------------------------------
def render_responsive_svg(
    mol: Chem.Mol,
    width: int = 320,
    height: int = 220,
    highlight_atoms: Optional[List[int]] = None,
    dark_mode: bool = False,
) -> str:
    """Renders a responsive ACS-styled SVG string for web UIs."""
    rdDepictor.SetPreferCoordGen(True)
    mol_copy = Chem.Mol(mol)
    if not mol_copy.GetNumConformers():
        rdDepictor.Compute2DCoords(mol_copy)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    if dark_mode:
        rdMolDraw2D.SetDarkMode(opts)
    opts.clearBackground = not dark_mode
    opts.bondLineWidth = 1.8
    opts.multipleBondOffset = 0.18
    opts.padding = 0.08
    drawer.DrawMolecule(mol_copy, highlightAtoms=highlight_atoms or [])
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    svg = re.sub(r"<\?xml[^>]*\?>\s*", "", svg)
    svg = re.sub(r"<!DOCTYPE[^>]*>\s*", "", svg)
    svg = re.sub(
        r'(<svg\b[^>]*?)\s+width=[\'"][^\'"]*[\'"]\s+height=[\'"][^\'"]*[\'"]',
        r'\1 width="100%" height="100%"',
        svg,
        count=1,
    )
    return svg.strip()


# ---------------------------------------------------------------------------
# 3. Functional Group SMARTS Catalog for Chemle Clues
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
        if pattern is not None and safe_substruct_match(mol, pattern):
            matched.append(name)
    return matched


# ---------------------------------------------------------------------------
# 4. Curated CHEMLE Compound Library (36 Molecules)
# ---------------------------------------------------------------------------
CHEMLE_RAW_COMPOUNDS: List[Dict[str, str]] = [
    # --- Category 1: Core Cellular Metabolites & Neurotransmitters ---
    {
        "id": "caffeine",
        "name": "Caffeine",
        "Category": "Alkaloid / Metabolite",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "hint_role": "Adenosine A1/A2A receptor antagonist & methylxanthine CNS stimulant.",
        "hint_origin": "Biosynthesized from xanthosine in Coffea arabica and Camellia sinensis.",
    },
    {
        "id": "dopamine",
        "name": "Dopamine",
        "Category": "Neurotransmitter",
        "smiles": "NCCc1ccc(O)c(O)c1",
        "hint_role": "Catecholamine neurotransmitter governing reward, motivation, and motor control.",
        "hint_origin": "Synthesized from L-DOPA via aromatic L-amino acid decarboxylase (DOPA decarboxylase).",
    },
    {
        "id": "serotonin",
        "name": "Serotonin",
        "Category": "Neurotransmitter",
        "smiles": "NCCc1c[nH]c2ccc(O)cc12",
        "hint_role": "Indolamine monoamine (5-HT) regulating mood, sleep, and GI motility.",
        "hint_origin": "Derived from L-tryptophan via tryptophan hydroxylase and decarboxylation.",
    },
    {
        "id": "adrenaline",
        "name": "Adrenaline",
        "Category": "Hormone / Neurotransmitter",
        "smiles": "CNC[C@H](O)c1ccc(O)c(O)c1",
        "hint_role": "Fight-or-flight catecholamine hormone acting on alpha- and beta-adrenergic receptors.",
        "hint_origin": "Produced by adrenal medulla chromaffin cells via PNMT methylation of noradrenaline.",
    },
    {
        "id": "melatonin",
        "name": "Melatonin",
        "Category": "Hormone / Metabolite",
        "smiles": "CC(=O)NCCc1c[nH]c2ccc(OC)cc12",
        "hint_role": "Pineal gland methoxy-indole hormone synchronizing circadian rhythms.",
        "hint_origin": "Formed by N-acetylation and O-methylation of serotonin in the dark phase.",
    },
    {
        "id": "histamine",
        "name": "Histamine",
        "Category": "Biogenic Amine",
        "smiles": "NCCc1c[nH]cn1",
        "hint_role": "Imidazole biogenic amine mediating allergic mast-cell degranulation and gastric acid secretion.",
        "hint_origin": "Produced by histidine decarboxylase removal of the carboxyl group from L-histidine.",
    },
    {
        "id": "tryptophan",
        "name": "L-Tryptophan",
        "Category": "Primary Metabolite",
        "smiles": "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",
        "hint_role": "Largest standard proteinogenic amino acid and precursor to auxin, NAD+, and serotonin.",
        "hint_origin": "Encoded by the single codon UGG; contains an electron-rich indole side chain.",
    },
    {
        "id": "ascorbic_acid",
        "name": "Ascorbic Acid",
        "Category": "Cofactor / Vitamin",
        "smiles": "O=C1O[C@H]([C@@H](O)CO)C(O)=C1O",
        "hint_role": "Redox-active enediol lactone cofactor (Vitamin C) required for prolyl/lysyl hydroxylases.",
        "hint_origin": "Essential dietary micronutrient in primates due to loss of L-gulonolactone oxidase (GULO).",
    },
    {
        "id": "biotin",
        "name": "Biotin",
        "Category": "Cofactor / Vitamin",
        "smiles": "O=C(O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12",
        "hint_role": "Covalent CO2-carrier cofactor (Vitamin B7) for acetyl-CoA and pyruvate carboxylases.",
        "hint_origin": "Fused ureido-tetrahydrothiophene bicycle with ultra-high affinity for avidin/streptavidin.",
    },
    {
        "id": "citric_acid",
        "name": "Citric Acid",
        "Category": "Primary Metabolite",
        "smiles": "O=C(O)CC(O)(CC(=O)O)C(=O)O",
        "hint_role": "Achiral prochiral tricarboxylic acid at the entry point of the Krebs / TCA cycle.",
        "hint_origin": "Formed by citrate synthase condensation of oxaloacetate and acetyl-CoA.",
    },
    {
        "id": "shikimic_acid",
        "name": "Shikimic Acid",
        "Category": "Primary Metabolite",
        "smiles": "O=C(O)C1=C[C@@H](O)[C@@H](O)[C@H](O)C1",
        "hint_role": "Central cyclohexene intermediate in aromatic amino acid & phenylpropanoid biosynthesis.",
        "hint_origin": "Isolated from Illicium verum (star anise); industrial starting material for oseltamivir.",
    },
    # --- Category 2: Hydrocarbon Degradation & Microbial Metabolites (CANT-HYD / NMZ nod!) ---
    {
        "id": "catechol",
        "name": "Catechol",
        "Category": "Hydrocarbon Degradation",
        "smiles": "Oc1ccccc1O",
        "hint_role": "Central 1,2-dihydroxybenzene ring-cleavage hub in aerobic aromatic hydrocarbon degradation.",
        "hint_origin": "Cleaved via ortho-fission (CatA, intradiol) or meta-fission (DmpB/XylE, extradiol) dioxygenases.",
    },
    {
        "id": "salicylic_acid",
        "name": "Salicylic Acid",
        "Category": "Hydrocarbon Degradation",
        "smiles": "O=C(O)c1ccccc1O",
        "hint_role": "Key intermediate in upper naphthalene degradation (NahG salicylate hydroxylase) & plant defense hormone.",
        "hint_origin": "Converted to catechol by FAD-dependent salicylate 1-monooxygenase in Pseudomonas NAH7 plasmids.",
    },
    {
        "id": "protocatechuic_acid",
        "name": "Protocatechuic Acid",
        "Category": "Hydrocarbon Degradation",
        "smiles": "O=C(O)c1ccc(O)c(O)c1",
        "hint_role": "3,4-Dihydroxybenzoate funneling hub of the bacterial beta-ketoadipate pathway (PcaHG).",
        "hint_origin": "Central catabolic node for lignin-derived aromatics, phthalate, and p-hydroxybenzoate.",
    },
    {
        "id": "gentisic_acid",
        "name": "Gentisic Acid",
        "Category": "Hydrocarbon Degradation",
        "smiles": "O=C(O)c1cc(O)ccc1O",
        "hint_role": "2,5-Dihydroxybenzoate intermediate cleaved to maleylpyruvate by gentisate 1,2-dioxygenase.",
        "hint_origin": "Alternative aerobic catabolic route for naphthalene, salicylate, and 3-hydroxybenzoate.",
    },
    {
        "id": "naphthalene",
        "name": "Naphthalene",
        "Category": "Hydrocarbon Degradation",
        "smiles": "c1ccc2ccccc2c1",
        "hint_role": "Simplest fused bicyclic polycyclic aromatic hydrocarbon (PAH) model substrate.",
        "hint_origin": "Oxidized to cis-1,2-dihydroxy-1,2-dihydronaphthalene by Rieske non-heme iron dioxygenase (NahAcAd).",
    },
    {
        "id": "phenanthrene",
        "name": "Phenanthrene",
        "Category": "Hydrocarbon Degradation",
        "smiles": "c1ccc2c(c1)ccc3ccccc23",
        "hint_role": "Angular tricyclic PAH featuring the reactive K-region (C9=C10) and bay-region architecture.",
        "hint_origin": "Signature CANT-HYD / marine oil-spill biodegradation benchmark substrate.",
    },
    {
        "id": "camphor",
        "name": "Camphor",
        "Category": "Hydrocarbon Degradation",
        "smiles": "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2",
        "hint_role": "Bicyclic monoterpene ketone and canonical substrate of cytochrome P450cam (CYP101A1).",
        "hint_origin": "Regio- and stereoselectively hydroxylated at C5-exo by Pseudomonas putida P450cam.",
    },
    {
        "id": "pyocyanin",
        "name": "Pyocyanin",
        "Category": "Microbial Natural Product",
        "smiles": "Cn1c2ccccc2nc2ccc(=O)cc21",
        "hint_role": "Blue redox-active tricyclic phenazine zwitterion & quorum-sensing virulence factor.",
        "hint_origin": "Secreted by Pseudomonas aeruginosa; generates reactive oxygen species via NAD(P)H oxidation.",
    },
    {
        "id": "prodigiosin",
        "name": "Prodigiosin",
        "Category": "Microbial Natural Product",
        "smiles": "CCCCC1=CC(=C(N1)C=C2C=CC(=N2)C3=CC(=CN3)OC)C",
        "hint_role": "Blood-red tripyrrole bacterial secondary metabolite with H+/Cl- symporter & antimicrobial activity.",
        "hint_origin": "Biosynthesized via hybrid NRPS-PKS pig/red gene cluster in Serratia marcescens and marine Vibrio.",
    },
    {
        "id": "violacein",
        "name": "Violacein",
        "Category": "Microbial Natural Product",
        "smiles": "O=C1NC(=C(c2c[nH]c3ccccc23)C1=C1Nc2ccccc2C1=O)c1c[nH]c2ccc(O)cc12",
        "hint_role": "Deep-violet bis-indole pigment assembled by the 5-enzyme VioABCDE pathway from two L-tryptophans.",
        "hint_origin": "Signature metabolite of Chromobacterium violaceum and marine Pseudoalteromonas.",
    },
    {
        "id": "geosmin",
        "name": "Geosmin",
        "Category": "Microbial Natural Product",
        "smiles": "C[C@H]1CCCC2(C)[C@H]1CC[C@@]2(C)O",
        "hint_role": "Trans-decalin sesquiterpene alcohol responsible for the scent of wet earth (petrichor).",
        "hint_origin": "Formed from farnesyl diphosphate by bifunctional geosmin synthase (Sco6073) in Streptomyces.",
    },
    # --- Category 3: Plant Natural Products & Secondary Metabolites ---
    {
        "id": "vanillin",
        "name": "Vanillin",
        "Category": "Plant Natural Product",
        "smiles": "COc1cc(C=O)ccc1O",
        "hint_role": "Phenolic benzaldehyde derivative and primary flavor molecule of cured vanilla beans.",
        "hint_origin": "Also formed during microbial ferulic acid and lignin depolymerization.",
    },
    {
        "id": "resveratrol",
        "name": "Resveratrol",
        "Category": "Plant Natural Product",
        "smiles": "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1",
        "hint_role": "Trans-stilbene phytoalexin polyphenol synthesised by stilbene synthase in grapevines.",
        "hint_origin": "Combines one p-coumaroyl-CoA and three malonyl-CoA extender units via type III PKS.",
    },
    {
        "id": "curcumin",
        "name": "Curcumin",
        "Category": "Plant Natural Product",
        "smiles": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
        "hint_role": "Symmetric diarylheptanoid beta-diketone pigment giving turmeric its bright yellow fluorescence.",
        "hint_origin": "Classic PAINS / Michael acceptor benchmark compound in medicinal chemistry.",
    },
    {
        "id": "capsaicin",
        "name": "Capsaicin",
        "Category": "Plant Natural Product",
        "smiles": "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O",
        "hint_role": "Vanilloid amide agonist of the heat-activated TRPV1 ion channel nociceptor.",
        "hint_origin": "Synthesized in Capsicum placenta from vanillylamine and a branched-chain fatty acyl-CoA.",
    },
    {
        "id": "nicotine",
        "name": "Nicotine",
        "Category": "Alkaloid / Metabolite",
        "smiles": "CN1CCC[C@H]1c2cccnc2",
        "hint_role": "Pyridine-pyrrolidine alkaloid agonist of nicotinic acetylcholine receptors (nAChRs).",
        "hint_origin": "Biosynthesized in Nicotiana roots by coupling ornithine-derived N-methylpyrrolinium with nicotinic acid.",
    },
    {
        "id": "quinine",
        "name": "Quinine",
        "Category": "Alkaloid / Metabolite",
        "smiles": "COc1ccc2nccc([C@@H](O)[C@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1",
        "hint_role": "Cinchona quinoline-quinuclidine alkaloid that inhibits Plasmodium hemozoin biocrystallization.",
        "hint_origin": "Derived from strictosidine via tryptophan-terpenoid indole-to-quinoline rearrangement.",
    },
    {
        "id": "artemisinin",
        "name": "Artemisinin",
        "Category": "Plant Natural Product",
        "smiles": "C[C@@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@]23OO4",
        "hint_role": "Sesquiterpene lactone containing an essential endoperoxide bridge cleaved by heme iron(II).",
        "hint_origin": "Isolated from Artemisia annua by Tu Youyou (2015 Nobel Prize in Physiology or Medicine).",
    },
    # --- Category 4: Classic Therapeutics & Bioactive Probes ---
    {
        "id": "aspirin",
        "name": "Aspirin",
        "Category": "Therapeutic Drug",
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "hint_role": "Irreversible covalent inhibitor acetylating Ser530/Ser516 of COX-1 and COX-2.",
        "hint_origin": "O-Acetyl derivative of salicylic acid introduced clinically in 1899.",
    },
    {
        "id": "ibuprofen",
        "name": "Ibuprofen",
        "Category": "Therapeutic Drug",
        "smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
        "hint_role": "2-Arylpropionic acid reversible non-selective cyclooxygenase (COX) NSAID inhibitor.",
        "hint_origin": "Features an isobutyl-substituted phenyl ring and a chiral alpha-methyl propanoate group.",
    },
    {
        "id": "acetaminophen",
        "name": "Acetaminophen",
        "Category": "Therapeutic Drug",
        "smiles": "CC(=O)Nc1ccc(O)cc1",
        "hint_role": "Para-aminophenol analgesic/antipyretic (Paracetamol); oxidized by CYP2E1 to reactive NAPQI.",
        "hint_origin": "Minimalist 151 Da acylated phenol with zero rotatable carbon-carbon bonds.",
    },
    {
        "id": "penicillin_g",
        "name": "Penicillin G",
        "Category": "Therapeutic Drug",
        "smiles": "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
        "hint_role": "Fused beta-lactam / thiazolidine antibiotic mimicking D-Ala-D-Ala to acylate bacterial PBPs.",
        "hint_origin": "Assembled by ACV synthetase (NRPS) and isopenicillin N synthase in Penicillium chrysogenum.",
    },
    {
        "id": "ciprofloxacin",
        "name": "Ciprofloxacin",
        "Category": "Therapeutic Drug",
        "smiles": "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
        "hint_role": "Second-generation fluoroquinolone trapping bacterial DNA gyrase and topoisomerase IV.",
        "hint_origin": "Features a cyclopropyl N1 substituent, C6 fluorine, and C7 piperazine ring.",
    },
    {
        "id": "fluoxetine",
        "name": "Fluoxetine",
        "Category": "Therapeutic Drug",
        "smiles": "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",
        "hint_role": "Archetypal selective serotonin reuptake inhibitor (SSRI) blocking the SERT transporter.",
        "hint_origin": "Para-trifluoromethylphenoxy phenylpropylamine with high CNS MPO blood-brain barrier permeability.",
    },
    {
        "id": "cholesterol",
        "name": "Cholesterol",
        "Category": "Primary Metabolite",
        "smiles": "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(C)CC[C@H]3[C@H]2CC=C4[C@@]3(C)CC[C@@H](O)C4",
        "hint_role": "Tetracyclic cyclopentanoperhydrophenanthrene sterol modulating eukaryotic membrane fluidity.",
        "hint_origin": "Formed by oxidosqualene cyclase chair-boat-chair-boat cyclization of 2,3-oxidosqualene.",
    },
]


# ---------------------------------------------------------------------------
# 5. Curated SEQLE 5-Residue Functional Motif Library (30 Motifs)
# ---------------------------------------------------------------------------
AA_PROPERTIES: Dict[str, Dict[str, Any]] = {
    "A": {"name": "Alanine", "three": "Ala", "group": "Aliphatic", "hydropathy": 1.8, "charge": 0.0, "mw": 89.09},
    "V": {"name": "Valine", "three": "Val", "group": "Aliphatic", "hydropathy": 4.2, "charge": 0.0, "mw": 117.15},
    "I": {"name": "Isoleucine", "three": "Ile", "group": "Aliphatic", "hydropathy": 4.5, "charge": 0.0, "mw": 131.17},
    "L": {"name": "Leucine", "three": "Leu", "group": "Aliphatic", "hydropathy": 3.8, "charge": 0.0, "mw": 131.17},
    "M": {"name": "Methionine", "three": "Met", "group": "Aliphatic", "hydropathy": 1.9, "charge": 0.0, "mw": 149.21},
    "F": {"name": "Phenylalanine", "three": "Phe", "group": "Aromatic", "hydropathy": 2.8, "charge": 0.0, "mw": 165.19},
    "Y": {"name": "Tyrosine", "three": "Tyr", "group": "Aromatic", "hydropathy": -1.3, "charge": 0.0, "mw": 181.19},
    "W": {"name": "Tryptophan", "three": "Trp", "group": "Aromatic", "hydropathy": -0.9, "charge": 0.0, "mw": 204.23},
    "S": {"name": "Serine", "three": "Ser", "group": "Polar", "hydropathy": -0.8, "charge": 0.0, "mw": 105.09},
    "T": {"name": "Threonine", "three": "Thr", "group": "Polar", "hydropathy": -0.7, "charge": 0.0, "mw": 119.12},
    "N": {"name": "Asparagine", "three": "Asn", "group": "Polar", "hydropathy": -3.5, "charge": 0.0, "mw": 132.12},
    "Q": {"name": "Glutamine", "three": "Gln", "group": "Polar", "hydropathy": -3.5, "charge": 0.0, "mw": 146.15},
    "K": {"name": "Lysine", "three": "Lys", "group": "Positive", "hydropathy": -3.9, "charge": 1.0, "mw": 146.19},
    "R": {"name": "Arginine", "three": "Arg", "group": "Positive", "hydropathy": -4.5, "charge": 1.0, "mw": 174.20},
    "H": {"name": "Histidine", "three": "His", "group": "Positive", "hydropathy": -3.2, "charge": 0.1, "mw": 155.16},
    "D": {"name": "Aspartate", "three": "Asp", "group": "Negative", "hydropathy": -3.5, "charge": -1.0, "mw": 133.10},
    "E": {"name": "Glutamate", "three": "Glu", "group": "Negative", "hydropathy": -3.5, "charge": -1.0, "mw": 147.13},
    "C": {"name": "Cysteine", "three": "Cys", "group": "Special", "hydropathy": 2.5, "charge": 0.0, "mw": 121.16},
    "G": {"name": "Glycine", "three": "Gly", "group": "Special", "hydropathy": -0.4, "charge": 0.0, "mw": 75.07},
    "P": {"name": "Proline", "three": "Pro", "group": "Special", "hydropathy": -1.6, "charge": 0.0, "mw": 115.13},
}

SEQLE_RAW_MOTIFS: List[Dict[str, str]] = [
    {
        "seq": "GHSMG",
        "title": "Serine Hydrolase / Cutinase GXSXG Elbow",
        "pattern": "G-X-S-X-G",
        "family": "Catalytic Active Site",
        "hint_role": "Nucleophilic elbow of alpha/beta-hydrolases, lipases, PETases, and esterases.",
        "hint_detail": "Flanking glycines prevent steric clashes at the sharp beta-turn carrying the catalytic Ser.",
    },
    {
        "seq": "KDELR",
        "title": "ER Luminal Retention & Retrieval Signal",
        "pattern": "K-D-E-L-[R/K]",
        "family": "Subcellular Targeting",
        "hint_role": "Recognized in the cis-Golgi by KDEL receptors to retrograde-transport escaped ER chaperones.",
        "hint_detail": "Combines a basic N-terminal Lys with an acidic Asp-Glu pair and hydrophobic Leu.",
    },
    {
        "seq": "DEVDG",
        "title": "Caspase-3 Apoptotic Protease Cleavage Site",
        "pattern": "D-E-V-D-G",
        "family": "Protease Cleavage Motif",
        "hint_role": "Canonical P4-P1' recognition sequence cleaved after the second Asp in PARP-1 during apoptosis.",
        "hint_detail": "Strictly requires Asp at both P4 and P1 positions followed by a small P1' Gly.",
    },
    {
        "seq": "RGDVF",
        "title": "Fibronectin Integrin-Binding Motif",
        "pattern": "R-G-D-V-F",
        "family": "Cell Adhesion",
        "hint_role": "Extracellular matrix tripeptide core recognized by alphaV-beta3 and alpha5-beta1 integrins.",
        "hint_detail": "Starts with a cationic Arg, flexible Gly spacer, and anionic Asp; basis of cilengitide.",
    },
    {
        "seq": "LPETG",
        "title": "Sortase A Cell-Wall Anchoring Motif",
        "pattern": "L-P-X-T-G",
        "family": "Transpeptidation Tag",
        "hint_role": "Cleaved between Thr and Gly by Staphylococcus aureus Sortase A and ligated to pentaglycine.",
        "hint_detail": "Widely used for site-specific enzymatic protein bioconjugation ('sortagging').",
    },
    {
        "seq": "HEIGH",
        "title": "Zincin Metalloprotease Catalytic Core",
        "pattern": "H-E-X-X-H",
        "family": "Metal-Binding Catalytic Site",
        "hint_role": "Two histidines coordinate the catalytic Zn2+ ion while glutamate acts as general base.",
        "hint_detail": "Found in thermolysin, matrix metalloproteinases (MMPs), and botulinum neurotoxin.",
    },
    {
        "seq": "CCHHC",
        "title": "Retrotransposon / Gag Zinc Knuckle Core",
        "pattern": "C-C-H-H-C",
        "family": "Nucleic Acid Binding",
        "hint_role": "Tetrahedral Zn2+-coordinating Cys/His cluster that binds single-stranded viral RNA.",
        "hint_detail": "Dominated by special thiol-bearing Cys and imidazole His residues.",
    },
    {
        "seq": "YIGSR",
        "title": "Laminin B1 Chain Receptor-Binding Pentapeptide",
        "pattern": "Y-I-G-S-R",
        "family": "Cell Adhesion",
        "hint_role": "Basement membrane laminin-1 motif promoting epithelial cell attachment and inhibiting metastasis.",
        "hint_detail": "Begins with an aromatic Tyr and hydrophobic Ile, ending with a basic Arg.",
    },
    {
        "seq": "NPVYK",
        "title": "Clathrin Endocytosis NPXY Internalization Signal",
        "pattern": "N-P-X-Y-[K/R]",
        "family": "Membrane Trafficking",
        "hint_role": "Cytoplasmic tail motif of LDL, insulin, and integrin receptors recruiting PTB-domain adaptors.",
        "hint_detail": "Asn-Pro forms a tight beta-turn positioning the aromatic Tyr into the adaptor pocket.",
    },
    {
        "seq": "CVIMK",
        "title": "K-Ras4B C-Terminal CAAX Prenylation Box",
        "pattern": "C-A-A-X",
        "family": "Post-Translational Lipidation",
        "hint_role": "Farnesyltransferase attaches a C15 isoprenoid lipid to the conserved Cys thiol.",
        "hint_detail": "Terminal Met directs high-affinity recognition by farnesyltransferase over GGTase-I.",
    },
    {
        "seq": "DEADR",
        "title": "DEAD-Box RNA Helicase Walker B Motif",
        "pattern": "D-E-A-D-[R/K]",
        "family": "NTPase / Helicase Core",
        "hint_role": "Coordinates Mg2+ and polarizes the attacking water molecule for ATP hydrolysis in RNA unwinding.",
        "hint_detail": "Iconic sequence motif that gave the entire DEAD-box helicase superfamily its name.",
    },
    {
        "seq": "GAGTG",
        "title": "Walker A / P-Loop Nucleotide-Binding Core",
        "pattern": "G-X-G-X-G",
        "family": "NTPase / Helicase Core",
        "hint_role": "Glycine-rich phosphate-binding loop wrapping around beta/gamma phosphates of ATP or GTP.",
        "hint_detail": "Alternating glycines allow backbone amides to form a tight anion-binding nest.",
    },
    {
        "seq": "VDFGL",
        "title": "Protein Kinase Activation Loop DFG Motif",
        "pattern": "[V/I]-D-F-G-[L/M]",
        "family": "Kinase Switch",
        "hint_role": "Conformational switch (DFG-in vs. DFG-out) governing Mg2+-ATP coordination in all eukaryotic kinases.",
        "hint_detail": "Target of Type II kinase inhibitors like imatinib; contains Asp-Phe-Gly at positions 2-4.",
    },
    {
        "seq": "HRDLA",
        "title": "Protein Kinase Catalytic Loop HRD Motif",
        "pattern": "H-R-D-L-A",
        "family": "Catalytic Active Site",
        "hint_role": "Catalytic Asp acts as proton acceptor for the substrate hydroxyl group during phosphoryl transfer.",
        "hint_detail": "Arg neutralizes the phosphorylated activation loop in RD-class tyrosine and Ser/Thr kinases.",
    },
    {
        "seq": "WSEWS",
        "title": "Type I Cytokine Receptor WSXWS Box",
        "pattern": "W-S-X-W-S",
        "family": "Structural Folding Motif",
        "hint_role": "Forms a Tryptophan-Arginine zipper ('Trp-stack') essential for folding EPO, IL-2, and GH receptors.",
        "hint_detail": "Symmetric pair of Trp-Ser dipeptides separated by a central acidic Glu.",
    },
    {
        "seq": "PPPVP",
        "title": "Polyproline Type II Helix SH3 Ligand",
        "pattern": "P-X-X-P",
        "family": "Protein-Protein Interaction",
        "hint_role": "Left-handed PPII helix docking into aromatic grooves of Src Homology 3 (SH3) signaling domains.",
        "hint_detail": "Four rigid cyclic Pro residues flanking a single branched aliphatic Val.",
    },
    {
        "seq": "INFEK",
        "title": "Ovalbumin SIINFEKL MHC Class I Epitope Core",
        "pattern": "I-N-F-E-K",
        "family": "Immunology Epitope",
        "hint_role": "Central TCR-contact and H-2Kb anchor core of the benchmark CD8+ T-cell octapeptide SIINFEKL.",
        "hint_detail": "Aromatic Phe at center with C-terminal acidic Glu and basic Lys.",
    },
    {
        "seq": "YQRLN",
        "title": "Trans-Golgi Network Tyrosine Sorting Motif",
        "pattern": "Y-X-X-[L/I/V]",
        "family": "Membrane Trafficking",
        "hint_role": "Binds the mu-2 subunit of AP-2 clathrin adaptor complexes for rapid endocytic sorting.",
        "hint_detail": "Starts with an essential unphosphorylated Tyr and places hydrophobic Leu at Y+3.",
    },
    {
        "seq": "KKMPQ",
        "title": "COPI Coatomer Dilysine ER Retrieval Signal",
        "pattern": "K-K-X-X",
        "family": "Subcellular Targeting",
        "hint_role": "C-terminal cytosolic tail motif of ER membrane proteins recognized by alpha/beta'-COP.",
        "hint_detail": "Starts with a twin positively charged Lys-Lys dyad.",
    },
    {
        "seq": "CAACH",
        "title": "Cytochrome c CxxCH Heme-Binding Motif",
        "pattern": "C-X-X-C-H",
        "family": "Cofactor Covalent Attachment",
        "hint_role": "Two cysteines form thioether bonds to heme vinyl groups while His serves as axial Fe ligand.",
        "hint_detail": "Signature motif of c-type cytochromes in respiratory and marine electron transport chains.",
    },
    {
        "seq": "NLTGY",
        "title": "Eukaryotic N-Linked Glycosylation Sequon",
        "pattern": "N-X-[S/T]",
        "family": "Post-Translational Modification",
        "hint_role": "Oligosaccharyltransferase transfers Glc3Man9GlcNAc2 onto the Asn amide nitrogen in the ER.",
        "hint_detail": "Requires Asn at position 1, any non-Pro at position 2, and Thr/Ser at position 3.",
    },
    {
        "seq": "DSGLS",
        "title": "Beta-TrCP SCF Ubiquitin Ligase Phosphodegron",
        "pattern": "D-S-G-X-X-S",
        "family": "Proteasomal Degradation",
        "hint_role": "Dual phosphorylation of the Ser residues triggers ubiquitination of I-kappa-B and beta-catenin.",
        "hint_detail": "Contains an acidic Asp followed by Ser-Gly-Leu-Ser.",
    },
    {
        "seq": "KENKR",
        "title": "APC/C Cdh1 Mitotic Exit KEN-Box Degron",
        "pattern": "K-E-N-X-R",
        "family": "Proteasomal Degradation",
        "hint_role": "Targets cell-cycle kinases (Aurora B, Nek2, Cdc6) for destruction during late mitosis and G1.",
        "hint_detail": "Named directly after its N-terminal Lys-Glu-Asn tripeptide.",
    },
    {
        "seq": "PKKKR",
        "title": "SV40 Large T-Antigen Nuclear Localization Signal",
        "pattern": "P-K-K-K-R",
        "family": "Subcellular Targeting",
        "hint_role": "Prototype monopartite classical NLS recognized by Importin-alpha for nuclear pore import.",
        "hint_detail": "Initiated by a helix-breaking Pro followed by four consecutive cationic residues.",
    },
    {
        "seq": "FGGFD",
        "title": "Nociceptin / Orphanin FQ Neuropeptide Core",
        "pattern": "F-G-G-F-[D/T]",
        "family": "Neuropeptide GPCR Message",
        "hint_role": "N-terminal 'message' domain activating the NOP opioid-related GPCR (vs. YGGFM in enkephalin).",
        "hint_detail": "Two Phe rings separated by a flexible Gly-Gly spacer.",
    },
    {
        "seq": "YGGFM",
        "title": "Met-Enkephalin Endogenous Opioid Pentapeptide",
        "pattern": "Y-G-G-F-M",
        "family": "Neuropeptide GPCR Message",
        "hint_role": "Endogenous delta/mu-opioid receptor agonist pentapeptide cleaved from proenkephalin.",
        "hint_detail": "N-terminal Tyr1 phenol and amine mimic the tyramine moiety of morphine.",
    },
    {
        "seq": "SYGYS",
        "title": "FUS / FET Low-Complexity Phase-Separation Motif",
        "pattern": "[S/G]-Y-[G/S]-Y-S",
        "family": "Biomolecular Condensate",
        "hint_role": "Pi-pi and cation-pi sticker-and-spacer repeat driving liquid-liquid phase separation of stress granules.",
        "hint_detail": "Palindrome of polar Ser, flexible Gly, and aromatic Tyr stickers.",
    },
    {
        "seq": "LDLAS",
        "title": "CRM1 / Exportin-1 Leucine-Rich Nuclear Export Core",
        "pattern": "L-X-L-A-S",
        "family": "Subcellular Targeting",
        "hint_role": "Hydrophobic ridge docking into the NES cleft of Exportin-1 (XPO1) in the presence of RanGTP.",
        "hint_detail": "Two Leu residues separated by an acidic Asp.",
    },
    {
        "seq": "NPLYR",
        "title": "Integrin Beta-Tail Talin-Binding PTB Motif",
        "pattern": "N-P-L-Y-R",
        "family": "Cell Adhesion / Inside-Out Signaling",
        "hint_role": "Binds the FERM domain of talin to switch integrins from bent-inactive to extended-active states.",
        "hint_detail": "Classic Asn-Pro-Leu-Tyr beta-turn capped by a C-terminal Arg.",
    },
    {
        "seq": "HPHPH",
        "title": "Histidine-Proline Cu2+/Zn2+ Antimicrobial Repeat",
        "pattern": "H-P-H-P-H",
        "family": "Metal Chelation / Host Defense",
        "hint_role": "Conformationally restricted metal-chelating motif found in histidine-rich glycoproteins.",
        "hint_detail": "Strictly alternating His and Pro residues.",
    },
]


def compute_peptide_biochem(seq: str) -> Dict[str, Any]:
    """Computes sequence-level biochemical properties for a 5-AA peptide."""
    gravy = sum(AA_PROPERTIES[aa]["hydropathy"] for aa in seq) / len(seq)
    # Net charge at pH 7.4 (including N-terminal +0.75 and C-terminal -0.99 or sidechain net)
    sidechain_charge = sum(AA_PROPERTIES[aa]["charge"] for aa in seq)
    aromatic_count = sum(1 for aa in seq if AA_PROPERTIES[aa]["group"] == "Aromatic")
    polar_count = sum(
        1 for aa in seq if AA_PROPERTIES[aa]["group"] in ("Polar", "Positive", "Negative")
    )
    # Approximate isoelectric point (pI) via Henderson-Hasselbalch bisection
    pka_nterm = 8.0
    pka_cterm = 3.1
    pka_side = {"D": (3.9, -1), "E": (4.1, -1), "C": (8.3, -1), "Y": (10.1, -1), "H": (6.0, 1), "K": (10.5, 1), "R": (12.5, 1)}

    def net_q(ph: float) -> float:
        q = (1.0 / (1.0 + 10 ** (ph - pka_nterm))) - (1.0 / (1.0 + 10 ** (pka_cterm - ph)))
        for aa in seq:
            if aa in pka_side:
                pk, sign = pka_side[aa]
                if sign > 0:
                    q += 1.0 / (1.0 + 10 ** (ph - pk))
                else:
                    q -= 1.0 / (1.0 + 10 ** (pk - ph))
        return q

    lo, hi = 1.0, 13.5
    for _ in range(30):
        mid = 0.5 * (lo + hi)
        if net_q(mid) > 0:
            lo = mid
        else:
            hi = mid
    pi_val = round(0.5 * (lo + hi), 2)

    return {
        "gravy": round(gravy, 2),
        "net_charge_ph74": round(sidechain_charge, 1),
        "isoelectric_point": pi_val,
        "aromatic_count": aromatic_count,
        "polar_count": polar_count,
        "residue_groups": [AA_PROPERTIES[aa]["group"] for aa in seq],
    }


def build_all() -> Dict[str, Any]:
    """Compiles and validates the complete Chemle + Seqle puzzle database."""
    chemle_records: List[Dict[str, Any]] = []
    chemle_fps: List[Any] = []

    print(f"[BioWordle Build] Validating & compiling {len(CHEMLE_RAW_COMPOUNDS)} Chemle molecules with RDKit...")
    for item in CHEMLE_RAW_COMPOUNDS:
        raw_mol = Chem.MolFromSmiles(item["smiles"])
        if raw_mol is None:
            raise ValueError(f"Invalid SMILES for {item['id']}: {item['smiles']}")
        mol = standardize_molecule(raw_mol, uncharge=False)
        if mol is None:
            mol = raw_mol

        desc = compute_standard_descriptors(mol, style="short", include_sa_score=True)
        cns_mpo = calculate_cns_mpo(mol)
        ro5_viols, ro5_list = calculate_rule_of_5(mol)
        alerts = screen_structural_alerts(mol, pains=True, brenk=True, nih=False)
        fg_list = detect_functional_groups(mol)

        # Full 2D ACS depiction
        svg_full = render_responsive_svg(mol, width=320, height=220, dark_mode=False)

        # Bemis-Murcko Scaffold depiction for progressive clue
        scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
        if scaffold_mol is not None and scaffold_mol.GetNumHeavyAtoms() >= 3:
            svg_scaffold = render_responsive_svg(scaffold_mol, width=320, height=220, dark_mode=False)
            has_scaffold = True
        else:
            svg_scaffold = svg_full
            has_scaffold = False

        fp = get_morgan_fingerprint(mol, radius=2, n_bits=2048)
        chemle_fps.append(fp)

        record = {
            "id": item["id"],
            "name": item["name"],
            "Category": item["Category"],
            "canonical_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "exact_mass": round(float(Descriptors.ExactMolWt(mol)), 2),
            "MW": round(float(desc["MW"]), 1),
            "cLogP": round(float(desc["cLogP"]), 2),
            "TPSA": round(float(desc["TPSA"]), 1),
            "HBD": int(desc["HBD"]),
            "HBA": int(desc["HBA"]),
            "RotB": int(desc["RotB"]),
            "Fsp3": round(float(desc["Fsp3"]), 2),
            "HeavyAtoms": int(desc["HeavyAtoms"]),
            "Rings": int(desc["Rings"]),
            "AromaticRings": int(desc["AromaticRings"]),
            "AliphaticRings": int(desc["AliphaticRings"]),
            "QED": round(float(desc["QED"]), 2),
            "SA_Score": round(float(desc["SA_Score"]), 2),
            "CNS_MPO": round(float(cns_mpo), 2),
            "Ro5_Violations": int(ro5_viols),
            "Ro5_Details": ro5_list,
            "PAINS_Alerts": alerts["pains"],
            "Brenk_Alerts": alerts["brenk"],
            "functional_groups": fg_list,
            "hint_role": item["hint_role"],
            "hint_origin": item["hint_origin"],
            "has_scaffold": has_scaffold,
            "svg_scaffold": svg_scaffold,
            "svg_full": svg_full,
        }
        chemle_records.append(record)

    # Compute N x N pairwise Morgan Tanimoto similarity matrix (0 - 100%)
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

    print(f"[BioWordle Build] Validating & compiling {len(SEQLE_RAW_MOTIFS)} Seqle 5-AA peptide motifs with RDKit...")
    seqle_records: List[Dict[str, Any]] = []
    for item in SEQLE_RAW_MOTIFS:
        seq = item["seq"].strip().upper()
        if len(seq) != 5 or any(ch not in AA_PROPERTIES for ch in seq):
            raise ValueError(f"Invalid 5-AA sequence: {seq}")
        pep_mol = Chem.MolFromSequence(seq)
        if pep_mol is None:
            raise ValueError(f"RDKit MolFromSequence failed for: {seq}")

        pep_desc = compute_standard_descriptors(pep_mol, style="short", include_sa_score=False)
        biochem = compute_peptide_biochem(seq)
        svg_peptide = render_responsive_svg(pep_mol, width=380, height=210, dark_mode=False)

        seqle_records.append(
            {
                "seq": seq,
                "title": item["title"],
                "pattern": item["pattern"],
                "family": item["family"],
                "hint_role": item["hint_role"],
                "hint_detail": item["hint_detail"],
                "formula": rdMolDescriptors.CalcMolFormula(pep_mol),
                "MW": round(float(pep_desc["MW"]), 1),
                "cLogP": round(float(pep_desc["cLogP"]), 2),
                "TPSA": round(float(pep_desc["TPSA"]), 1),
                "HBD": int(pep_desc["HBD"]),
                "HBA": int(pep_desc["HBA"]),
                "HeavyAtoms": int(pep_desc["HeavyAtoms"]),
                **biochem,
                "svg_peptide": svg_peptide,
            }
        )

    return {
        "metadata": {
            "generator": "RDKit Build-Time Puzzle Compiler (rdkit-general + cheminformatics-core)",
            "rdkit_version": Chem.rdBase.rdkitVersion,
            "chemle_count": len(chemle_records),
            "seqle_count": len(seqle_records),
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

    # Atomic writes as recommended in rdkit_applications_and_visualizers.md
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
        f"[BioWordle Build] Done! Wrote {len(payload['chemle'])} Chemle compounds "
        f"and {len(payload['seqle'])} Seqle motifs to:\n"
        f"  - {json_path}\n"
        f"  - {js_path}"
    )

    # Also emit a single-file standalone HTML bundle if index.html, styles.css, and app.js exist
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

        # If --artifact-output is supplied, atomically update that file too
        if len(sys.argv) > 1 and sys.argv[1] == "--artifact-output" and len(sys.argv) > 2:
            art_path = sys.argv[2]
            tmp_art = art_path + ".tmp"
            with open(tmp_art, "w", encoding="utf-8") as f:
                f.write(bundled)
            os.replace(tmp_art, art_path)
            print(f"  - {art_path} (Artifact preview updated)")



if __name__ == "__main__":
    main()
