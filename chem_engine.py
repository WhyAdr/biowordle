"""Project-owned RDKit Cheminformatics, Sequence Biophysics & Information-Theory Engine.

Zero external plugin dependencies: uses only public RDKit APIs and the Python
standard library so builds are 100% identical across developer machines, CI, and HPC.
"""

from __future__ import annotations

import math
import os
import re
import sys
from collections import Counter
from typing import Any, Dict, List, Optional, Tuple

import rdkit
from rdkit import Chem, DataStructs, RDConfig
from rdkit.Chem import (
    Crippen,
    Descriptors,
    Lipinski,
    QED,
    rdDepictor,
    rdFingerprintGenerator,
    rdMolDescriptors,
)
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

# Load RDKit's bundled Ertl & Schuffenhauer SA_Score contrib module
_CONTRIB_SA = os.path.join(RDConfig.RDContribDir, "SA_Score")
if os.path.isdir(_CONTRIB_SA) and _CONTRIB_SA not in sys.path:
    sys.path.insert(0, _CONTRIB_SA)
try:
    import sascorer as _sascorer  # type: ignore
except Exception:
    _sascorer = None

# ---------------------------------------------------------------------------
# 1. Explicit Fingerprint & Standardization Specifications
# ---------------------------------------------------------------------------
FINGERPRINT_SPEC: Dict[str, Any] = {
    "algorithm": "Morgan (ECFP4-like circular fingerprint)",
    "radius": 2,
    "fpSize": 2048,
    "includeChirality": True,
    "useBondTypes": True,
    "metric": "Tanimoto coefficient T(A,B) = |A cap B| / |A cup B|",
}

STANDARDIZATION_SPEC: Dict[str, Any] = {
    "cleanup": True,
    "fragment_parent": True,
    "uncharge": False,
    "canonical_tautomer": False,
    "preserve_isomeric_stereo": True,
    "rationale": (
        "Validates valence, strips counterions if present, and canonicalizes aromaticity "
        "while preserving curated keto/enol states and stereocenters (R/S, E/Z)."
    ),
}

PEPTIDE_PKA_SPEC: Dict[str, Any] = {
    "model": "Henderson-Hasselbalch Monoprotic Additive Model at pH 7.4",
    "ph": 7.4,
    "n_terminus_pka": 8.0,
    "c_terminus_pka": 3.1,
    "sidechain_pka": {
        "D": {"pka": 3.9, "charge_state": -1},
        "E": {"pka": 4.1, "charge_state": -1},
        "C": {"pka": 8.3, "charge_state": -1},
        "Y": {"pka": 10.1, "charge_state": -1},
        "H": {"pka": 6.0, "charge_state": +1},
        "K": {"pka": 10.5, "charge_state": +1},
        "R": {"pka": 12.5, "charge_state": +1},
    },
    "note": (
        "Computes both side-chain formal charge sum and full Henderson-Hasselbalch net charge "
        "(including free N-terminal alpha-ammonium and C-terminal alpha-carboxylate groups)."
    ),
}

_MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(
    radius=FINGERPRINT_SPEC["radius"],
    fpSize=FINGERPRINT_SPEC["fpSize"],
    includeChirality=FINGERPRINT_SPEC["includeChirality"],
    useBondTypes=FINGERPRINT_SPEC["useBondTypes"],
)

_PAINS_PARAMS = FilterCatalogParams()
_PAINS_PARAMS.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
_PAINS_PARAMS.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
_PAINS_PARAMS.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
_PAINS_CATALOG = FilterCatalog(_PAINS_PARAMS)

_BRENK_PARAMS = FilterCatalogParams()
_BRENK_PARAMS.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
_BRENK_CATALOG = FilterCatalog(_BRENK_PARAMS)


def standardize_and_audit(smiles: str) -> Tuple[Chem.Mol, Dict[str, Any]]:
    """Standardizes a curated SMILES string transparently and records provenance."""
    raw_mol = Chem.MolFromSmiles(smiles)
    if raw_mol is None:
        raise ValueError(f"RDKit failed to parse SMILES: {smiles!r}")

    curated_canonical = Chem.MolToSmiles(raw_mol, isomericSmiles=True)
    curated_inchikey = Chem.MolToInchiKey(raw_mol)

    clean_mol = rdMolStandardize.Cleanup(raw_mol)
    parent_mol = rdMolStandardize.FragmentParent(clean_mol)
    Chem.AssignStereochemistry(parent_mol, cleanIt=False, force=True)

    std_canonical = Chem.MolToSmiles(parent_mol, isomericSmiles=True)
    std_inchikey = Chem.MolToInchiKey(parent_mol)

    audit = {
        "curated_smiles": smiles,
        "curated_canonical_smiles": curated_canonical,
        "curated_inchikey": curated_inchikey,
        "canonical_smiles": std_canonical,
        "inchikey": std_inchikey,
        "modified_by_standardization": curated_inchikey != std_inchikey,
    }
    return parent_mol, audit


def compute_chiral_morgan_fp(mol: Chem.Mol) -> Any:
    """Computes a 2048-bit radius-2 chiral Morgan fingerprint."""
    return _MORGAN_GEN.GetFingerprint(mol)


def compute_descriptors_and_insights(mol: Chem.Mol) -> Dict[str, Any]:
    """Computes validated 2D RDKit descriptors and post-game medicinal chemistry insights."""
    mw = float(Descriptors.MolWt(mol))
    exact_mass = float(Descriptors.ExactMolWt(mol))
    clogp = float(Crippen.MolLogP(mol))
    tpsa = float(rdMolDescriptors.CalcTPSA(mol))
    hbd = int(Lipinski.NumHDonors(mol))
    hba = int(Lipinski.NumHAcceptors(mol))
    rotb = int(Lipinski.NumRotatableBonds(mol))
    fsp3 = float(rdMolDescriptors.CalcFractionCSP3(mol))
    heavy = int(mol.GetNumHeavyAtoms())
    rings = int(rdMolDescriptors.CalcNumRings(mol))
    arom_rings = int(rdMolDescriptors.CalcNumAromaticRings(mol))
    aliph_rings = int(rdMolDescriptors.CalcNumAliphaticRings(mol))
    chiral_centers = len(
        Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    )
    qed_val = float(QED.qed(mol))

    if _sascorer is None:
        raise RuntimeError("RDKit SA_Score contrib module is required for reproducible builds.")
    sa_val = float(_sascorer.calculateScore(mol))

    # Heuristic CNS-MPO desirability (Wager et al. 2010 monotonic ramps with SMARTS amine pKa proxy)
    basic_amine = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=[O,S]);!$(NS=O);!$(Nc);!$(N#*)]")
    has_basic_amine = bool(mol.HasSubstructMatch(basic_amine)) if basic_amine else False
    est_pka = 9.5 if has_basic_amine else 5.0
    est_clogd = clogp - (est_pka - 7.4) if est_pka > 7.4 else clogp

    def _ramp(x: float, lo: float, hi: float) -> float:
        return min(1.0, max(0.0, (hi - x) / (hi - lo)))

    s_tpsa = (
        1.0
        if 40.0 <= tpsa <= 90.0
        else ((tpsa - 20.0) / 20.0 if 20.0 < tpsa < 40.0 else (_ramp(tpsa, 90.0, 120.0) if 90.0 < tpsa < 120.0 else 0.0))
    )
    cns_mpo = round(
        _ramp(clogp, 3.0, 5.0)
        + _ramp(est_clogd, 2.0, 4.0)
        + _ramp(mw, 360.0, 500.0)
        + s_tpsa
        + _ramp(float(hbd), 0.5, 3.5)
        + _ramp(est_pka, 8.0, 10.0),
        2,
    )

    ro5_violations: List[str] = []
    if mw > 500.0:
        ro5_violations.append(f"MW > 500 ({mw:.1f})")
    if clogp > 5.0:
        ro5_violations.append(f"cLogP > 5.0 ({clogp:.2f})")
    if hbd > 5:
        ro5_violations.append(f"HBD > 5 ({hbd})")
    if hba > 10:
        ro5_violations.append(f"HBA > 10 ({hba})")

    pains_hits = [m.GetDescription() for m in _PAINS_CATALOG.GetMatches(mol)]
    brenk_hits = [m.GetDescription() for m in _BRENK_CATALOG.GetMatches(mol)]

    values = [mw, exact_mass, clogp, tpsa, fsp3, qed_val, sa_val, cns_mpo]
    if any(math.isnan(v) or math.isinf(v) for v in values):
        raise ValueError(f"NaN or Inf descriptor encountered for {Chem.MolToSmiles(mol)}")

    return {
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "exact_mass": round(exact_mass, 2),
        "MW": round(mw, 1),
        "cLogP": round(clogp, 2),
        "TPSA": round(tpsa, 1),
        "HBD": hbd,
        "HBA": hba,
        "RotB": rotb,
        "Fsp3": round(fsp3, 2),
        "HeavyAtoms": heavy,
        "Rings": rings,
        "AromaticRings": arom_rings,
        "AliphaticRings": aliph_rings,
        "ChiralCenters": chiral_centers,
        "QED": round(qed_val, 2),
        "SA_Score": round(sa_val, 2),
        "CNS_MPO": cns_mpo,
        "CNS_MPO_Method": "Heuristic Wager 6-term desirability (SMARTS basic-amine pKa proxy)",
        "Ro5_Violations": len(ro5_violations),
        "Ro5_Details": ro5_violations,
        "PAINS_Alerts": pains_hits,
        "Brenk_Alerts": brenk_hits,
    }


def render_responsive_svg(
    mol: Chem.Mol,
    width: int = 320,
    height: int = 220,
    highlight_atoms: Optional[List[int]] = None,
    dark_mode: bool = False,
) -> str:
    """Renders an ACS-1996 styled responsive SVG string."""
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


def build_scaffold_clue_svg(mol: Chem.Mol) -> Tuple[str, str, str]:
    """Builds a progressive intermediate 2D clue SVG without leaking acyclic structures.

    Returns:
      (svg_scaffold, scaffold_smiles, scaffold_mode)
      where `scaffold_mode` is:
        - "bemis_murcko": Ring-containing molecule with a distinct Bemis-Murcko scaffold
        - "ring_core": Molecule equals its own unsubstituted ring framework (e.g. Naphthalene)
        - "carbon_skeleton": Acyclic molecule (e.g. Citric acid, Pyruvic acid); strips
          terminal heteroatoms/substituents to show the stripped carbon chain skeleton
          rather than revealing the full functional molecule!
    """
    scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold_mol is not None and scaffold_mol.GetNumHeavyAtoms() >= 3:
        scaf_smi = Chem.MolToSmiles(scaffold_mol, isomericSmiles=False)
        full_flat_smi = Chem.MolToSmiles(mol, isomericSmiles=False)
        mode = "bemis_murcko" if scaf_smi != full_flat_smi else "ring_core"
        return render_responsive_svg(scaffold_mol, width=320, height=220), scaf_smi, mode

    # Acyclic fallback: build an unfunctionalized carbon-chain skeleton so Clue 2
    # shows chain connectivity without leaking carbonyls/hydroxyls
    rw = Chem.RWMol()
    old_to_new: Dict[int, int] = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            old_to_new[atom.GetIdx()] = rw.AddAtom(Chem.Atom(6))
    for bond in mol.GetBonds():
        b_idx = bond.GetBeginAtomIdx()
        e_idx = bond.GetEndAtomIdx()
        if b_idx in old_to_new and e_idx in old_to_new:
            rw.AddBond(old_to_new[b_idx], old_to_new[e_idx], Chem.BondType.SINGLE)
    skel_mol = rw.GetMol()
    if skel_mol.GetNumHeavyAtoms() >= 2:
        Chem.SanitizeMol(skel_mol)
        skel_smi = Chem.MolToSmiles(skel_mol)
        return render_responsive_svg(skel_mol, width=320, height=220), skel_smi, "carbon_skeleton"

    return render_responsive_svg(mol, width=320, height=220), "", "minimal"


# ---------------------------------------------------------------------------
# 2. Seqle Two-Pass Wordle Evaluator & Documented Peptide Biophysics
# ---------------------------------------------------------------------------
def evaluate_seqle_guess(
    guess_seq: str, target_seq: str, aa_props: Dict[str, Dict[str, Any]]
) -> List[str]:
    """Deterministic 3-tier Seqle feedback with strict two-pass duplicate handling.

    Precedence:
      1. 'exact'   — guess[i] == target[i] (consumes target[i])
      2. 'present' — guess[i] occurs at an unconsumed position in target (consumes 1 instance)
      3. 'class'   — guess[i] is not 'exact' or 'present', and primary biochemical group
                     of guess[i] equals primary biochemical group of target[i]
      4. 'absent'  — no match
    """
    if len(guess_seq) != 5 or len(target_seq) != 5:
        raise ValueError("Both guess_seq and target_seq must be 5 amino acids.")

    evals = ["absent"] * 5
    remaining: List[Optional[str]] = list(target_seq)

    # Pass 1: Exact position matches
    for i in range(5):
        if guess_seq[i] == target_seq[i]:
            evals[i] = "exact"
            remaining[i] = None

    # Pass 2: Misplaced residue matches (respecting multiplicity)
    for i in range(5):
        if evals[i] != "absent":
            continue
        ch = guess_seq[i]
        if ch in remaining:
            idx = remaining.index(ch)
            evals[i] = "present"
            remaining[idx] = None

    # Pass 3: Positional primary biochemical class match
    for i in range(5):
        if evals[i] == "absent":
            if aa_props[guess_seq[i]]["group"] == aa_props[target_seq[i]]["group"]:
                evals[i] = "class"

    return evals


def compute_peptide_biochem(seq: str, aa_props: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
    """Computes sequence biophysics under the documented Henderson-Hasselbalch model."""
    gravy = sum(aa_props[aa]["hydropathy"] for aa in seq) / len(seq)
    formal_sidechain_q = sum(aa_props[aa]["charge"] for aa in seq)
    aromatic_count = sum(1 for aa in seq if "Aromatic" in aa_props[aa]["traits"])
    polar_count = sum(
        1 for aa in seq if aa_props[aa]["group"] in ("Polar", "Positive", "Negative")
    )

    pka_nterm = PEPTIDE_PKA_SPEC["n_terminus_pka"]
    pka_cterm = PEPTIDE_PKA_SPEC["c_terminus_pka"]
    pka_side = PEPTIDE_PKA_SPEC["sidechain_pka"]

    def net_q(ph: float) -> float:
        q = (1.0 / (1.0 + 10.0 ** (ph - pka_nterm))) - (
            1.0 / (1.0 + 10.0 ** (pka_cterm - ph))
        )
        for aa in seq:
            if aa in pka_side:
                pk = pka_side[aa]["pka"]
                sign = pka_side[aa]["charge_state"]
                if sign > 0:
                    q += 1.0 / (1.0 + 10.0 ** (ph - pk))
                else:
                    q -= 1.0 / (1.0 + 10.0 ** (pk - ph))
        return q

    charge_ph74 = round(net_q(7.4), 2)

    lo, hi = 1.0, 13.5
    for _ in range(35):
        mid = 0.5 * (lo + hi)
        if net_q(mid) > 0:
            lo = mid
        else:
            hi = mid
    pi_val = round(0.5 * (lo + hi), 2)

    return {
        "gravy": round(gravy, 2),
        "net_charge_ph74": charge_ph74,
        "formal_sidechain_charge": round(formal_sidechain_q, 1),
        "isoelectric_point": pi_val,
        "aromatic_count": aromatic_count,
        "polar_count": polar_count,
        "residue_groups": [aa_props[aa]["group"] for aa in seq],
    }


# ---------------------------------------------------------------------------
# 3. Build-Time Shannon Entropy Information Gain & Puzzle Difficulty Evaluator
# ---------------------------------------------------------------------------
def _delta_code(g_val: float, t_val: float, exact_tol: float, close_tol: float) -> str:
    diff = abs(g_val - t_val)
    if diff <= exact_tol:
        return "E"
    bucket = "C" if diff <= close_tol else "F"
    direction = "U" if t_val > g_val else "D"
    return bucket + direction


def chemle_feedback_signature(
    guess: Dict[str, Any], target: Dict[str, Any], sim_pct: int
) -> Tuple[Any, ...]:
    """Returns the discrete feedback tuple f(g, t) observed by a player in Chemle."""
    sim_bucket = "100" if sim_pct == 100 else ("HI" if sim_pct >= 35 else ("MD" if sim_pct >= 18 else "LO"))
    return (
        guess["Category"] == target["Category"],
        sim_bucket,
        _delta_code(guess["MW"], target["MW"], 5.0, 35.0),
        _delta_code(guess["cLogP"], target["cLogP"], 0.25, 1.0),
        _delta_code(guess["TPSA"], target["TPSA"], 5.0, 25.0),
        _delta_code(guess["HBD"], target["HBD"], 0, 1),
        _delta_code(guess["HBA"], target["HBA"], 0, 1),
        _delta_code(guess["Rings"], target["Rings"], 0, 1),
        _delta_code(guess["AromaticRings"], target["AromaticRings"], 0, 1),
    )


def compute_chemle_information_metrics(
    records: List[Dict[str, Any]], sim_matrix: Dict[str, Dict[str, int]]
) -> Dict[str, Dict[str, Any]]:
    """Computes Shannon entropy information gain IG(g) = H(T) - E_f[H(T|f)] for every molecule."""
    n = len(records)
    h_prior = math.log2(n) if n > 1 else 0.0
    metrics: Dict[str, Dict[str, Any]] = {}

    for g in records:
        g_id = g["id"]
        sig_counts: Counter[Tuple[Any, ...]] = Counter()
        for t in records:
            t_id = t["id"]
            sim = sim_matrix[g_id][t_id]
            sig = chemle_feedback_signature(g, t, sim)
            sig_counts[sig] += 1

        # Expected posterior entropy E_f[H(T|f)] = sum_f (|T_f|/N) * log2(|T_f|)
        expected_post_h = sum((cnt / n) * math.log2(cnt) for cnt in sig_counts.values() if cnt > 1)
        ig_bits = round(h_prior - expected_post_h, 2)
        expected_remaining = round(sum((cnt * cnt) / n for cnt in sig_counts.values()), 2)

        # Nearest neighbor Tanimoto similarity (excluding self)
        other_sims = sorted(
            (sim_matrix[g_id][other["id"]], other["name"])
            for other in records
            if other["id"] != g_id
        )
        max_sim, nearest_name = other_sims[-1] if other_sims else (0, "None")

        # Difficulty classification based on chemical-space neighborhood & uniqueness
        if max_sim >= 45:
            difficulty = "Hard (Close Analogues)"
        elif max_sim >= 25:
            difficulty = "Balanced"
        else:
            difficulty = "Distinctive"

        metrics[g_id] = {
            "prior_entropy_bits": round(h_prior, 2),
            "information_gain_bits": ig_bits,
            "expected_remaining_candidates": expected_remaining,
            "nearest_neighbor": nearest_name,
            "nearest_neighbor_tanimoto": max_sim,
            "difficulty_tier": difficulty,
        }
    return metrics


def compute_seqle_information_metrics(
    records: List[Dict[str, Any]], aa_props: Dict[str, Dict[str, Any]]
) -> Dict[str, Dict[str, Any]]:
    """Computes Shannon entropy information gain IG(g) for every 5-AA motif in Seqle."""
    n = len(records)
    h_prior = math.log2(n) if n > 1 else 0.0
    metrics: Dict[str, Dict[str, Any]] = {}

    for g in records:
        g_seq = g["seq"]
        sig_counts: Counter[Tuple[str, ...]] = Counter()
        for t in records:
            sig = tuple(evaluate_seqle_guess(g_seq, t["seq"], aa_props))
            sig_counts[sig] += 1

        expected_post_h = sum((cnt / n) * math.log2(cnt) for cnt in sig_counts.values() if cnt > 1)
        ig_bits = round(h_prior - expected_post_h, 2)
        expected_remaining = round(sum((cnt * cnt) / n for cnt in sig_counts.values()), 2)
        unique_residues = len(set(g_seq))

        difficulty = (
            "Hard (Repeat Residues)"
            if unique_residues <= 3
            else ("Balanced" if unique_residues == 4 else "Standard")
        )
        metrics[g_seq] = {
            "information_gain_bits": ig_bits,
            "expected_remaining_candidates": expected_remaining,
            "unique_residues": unique_residues,
            "difficulty_tier": difficulty,
        }
    return metrics
