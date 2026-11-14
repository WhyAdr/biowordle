"""Project-owned RDKit Cheminformatics, Sequence Biophysics & Information-Theory Engine.

Zero external plugin dependencies: uses only public RDKit APIs and the Python
standard library so builds are 100% identical across developer machines, CI, and HPC.
"""

from __future__ import annotations

import hashlib
import itertools
import math
import os
import re
import statistics
import sys
from collections import Counter
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from rdkit import Chem, DataStructs, RDConfig
from rdkit.Chem import (
    Crippen,
    Descriptors,
    Draw,
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
# 1. Explicit Fingerprint, Standardization & Peptide Specifications
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

# ---------------------------------------------------------------------------
# 2. Split Non-Overlapping Functional Group SMARTS Catalog (P1-5)
# ---------------------------------------------------------------------------
FUNCTIONAL_GROUP_SMARTS: List[Tuple[str, str]] = [
    ("Aliphatic Alcohol (-OH)", "[OX2H1][CX4]"),
    ("Phenol / Enol (-OH)", "[OX2H1][c,CX3]"),
    ("Carboxylic Acid (-COOH)", "[CX3](=O)[OX2H1,OX1-]"),
    ("Aldehyde (-CHO)", "[CX3H1](=O)[#6,#1]"),
    ("Ketone (C-C(=O)-C)", "[#6][CX3](=O)[#6]"),
    ("Ester / Lactone", "[#6,#1][CX3](=O)[OX2H0][#6]"),
    ("Ether / Methoxy (C-O-C)", "[OD2]([#6;!$(C=O)])[#6;!$(C=O)]"),
    ("Primary / Secondary Amine", "[NX3;H2,H1;!$(NC=O)]"),
    ("Tertiary Amine", "[NX3;H0;!$(NC=O);!$(N=O)]"),
    ("Amide / Lactam", "[NX3][CX3](=[OX1])[#6,#1]"),
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
    """Detects functional groups present in `mol` using split non-overlapping SMARTS."""
    matched: List[str] = []
    for name, pattern in COMPILED_FG_SMARTS:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            matched.append(name)
    return matched


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
        else (
            (tpsa - 20.0) / 20.0
            if 20.0 < tpsa < 40.0
            else (_ramp(tpsa, 90.0, 120.0) if 90.0 < tpsa < 120.0 else 0.0)
        )
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

    cns_mpo_4prop = round(
        _ramp(clogp, 3.0, 5.0)
        + _ramp(mw, 360.0, 500.0)
        + s_tpsa
        + _ramp(float(hbd), 0.5, 3.5),
        2,
    )
    formula_str = rdMolDescriptors.CalcMolFormula(mol)
    combined_alerts = [f"PAINS: {p}" for p in pains_hits] + [f"Brenk: {b}" for b in brenk_hits]

    return {
        "formula": formula_str,
        "Formula": formula_str,
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
        "CNS_MPO_4prop": cns_mpo_4prop,
        "CNS_MPO_Method": "4-parameter RDKit CNS desirability (cLogP, MW, TPSA, HBD; max 4.0) + 6-term heuristic proxy (max 6.0)",
        "Ro5_Violations": len(ro5_violations),
        "Lipinski_Violations": len(ro5_violations),
        "Ro5_Details": ro5_violations,
        "PAINS_Alerts": pains_hits,
        "Brenk_Alerts": brenk_hits,
        "Alerts": combined_alerts,
    }


def _normalize_svg_floats(svg: str) -> str:
    """Rounds floating-point SVG coordinates to 1 decimal place to suppress platform float noise."""
    return re.sub(r"(-?\d+\.\d)\d+", lambda m: f"{float(m.group(0)):.1f}", svg)


def render_responsive_svg(
    mol: Chem.Mol,
    width: int = 320,
    height: int = 220,
    highlight_atoms: Optional[List[int]] = None,
    dark_mode: bool = False,
) -> str:
    """Renders an authentic ACS-1996 styled responsive SVG string (P1-5)."""
    rdDepictor.SetPreferCoordGen(True)
    mol_copy = Chem.Mol(mol)
    if not mol_copy.GetNumConformers():
        rdDepictor.Compute2DCoords(mol_copy)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    mean_bl = Draw.MeanBondLength(mol_copy)
    if mean_bl > 0:
        rdMolDraw2D.SetACS1996Mode(opts, mean_bl)
    if dark_mode:
        rdMolDraw2D.SetDarkMode(opts)
    opts.clearBackground = not dark_mode
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
    return _normalize_svg_floats(svg.strip())


def build_scaffold_clue_svg(mol: Chem.Mol) -> Tuple[str, str, str]:
    """Builds a progressive intermediate 2D clue SVG without leaking acyclic structures."""
    scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold_mol is not None and scaffold_mol.GetNumHeavyAtoms() >= 3:
        scaf_smi = Chem.MolToSmiles(scaffold_mol, isomericSmiles=False)
        full_flat_smi = Chem.MolToSmiles(mol, isomericSmiles=False)
        mode = "bemis_murcko" if scaf_smi != full_flat_smi else "ring_core"
        return render_responsive_svg(scaffold_mol, width=320, height=220), scaf_smi, mode

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
# 3. Authoritative Spec-Driven Feedback Engine (P1-1 & P0-2)
# ---------------------------------------------------------------------------
def evaluate_delta_from_field_spec(
    guess_val: float, target_val: float, field_spec: Dict[str, Any]
) -> Dict[str, str]:
    """Evaluates a single numeric descriptor comparison strictly from `feedback_spec.json`."""
    diff = abs(float(guess_val) - float(target_val))
    exact_tol = float(field_spec["exact_tol"])
    close_tol = float(field_spec["close_tol"])
    close_dir = bool(field_spec.get("close_directional", False))
    far_dir = bool(field_spec.get("far_directional", True))

    if diff <= exact_tol + 1e-9:
        return {"code": "E", "cls": "match-exact", "arrow": "✓"}

    direction = "U" if float(target_val) > float(guess_val) else "D"
    dir_arrow = "↑" if direction == "U" else "↓"

    if diff <= close_tol + 1e-9:
        if close_dir:
            return {"code": f"C{direction}", "cls": "match-close", "arrow": dir_arrow}
        return {"code": "C", "cls": "match-close", "arrow": "≈"}

    if far_dir:
        return {"code": f"F{direction}", "cls": "match-far", "arrow": dir_arrow}
    return {"code": "F", "cls": "match-far", "arrow": "≠"}


def evaluate_chemle_feedback_from_spec(
    guess: Dict[str, Any],
    target: Dict[str, Any],
    sim_pct: int,
    profile_spec: Dict[str, Any],
) -> Dict[str, Any]:
    """Computes the exact feedback signature and display cells for a Chemle guess."""
    is_win = guess["inchikey"] == target["inchikey"] or guess["id"] == target["id"]
    t_spec = profile_spec["tanimoto_bucket"]
    edges = t_spec["edges"]
    codes = t_spec.get("codes", ["LO", "MD", "HI"])
    labels = t_spec.get("labels", ["Low", "Moderate", "High"])

    if is_win or sim_pct == 100:
        sim_code = "WIN"
        sim_label = "Exact (100%)"
        sim_badge_cls = "badge-emerald"
    elif sim_pct >= edges[1]:
        sim_code = codes[2]
        sim_label = labels[2]
        sim_badge_cls = "badge-emerald"
    elif sim_pct >= edges[0]:
        sim_code = codes[1]
        sim_label = labels[1]
        sim_badge_cls = "badge-amber"
    else:
        sim_code = codes[0]
        sim_label = labels[0]
        sim_badge_cls = "badge-indigo"

    sig_parts: List[str] = []
    if is_win:
        sig_parts.append("WIN")
    else:
        if profile_spec.get("show_category_match_live", False):
            sig_parts.append("CAT_EQ" if guess["Category"] == target["Category"] else "CAT_NE")
        if t_spec.get("shown", True):
            sig_parts.append(sim_code)

    cells: List[Dict[str, Any]] = []
    for f_spec in profile_spec["chemle_fields"]:
        if not f_spec.get("shown", True):
            continue
        col = f_spec["name"]
        res = evaluate_delta_from_field_spec(guess[col], target[col], f_spec)
        if not is_win:
            sig_parts.append(res["code"])
        cells.append(
            {
                "name": col,
                "label": f_spec["label"],
                "decimals": f_spec.get("decimals", 1),
                "value": guess[col],
                **res,
            }
        )

    if not is_win:
        if profile_spec.get("show_exact_tanimoto"):
            sig_parts.append(f"SIM_PCT:{int(sim_pct)}")
        if profile_spec.get("show_functional_groups_live"):
            shared_fg = sorted(
                set(guess.get("functional_groups", []))
                & set(target.get("functional_groups", []))
            )
            sig_parts.append(f"FG:{','.join(shared_fg)}")

    return {
        "signature": "|".join(sig_parts),
        "sim_code": sim_code,
        "sim_label": sim_label,
        "sim_badge_cls": sim_badge_cls,
        "cells": cells,
    }


def _solve_chemle_pool(
    ids: List[str],
    by_id: Dict[str, Dict[str, Any]],
    sig_table: Dict[Tuple[str, str], str],
    profile_spec: Dict[str, Any],
    profile_id: str,
    condition_on_clues: bool,
) -> Dict[str, int]:
    """Runs greedy max-entropy search over Chemle IDs with optional clue conditioning."""
    sc_after = int(profile_spec.get("super_class_unlock_after_guesses", 4 if profile_id == "ranked" else 1))
    scaf_after = int(profile_spec.get("scaffold_unlock_after_guesses", 5 if profile_id == "ranked" else 3))

    def pick_probe(cands: List[str]) -> str:
        if len(cands) == 1:
            return cands[0]
        best_g = cands[0]
        best_key: Optional[Tuple[float, bool, str]] = None
        n_c = len(cands)
        for g_id in ids:
            parts = Counter(sig_table[(g_id, t_id)] for t_id in cands)
            h = -sum((c / n_c) * math.log2(c / n_c) for c in parts.values())
            key = (round(h, 9), g_id in cands, g_id)
            if best_key is None or key > best_key:
                best_g, best_key = g_id, key
        return best_g

    guess_counts: Dict[str, int] = {}
    for t_id in ids:
        target = by_id[t_id]
        ladder = (target.get("clue_ladder") or {}).get(profile_id) or []
        if condition_on_clues and ladder:
            lvl0_text = ladder[0].get("text")
            cands = [
                cid
                for cid in ids
                if (((by_id[cid].get("clue_ladder") or {}).get(profile_id) or [{}])[0].get("text") == lvl0_text)
            ]
        else:
            cands = list(ids)

        steps = 0
        while True:
            g_id = pick_probe(cands)
            steps += 1
            if g_id == t_id or steps >= 12:
                break
            obs_sig = sig_table[(g_id, t_id)]
            cands = [c for c in cands if sig_table[(g_id, c)] == obs_sig]
            if condition_on_clues:
                if steps >= sc_after:
                    cands = [c for c in cands if by_id[c].get("super_class") == target.get("super_class")]
                if steps >= scaf_after:
                    t_scaf = target.get("scaffold_smiles") or ""
                    cands = [c for c in cands if (by_id[c].get("scaffold_smiles") or "") == t_scaf]
        guess_counts[t_id] = steps
    return guess_counts


def _summarize_guess_counts(guess_counts: Dict[str, int]) -> Dict[str, Any]:
    vals = sorted(guess_counts.values())
    dist = dict(sorted(Counter(vals).items()))
    mean_v = round(statistics.mean(vals), 2)
    median_v = round(float(statistics.median(vals)), 2)
    p95_idx = min(len(vals) - 1, int(math.ceil(0.95 * len(vals))) - 1)
    p95_v = vals[p95_idx]
    max_v = max(vals)
    le2_pct = round(100.0 * sum(x <= 2 for x in vals) / len(vals), 2)
    le6_pct = round(100.0 * sum(x <= 6 for x in vals) / len(vals), 2)
    return {
        "distribution": {str(k): v for k, v in dist.items()},
        "mean": mean_v,
        "median": median_v,
        "p95": p95_v,
        "max": max_v,
        "pct_solved_le2": le2_pct,
        "pct_solved_le6": le6_pct,
    }


def run_chemle_greedy_solver(
    records: Sequence[Dict[str, Any]],
    sim_matrix: Dict[str, Dict[str, int]],
    profile_spec: Dict[str, Any],
    profile_id: str = "ranked",
    condition_on_clues: bool = True,
) -> Dict[str, Any]:
    """Runs the greedy max-entropy solver against `profile_spec` (with optional clue conditioning) and returns gate statistics."""
    ids = [r["id"] for r in records]
    by_id = {r["id"]: r for r in records}
    n_total = len(ids)

    sig_table: Dict[Tuple[str, str], str] = {}
    for g_id in ids:
        for t_id in ids:
            sig_table[(g_id, t_id)] = evaluate_chemle_feedback_from_spec(
                by_id[g_id], by_id[t_id], sim_matrix[g_id][t_id], profile_spec
            )["signature"]

    # Opening-probe information metrics for this profile
    h_prior = math.log2(n_total) if n_total > 1 else 0.0
    probe_metrics: Dict[str, Dict[str, float]] = {}
    min_exp_rem = float("inf")
    for g_id in ids:
        cnts = Counter(sig_table[(g_id, t_id)] for t_id in ids)
        exp_post_h = sum((c / n_total) * math.log2(c) for c in cnts.values() if c > 1)
        ig = round(h_prior - exp_post_h, 2)
        exp_rem = round(sum((c * c) / n_total for c in cnts.values()), 2)
        min_exp_rem = min(min_exp_rem, exp_rem)
        probe_metrics[g_id] = {
            "information_gain_bits": ig,
            "expected_remaining_candidates": exp_rem,
        }

    guess_counts = _solve_chemle_pool(
        ids, by_id, sig_table, profile_spec, profile_id, condition_on_clues
    )
    summary = _summarize_guess_counts(guess_counts)
    out: Dict[str, Any] = {
        **summary,
        "min_opening_expected_remaining": round(min_exp_rem, 2),
        "probe_metrics": probe_metrics,
        "per_target_guesses": guess_counts,
        "conditioned_on_clues": condition_on_clues,
    }
    if condition_on_clues:
        raw_counts = _solve_chemle_pool(
            ids, by_id, sig_table, profile_spec, profile_id, condition_on_clues=False
        )
        out["unconditioned_summary"] = _summarize_guess_counts(raw_counts)
    return out


def assert_solver_gate(solver_stats: Dict[str, Any], gate: Dict[str, Any], profile_id: str) -> None:
    """Asserts that `solver_stats` passes all thresholds in `gate` (P0-2)."""
    if not (gate["mean_min"] <= solver_stats["mean"] <= gate["mean_max"]):
        raise AssertionError(
            f"[{profile_id}] Solver mean {solver_stats['mean']} outside gate "
            f"[{gate['mean_min']}, {gate['mean_max']}]"
        )
    if solver_stats["pct_solved_le2"] >= gate["max_le2_pct"]:
        raise AssertionError(
            f"[{profile_id}] Solver <=2-guess rate {solver_stats['pct_solved_le2']}% "
            f"exceeds max {gate['max_le2_pct']}%"
        )
    if solver_stats["pct_solved_le6"] < gate["min_le6_pct"]:
        raise AssertionError(
            f"[{profile_id}] Solver <=6-guess rate {solver_stats['pct_solved_le6']}% "
            f"below min {gate['min_le6_pct']}%"
        )
    if solver_stats["min_opening_expected_remaining"] < gate["min_expected_remaining"]:
        raise AssertionError(
            f"[{profile_id}] Min opening expected remaining {solver_stats['min_opening_expected_remaining']} "
            f"below min {gate['min_expected_remaining']}"
        )


def assert_seqle_solver_gate(solver_stats: Dict[str, Any], gate: Dict[str, Any], profile_id: str) -> None:
    """Asserts that Seqle full-solver stats pass all thresholds in `gate`."""
    assert_solver_gate(solver_stats, gate, f"seqle:{profile_id}")
    max_dist = gate.get("max_single_probe_distinct_targets")
    if max_dist is not None and solver_stats["max_single_probe_distinct_targets"] > int(max_dist):
        raise AssertionError(
            f"[seqle:{profile_id}] Max single-probe distinct targets "
            f"{solver_stats['max_single_probe_distinct_targets']} exceeds max {max_dist}"
        )


# ---------------------------------------------------------------------------
# 4. Seqle Two-Pass Evaluator, Validator & Peptide Biophysics (P0-4)
# ---------------------------------------------------------------------------
def evaluate_seqle_guess(
    guess_seq: str, target_seq: str, aa_props: Dict[str, Dict[str, Any]]
) -> List[str]:
    """Deterministic 3-tier Seqle feedback with strict two-pass duplicate handling."""
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


def evaluate_seqle_feedback_from_spec(
    guess_biochem: Dict[str, Any],
    target_biochem: Dict[str, Any],
    guess_seq: str,
    target_seq: str,
    aa_props: Dict[str, Dict[str, Any]],
    profile_spec: Dict[str, Any],
) -> Dict[str, Any]:
    """Computes full Seqle tile + strip feedback from `feedback_spec.json`."""
    tiles = evaluate_seqle_guess(guess_seq, target_seq, aa_props)
    if not profile_spec.get("show_positional_class_match", True):
        tiles = ["absent" if t == "class" else t for t in tiles]
    strip_cells: List[Dict[str, Any]] = []
    sig_parts = ["".join(t[0].upper() for t in tiles)]

    for f_spec in profile_spec["seqle_strip_fields"]:
        if not f_spec.get("shown", True):
            continue
        col = f_spec["name"]
        res = evaluate_delta_from_field_spec(guess_biochem[col], target_biochem[col], f_spec)
        if guess_seq != target_seq:
            sig_parts.append(res["code"])
        strip_cells.append(
            {
                "name": col,
                "label": f_spec["label"],
                "decimals": f_spec.get("decimals", 2),
                "value": guess_biochem[col],
                **res,
            }
        )
    return {
        "tiles": tiles,
        "strip": strip_cells,
        "signature": "|".join(sig_parts),
    }


def run_seqle_greedy_solver(
    records: Sequence[Dict[str, Any]],
    aa_props: Dict[str, Dict[str, Any]],
    profile_spec: Dict[str, Any],
    profile_id: str = "ranked",
    condition_on_clues: bool = True,
) -> Dict[str, Any]:
    """Runs full multi-step greedy max-entropy solver for Seqle using tile + strip + clue feedback."""
    all_seqs = [r["seq"] for r in records]
    by_seq = {r["seq"]: r for r in records}
    target_seqs = [
        r["seq"] for r in records if profile_id in r.get("game_eligibility", ["learn"])
    ]
    n_t = len(target_seqs)

    sig_table: Dict[Tuple[str, str], str] = {}
    for g_seq in all_seqs:
        g_rec = by_seq[g_seq]
        for t_seq in all_seqs:
            sig_table[(g_seq, t_seq)] = evaluate_seqle_feedback_from_spec(
                g_rec, by_seq[t_seq], g_seq, t_seq, aa_props, profile_spec
            )["signature"]

    min_exp_rem = float("inf")
    max_distinct_targets = 0
    yigsr_distinct = 0
    for g_seq in all_seqs:
        cnts = Counter(sig_table[(g_seq, t_seq)] for t_seq in target_seqs)
        n_dist = len(cnts)
        if n_dist > max_distinct_targets:
            max_distinct_targets = n_dist
        if g_seq == "YIGSR":
            yigsr_distinct = n_dist
        exp_rem = round(sum((c * c) / n_t for c in cnts.values()), 2)
        if exp_rem < min_exp_rem:
            min_exp_rem = exp_rem

    def solve_pool(cond: bool) -> Dict[str, int]:
        def pick_probe(cands: List[str]) -> str:
            if len(cands) == 1:
                return cands[0]
            best_g = cands[0]
            best_key: Optional[Tuple[float, bool, str]] = None
            n_c = len(cands)
            for g_seq in all_seqs:
                parts = Counter(sig_table[(g_seq, t_seq)] for t_seq in cands)
                h = -sum((c / n_c) * math.log2(c / n_c) for c in parts.values())
                key = (round(h, 9), g_seq in cands, g_seq)
                if best_key is None or key > best_key:
                    best_g, best_key = g_seq, key
            return best_g

        counts: Dict[str, int] = {}
        for t_seq in target_seqs:
            target = by_seq[t_seq]
            ladder = (target.get("clue_ladder") or {}).get(profile_id) or []
            if cond and ladder:
                lvl0_text = ladder[0].get("text")
                cands = [
                    s
                    for s in target_seqs
                    if (((by_seq[s].get("clue_ladder") or {}).get(profile_id) or [{}])[0].get("text") == lvl0_text)
                ]
            else:
                cands = list(target_seqs)

            steps = 0
            while True:
                g_seq = pick_probe(cands)
                steps += 1
                if g_seq == t_seq or steps >= 12:
                    break
                obs_sig = sig_table[(g_seq, t_seq)]
                cands = [c for c in cands if sig_table[(g_seq, c)] == obs_sig]
                if cond and len(ladder) > 1:
                    lvl1 = ladder[1]
                    if steps >= int(lvl1.get("after_guesses", 99)):
                        lvl1_text = lvl1.get("text")
                        cands = [
                            c
                            for c in cands
                            if (((by_seq[c].get("clue_ladder") or {}).get(profile_id) or [{}, {}])[1].get("text") == lvl1_text)
                        ]
            counts[t_seq] = steps
        return counts

    guess_counts = solve_pool(condition_on_clues)
    summary = _summarize_guess_counts(guess_counts)
    out: Dict[str, Any] = {
        **summary,
        "eligible_count": n_t,
        "yigsr_distinct_signatures": yigsr_distinct,
        "max_single_probe_distinct_targets": max_distinct_targets,
        "min_opening_expected_remaining": round(min_exp_rem, 2),
        "per_target_guesses": guess_counts,
        "conditioned_on_clues": condition_on_clues,
    }
    if condition_on_clues:
        raw_counts = solve_pool(False)
        out["unconditioned_summary"] = _summarize_guess_counts(raw_counts)
    return out


def _consensus_tokens(consensus: str) -> List[str]:
    """Splits a PROSITE-style consensus like 'G-X-S-X-G' or '[V/I]-D-F-G-[L/M]' into tokens."""
    tokens: List[str] = []
    for part in consensus.split("-"):
        part = part.strip()
        if not part:
            continue
        m = re.match(r"^X\((\d+)(?:,(\d+))?\)$", part)
        if m:
            count = int(m.group(1))
            tokens.extend(["[A-Z]"] * count)
            continue
        m2 = re.match(r"^X(\d+)$", part)
        if m2:
            tokens.extend(["[A-Z]"] * int(m2.group(1)))
            continue
        if part == "X":
            tokens.append("[A-Z]")
        elif part.startswith("[") and part.endswith("]"):
            inner = part[1:-1].replace("/", "")
            tokens.append(f"[{inner}]")
        else:
            tokens.append(re.escape(part))
    return tokens


# Authoritative reference provenance registry verified against UniProtKB/Swiss-Prot
# FASTA and RCSB PDB polymer entities (BioWordle v0.2.1).
SEQLE_REFERENCE_REGISTRY: Dict[str, Dict[str, Any]] = {
    "GHSMG": {
        "status": "verified_uniprot",
        "accession": "G9BY57",
        "isoform": "canonical",
        "parent_residue_range": [160, 170],
        "motif_residue_range": [163, 167],
        "window_offset": 3,
        "verified_parent_sequence": "AVAGHSMGGGG",
        "pdb_ids": ["4EB0"],
        "notes": "Leaf-branch compost cutinase (LCC); previously mis-cited PETase A0A0K8P6T7 carries GWSMG.",
    },
    "EKDEL": {
        "status": "verified_uniprot",
        "accession": "P11021",
        "isoform": "canonical",
        "parent_residue_range": [648, 654],
        "motif_residue_range": [650, 654],
        "window_offset": 2,
        "verified_parent_sequence": "TAEKDEL",
        "pdb_ids": ["6I6H"],
        "notes": "Human HSPA5/BiP extreme C-terminal ER retention signal.",
    },
    "DEVDG": {
        "status": "verified_uniprot",
        "accession": "P09874",
        "isoform": "canonical",
        "parent_residue_range": [210, 219],
        "motif_residue_range": [211, 215],
        "window_offset": 1,
        "verified_parent_sequence": "GDEVDGVDEV",
        "pdb_ids": [],
        "notes": "Human PARP1 Caspase-3/7 cleavage site.",
    },
    "GRGDS": {
        "status": "verified_uniprot",
        "accession": "P02751",
        "isoform": "canonical",
        "parent_residue_range": [1612, 1622],
        "motif_residue_range": [1614, 1618],
        "window_offset": 2,
        "verified_parent_sequence": "VTGRGDSPASS",
        "pdb_ids": ["1FNF"],
        "notes": "Human Fibronectin FN1 canonical coordinates 1614-1618 (older mature fragment numbering 1492-1496).",
    },
    "LPETG": {
        "status": "verified_uniprot",
        "accession": "P02976",
        "isoform": "canonical",
        "parent_residue_range": [480, 489],
        "motif_residue_range": [482, 486],
        "window_offset": 2,
        "verified_parent_sequence": "QALPETGEEN",
        "pdb_ids": ["1T2W"],
        "notes": "Staphylococcus aureus Protein A Sortase-A cell-wall sorting signal.",
    },
    "HEIGH": {
        "status": "verified_uniprot",
        "accession": "P08254",
        "isoform": "canonical",
        "parent_residue_range": [215, 226],
        "motif_residue_range": [218, 222],
        "window_offset": 3,
        "verified_parent_sequence": "VAAHEIGHSLGL",
        "pdb_ids": ["1B3D"],
        "notes": "Human MMP3 Stromelysin-1 Zn-binding helix; previously mis-cited MMP1 P03956 carries HELGH.",
    },
    "CCHHC": {
        "status": "abstraction",
        "accession": "PROSITE:PS50158",
        "isoform": None,
        "parent_residue_range": None,
        "motif_residue_range": None,
        "window_offset": 0,
        "verified_parent_sequence": "CCHHC",
        "pdb_ids": ["1MFS"],
        "notes": "Learn-only pedagogical compression of non-contiguous C-X2-C-X4-H-X4-C zinc knuckle.",
    },
    "YIGSR": {
        "status": "verified_uniprot",
        "accession": "P07942",
        "isoform": "canonical",
        "parent_residue_range": [946, 956],
        "motif_residue_range": [950, 954],
        "window_offset": 4,
        "verified_parent_sequence": "CDPGYIGSRCD",
        "pdb_ids": [],
        "notes": "Human LAMB1 canonical coordinates 950-954 (mature chain omitting 21-aa signal peptide: 929-933).",
    },
    "NPVYK": {
        "status": "consensus_exemplar",
        "accession": "ELM:TRG_ENDOCYTIC_2",
        "isoform": None,
        "parent_residue_range": None,
        "motif_residue_range": None,
        "window_offset": 0,
        "verified_parent_sequence": "NPVYK",
        "pdb_ids": ["3SO6"],
        "notes": "NPXY[KR] consensus representative; native human LDLR P01130 (825-829) and PDB:3SO6 carry NPVYQ.",
    },
    "KCVIM": {
        "status": "verified_uniprot",
        "accession": "P01116-2",
        "isoform": "2B",
        "parent_residue_range": [175, 188],
        "motif_residue_range": [184, 188],
        "window_offset": 9,
        "verified_parent_sequence": "KKKKKKSKTKCVIM",
        "pdb_ids": ["1D8D"],
        "notes": "Human KRAS Isoform 2B polybasic C-terminal CaaX farnesylation box.",
    },
    "DEADE": {
        "status": "verified_uniprot",
        "accession": "P60842",
        "isoform": "canonical",
        "parent_residue_range": [180, 190],
        "motif_residue_range": [182, 186],
        "window_offset": 2,
        "verified_parent_sequence": "VLDEADEMLSR",
        "pdb_ids": ["2ZU6"],
        "notes": "Human eIF4A1 / DDX2A Walker-B DEAD-box helicase motif.",
    },
    "GVGKS": {
        "status": "verified_uniprot",
        "accession": "P01116",
        "isoform": "canonical",
        "parent_residue_range": [10, 17],
        "motif_residue_range": [13, 17],
        "window_offset": 3,
        "verified_parent_sequence": "GAGGVGKS",
        "pdb_ids": ["4OBE"],
        "notes": "Human KRAS GTPase Walker-A P-loop.",
    },
    "VDFGL": {
        "status": "verified_uniprot",
        "accession": "P41279",
        "isoform": "canonical",
        "parent_residue_range": [265, 277],
        "motif_residue_range": [269, 273],
        "window_offset": 4,
        "verified_parent_sequence": "KAVLVDFGLSVQM",
        "pdb_ids": ["4Y83"],
        "notes": "Human MAP3K8/COT activation-loop DFG switch; previously mis-cited ABL1 P00519 carries KVADFGLSR.",
    },
    "HRDLA": {
        "status": "verified_uniprot",
        "accession": "P10721",
        "isoform": "canonical",
        "parent_residue_range": [789, 797],
        "motif_residue_range": [790, 794],
        "window_offset": 1,
        "verified_parent_sequence": "IHRDLAARN",
        "pdb_ids": ["1PKG"],
        "notes": "Human KIT catalytic loop (also FGFR1 P11362 620-628); previously mis-cited SRC P12931 carries HRDLR.",
    },
    "WSEWS": {
        "status": "verified_uniprot",
        "accession": "P32927",
        "isoform": "canonical",
        "parent_residue_range": [421, 433],
        "motif_residue_range": [425, 429],
        "window_offset": 4,
        "verified_parent_sequence": "YNGIWSEWSEARS",
        "pdb_ids": ["1EGJ"],
        "notes": "Human CSF2RB cytokine receptor WSXWS box; canonical EPOR P19235 carries WSAWS (PDB:1EER has A235E variant).",
    },
    "PPPVP": {
        "status": "verified_uniprot",
        "accession": "Q07889",
        "isoform": "canonical",
        "parent_residue_range": [1149, 1158],
        "motif_residue_range": [1150, 1154],
        "window_offset": 1,
        "verified_parent_sequence": "VPPPVPPRRR",
        "pdb_ids": ["1GBR"],
        "notes": "Human SOS1 Class-II SH3-binding polyproline motif.",
    },
    "INFEK": {
        "status": "verified_uniprot",
        "accession": "P01012",
        "isoform": "canonical",
        "parent_residue_range": [258, 265],
        "motif_residue_range": [260, 264],
        "window_offset": 2,
        "verified_parent_sequence": "SIINFEKL",
        "pdb_ids": ["1VAC"],
        "notes": "Gallus gallus Ovalbumin H-2Kb epitope SIINFEKL (canonical coords 258-265).",
    },
    "YQRLN": {
        "status": "verified_uniprot",
        "accession": "P19814",
        "isoform": "canonical",
        "parent_residue_range": [348, 357],
        "motif_residue_range": [350, 354],
        "window_offset": 2,
        "verified_parent_sequence": "SDYQRLNLKL",
        "pdb_ids": ["1BXX"],
        "notes": "Rat TGN38 tyrosine-based AP2 internalization motif.",
    },
    "KKMPQ": {
        "status": "consensus_exemplar",
        "accession": "ELM:TRG_ER_diLys_1",
        "isoform": None,
        "parent_residue_range": None,
        "motif_residue_range": None,
        "window_offset": 0,
        "verified_parent_sequence": "KKMPQ",
        "pdb_ids": ["4J77"],
        "notes": "Synthetic C-terminal -KKxxx COPI retrieval consensus exemplar; PDB:4J77 co-crystal peptide is AKEKSD.",
    },
    "CAACH": {
        "status": "verified_uniprot",
        "accession": "P07143",
        "isoform": "canonical",
        "parent_residue_range": [97, 109],
        "motif_residue_range": [101, 105],
        "window_offset": 4,
        "verified_parent_sequence": "YREVCAACHSLDR",
        "pdb_ids": ["3CX5"],
        "notes": "Yeast mitochondrial Cytochrome c1 CXXCH heme site; previously mis-cited CYC1 P00044 carries CLQCH.",
    },
    "NLTGY": {
        "status": "consensus_exemplar",
        "accession": "ELM:MOD_N-GLC_1",
        "isoform": None,
        "parent_residue_range": None,
        "motif_residue_range": None,
        "window_offset": 0,
        "verified_parent_sequence": "NLTGY",
        "pdb_ids": ["3WAK"],
        "notes": "N-X-[ST] glycosylation sequon consensus exemplar; PDB:3WAK is apo Archaeoglobus fulgidus AglB (O29867).",
    },
    "DSGLD": {
        "status": "verified_uniprot",
        "accession": "P25963",
        "isoform": "canonical",
        "parent_residue_range": [27, 39],
        "motif_residue_range": [31, 35],
        "window_offset": 4,
        "verified_parent_sequence": "DDRHDSGLDSMKD",
        "pdb_ids": ["1P22"],
        "notes": "Human NFKBIA / I-kappa-B-alpha beta-TrCP phosphodegron.",
    },
    "KENIP": {
        "status": "consensus_exemplar",
        "accession": "ELM:DEG_APCC_KENBOX_2",
        "isoform": None,
        "parent_residue_range": None,
        "motif_residue_range": None,
        "window_offset": 0,
        "verified_parent_sequence": "KENIP",
        "pdb_ids": ["4BH6"],
        "notes": "APC/C KEN-X-P degron consensus exemplar; PDB:4BH6 Acm1 (Q08981) carries KENLS and Hsl1 (P34244) carries NKENE.",
    },
    "KKKRK": {
        "status": "verified_uniprot",
        "accession": "P03070",
        "isoform": "canonical",
        "parent_residue_range": [126, 132],
        "motif_residue_range": [127, 131],
        "window_offset": 1,
        "verified_parent_sequence": "PKKKRKV",
        "pdb_ids": ["1BK6"],
        "notes": "SV40 Large T-antigen monopartite nuclear localization signal.",
    },
    "FGGFT": {
        "status": "verified_uniprot",
        "accession": "Q13519",
        "isoform": "canonical",
        "parent_residue_range": [130, 146],
        "motif_residue_range": [130, 134],
        "window_offset": 0,
        "verified_parent_sequence": "FGGFTGARKSARKLANQ",
        "pdb_ids": [],
        "notes": "Human Prepronociceptin cleaved mature Nociceptin (130-146); N-terminus generated by dibasic KR cleavage.",
    },
    "YGGFM": {
        "status": "verified_uniprot",
        "accession": "P01210",
        "isoform": "canonical",
        "parent_residue_range": [100, 104],
        "motif_residue_range": [100, 104],
        "window_offset": 0,
        "verified_parent_sequence": "YGGFM",
        "pdb_ids": ["1PLW"],
        "notes": "Human Proenkephalin-A cleaved mature Met-enkephalin pentapeptide (100-104).",
    },
    "SYSGY": {
        "status": "verified_uniprot",
        "accession": "P35637",
        "isoform": "canonical",
        "parent_residue_range": [34, 43],
        "motif_residue_range": [37, 41],
        "window_offset": 3,
        "verified_parent_sequence": "GQQSYSGYSQ",
        "pdb_ids": ["5W3N"],
        "notes": "Human FUS low-complexity domain canonical repeat 37-41 (previously noted as 139-143).",
    },
    "LDLAS": {
        "status": "abstraction",
        "accession": "ELM:TRG_NES_CRM1_1",
        "isoform": None,
        "parent_residue_range": None,
        "motif_residue_range": None,
        "window_offset": 0,
        "verified_parent_sequence": "LDLAS",
        "pdb_ids": ["3NBY"],
        "notes": "Learn-only pedagogical compression of amphipathic CRM1 nuclear export ridge.",
    },
    "NPLYK": {
        "status": "verified_uniprot",
        "accession": "P05106",
        "isoform": "canonical",
        "parent_residue_range": [763, 788],
        "motif_residue_range": [770, 774],
        "window_offset": 7,
        "verified_parent_sequence": "AKWDTANNPLYKEATSTFTNITYRGT",
        "pdb_ids": ["1MK7"],
        "notes": "Human Integrin Beta-3 (ITGB3) canonical cytoplasmic talin-binding turn 770-774 (744-748 mature).",
    },
    "GHHPH": {
        "status": "verified_uniprot",
        "accession": "P04196",
        "isoform": "canonical",
        "parent_residue_range": [383, 392],
        "motif_residue_range": [383, 387],
        "window_offset": 0,
        "verified_parent_sequence": "GHHPHGHHPH",
        "pdb_ids": [],
        "notes": "Human HRG tandem HPR repeat 383-392 (motif 383-387; single-copy repeat also at 363-367).",
    },
}


def validate_seqle_motif_record(
    item: Dict[str, Any],
    aa_props: Dict[str, Dict[str, Any]],
    enforce_reference_registry: bool = True,
) -> None:
    """Validates a Seqle motif record against P0-4 schema, biological constraints, and reference provenance."""
    seq = item.get("seq", "")
    if len(seq) != 5 or any(ch not in aa_props for ch in seq):
        raise ValueError(f"Invalid 5-AA sequence: {seq!r}")

    kind = item.get("kind")
    if kind not in ("exact", "fragment", "abstraction"):
        raise ValueError(f"[{seq}] Invalid kind {kind!r}; expected exact|fragment|abstraction")

    terminal = item.get("terminal")
    if terminal not in ("N", "C", "none"):
        raise ValueError(f"[{seq}] Invalid terminal {terminal!r}; expected N|C|none")

    exemplar = item.get("exemplar_source", "")
    if not re.search(r"(UniProt|PDB|PROSITE|ELM):[A-Za-z0-9_-]+", exemplar):
        raise ValueError(f"[{seq}] exemplar_source {exemplar!r} must cite a UniProt/PDB/PROSITE/ELM accession")

    parent_seq = item.get("parent_sequence", "")
    offset = item.get("window_offset", 0)
    if not isinstance(parent_seq, str) or len(parent_seq) < 5:
        raise ValueError(f"[{seq}] parent_sequence must be at least 5 residues")
    if not isinstance(offset, int) or offset < 0 or offset + 5 > len(parent_seq):
        raise ValueError(f"[{seq}] window_offset {offset} out of bounds for parent_sequence {parent_seq!r}")
    if parent_seq[offset : offset + 5] != seq:
        raise ValueError(
            f"[{seq}] parent_sequence[{offset}:{offset+5}] = {parent_seq[offset:offset+5]!r} != {seq!r}"
        )

    if terminal == "C" and offset + 5 != len(parent_seq):
        raise ValueError(
            f"[{seq}] terminal='C' requires 5-mer window at the extreme C-terminus of parent_sequence {parent_seq!r}"
        )
    if terminal == "N" and offset != 0:
        raise ValueError(
            f"[{seq}] terminal='N' requires 5-mer window at offset 0 of parent_sequence {parent_seq!r}"
        )

    consensus = item.get("consensus", "")
    tokens = _consensus_tokens(consensus)
    if kind == "exact":
        if len(tokens) != 5:
            raise ValueError(
                f"[{seq}] kind='exact' requires a 5-position consensus, got {len(tokens)} in {consensus!r}"
            )
        regex = "^" + "".join(tokens) + "$"
        if not re.match(regex, seq):
            raise ValueError(f"[{seq}] Sequence does not match its exact consensus {consensus!r} ({regex})")
    elif kind == "fragment":
        if len(tokens) < 5:
            raise ValueError(f"[{seq}] kind='fragment' consensus {consensus!r} shorter than 5 positions")

    elig = item.get("game_eligibility", [])
    if not elig or not set(elig).issubset({"learn", "ranked"}):
        raise ValueError(f"[{seq}] Invalid game_eligibility {elig!r}")
    if kind == "abstraction" and "ranked" in elig:
        raise ValueError(f"[{seq}] kind='abstraction' cannot be eligible for 'ranked' mode")

    if enforce_reference_registry and seq in SEQLE_REFERENCE_REGISTRY:
        ref = SEQLE_REFERENCE_REGISTRY[seq]
        if parent_seq != ref["verified_parent_sequence"]:
            raise ValueError(
                f"[{seq}] parent_sequence {parent_seq!r} does not match verified reference "
                f"{ref['verified_parent_sequence']!r} ({ref['accession']})"
            )
        if offset != ref["window_offset"]:
            raise ValueError(
                f"[{seq}] window_offset {offset} != verified offset {ref['window_offset']} ({ref['accession']})"
            )
        status = ref["status"]
        if status == "verified_uniprot":
            acc = str(ref["accession"])
            if f"UniProt:{acc}" not in exemplar:
                raise ValueError(
                    f"[{seq}] exemplar_source {exemplar!r} must cite verified UniProt:{acc}"
                )
            p_rng = ref["parent_residue_range"]
            m_rng = ref["motif_residue_range"]
            if not (isinstance(p_rng, list) and len(p_rng) == 2 and p_rng[1] - p_rng[0] + 1 == len(parent_seq)):
                raise ValueError(f"[{seq}] Invalid parent_residue_range {p_rng!r} for len={len(parent_seq)}")
            expected_motif_start = p_rng[0] + offset
            if m_rng != [expected_motif_start, expected_motif_start + 4]:
                raise ValueError(
                    f"[{seq}] motif_residue_range {m_rng!r} inconsistent with parent range {p_rng!r} + offset {offset}"
                )
        elif status in ("consensus_exemplar", "abstraction"):
            if parent_seq != seq or offset != 0:
                raise ValueError(
                    f"[{seq}] {status} must not fabricate flanking residues (got parent={parent_seq!r}, offset={offset})"
                )


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


# ---------------------------------------------------------------------------
# 5. Opening-Clue k-Anonymity Selector & Hint Leakage Audit (P0-1 & P2-3)
# ---------------------------------------------------------------------------
def select_k_anonymous_level0_features(
    records: Sequence[Dict[str, Any]],
    candidate_extractors: List[Tuple[str, Callable[[Dict[str, Any]], str]]],
    min_k: int,
    max_features: int = 2,
) -> Tuple[List[str], Dict[int, str], int, List[int]]:
    """Searches feature subsets of size 1..max_features to maximize group count subject to k_min >= min_k."""
    best_choice: Optional[Tuple[int, int, List[str], Dict[int, str], List[int]]] = None

    for r_size in range(1, min(max_features, len(candidate_extractors)) + 1):
        for combo in itertools.combinations(candidate_extractors, r_size):
            names = [c[0] for c in combo]
            funcs = [c[1] for c in combo]
            labels_by_idx: Dict[int, str] = {}
            counts: Counter[str] = Counter()
            for idx, rec in enumerate(records):
                val = " · ".join(fn(rec) for fn in funcs)
                labels_by_idx[idx] = val
                counts[val] += 1
            k_min = min(counts.values()) if counts else 0
            if k_min >= min_k:
                num_groups = len(counts)
                score = (num_groups, k_min)
                if best_choice is None or score > (best_choice[0], best_choice[1]):
                    best_choice = (
                        num_groups,
                        k_min,
                        names,
                        labels_by_idx,
                        sorted(counts.values()),
                    )

    if best_choice is None:
        raise AssertionError(
            f"No Level-0 opening clue feature subset satisfies k_min >= {min_k} across {len(records)} records!"
        )
    return best_choice[2], best_choice[3], best_choice[1], best_choice[4]


def audit_hint_text_leakage(record: Dict[str, Any], text_fields: Sequence[str], id_key: str) -> None:
    """Asserts that no pre-reveal hint text contains the target's name, ID, sequence, or >=5-char stem (P2-3)."""
    generic_structural_words = {
        "acid",
        "acids",
    }
    forbidden_tokens: set[str] = set()
    raw_id = str(record.get(id_key, "")).lower()
    if len(raw_id) >= 4:
        forbidden_tokens.add(raw_id)

    if id_key == "seq":
        seq_up = str(record.get("seq", "")).upper()
        if len(seq_up) == 5:
            forbidden_tokens.add(seq_up.lower())
            forbidden_tokens.add(seq_up[:4].lower())
            forbidden_tokens.add(seq_up[1:].lower())
            aa3 = {
                "A": "ala", "C": "cys", "D": "asp", "E": "glu", "F": "phe",
                "G": "gly", "H": "his", "I": "ile", "K": "lys", "L": "leu",
                "M": "met", "N": "asn", "P": "pro", "Q": "gln", "R": "arg",
                "S": "ser", "T": "thr", "V": "val", "W": "trp", "Y": "tyr",
            }
            for i in range(3):
                tri = f"{aa3[seq_up[i]]}-{aa3[seq_up[i+1]]}-{aa3[seq_up[i+2]]}"
                forbidden_tokens.add(tri)
    else:
        for part in re.findall(r"[a-z0-9]{5,}", raw_id):
            if part not in generic_structural_words:
                forbidden_tokens.add(part)
                forbidden_tokens.add(part[:6] if len(part) >= 6 else part[:5])

        name = str(record.get("name") or "").lower()
        for word in re.findall(r"[a-z0-9]{5,}", name):
            if word not in generic_structural_words:
                forbidden_tokens.add(word)
                if len(word) >= 6:
                    forbidden_tokens.add(word[:6])

    for field in text_fields:
        val = str(record.get(field, "")).lower()
        for tok in forbidden_tokens:
            if re.search(rf"\b{re.escape(tok)}", val):
                raise AssertionError(
                    f"[P2-3 Hint Leakage] Record {record.get(id_key)!r} leaks token/stem {tok!r} in field {field!r}: {record.get(field)!r}"
                )
