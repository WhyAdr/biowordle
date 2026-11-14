"""Project-owned RDKit Cheminformatics, Sequence Biophysics & Information-Theory Engine.

Zero external plugin dependencies: uses only public RDKit APIs and the Python
standard library so builds are 100% identical across developer machines, CI, and HPC.
"""

from __future__ import annotations

import hashlib
import itertools
import math
import os
import random
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
    ("Phenol / Enol (-OH)", "[OX2H1][$([c]),$([CX3]=[CX3])]"),
    ("Carboxylic Acid (-COOH)", "[CX3;H1,$([CX3]-[#6])](=O)[OX2H1,OX1-]"),
    ("Aldehyde (-CHO)", "[$([CX3H1](=O)[#6]),$([CX3H2]=O)]"),
    ("Ketone (C-C(=O)-C)", "[#6][CX3](=O)[#6]"),
    ("Ester / Lactone", "[CX3;H1,$([CX3]-[#6])](=O)[OX2H0][#6]"),
    ("Ether / Methoxy (C-O-C)", "[OD2]([#6;!$(C=O)])[#6;!$(C=O)]"),
    ("Primary / Secondary Amine", "[NX3;H2,H1;!$(N[C,S,P]=[O,S,N]);!$(N=*)]"),
    ("Tertiary Amine", "[NX3;H0;!$(N[C,S,P]=[O,S,N]);!$(N=*)]"),
    ("Amide / Lactam", "[NX3][CX3]=[OX1]"),
    ("Indole / Pyrrole Ring", "[nH1]1cccc1"),
    ("Heteroaromatic N", "[n]"),
    ("Halogen (F/Cl/Br/I)", "[F,Cl,Br,I]"),
    ("Sulfur-containing group", "[#16]"),
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


def _ladder_clue_text(record: Dict[str, Any], profile_id: str, level: int) -> Optional[str]:
    ladder = (record.get("clue_ladder") or {}).get(profile_id) or []
    for entry in ladder:
        if entry.get("level") == level:
            return entry.get("text")
    return None


def _solve_chemle_pool(
    answer_ids: Sequence[str],
    probe_ids: Sequence[str],
    by_id: Dict[str, Dict[str, Any]],
    sig_table: Dict[Tuple[str, str], str],
    profile_spec: Dict[str, Any],
    profile_id: str,
    condition_on_clues: bool,
    max_steps: Optional[int] = None,
) -> Dict[str, Tuple[int, bool]]:
    """Runs greedy max-entropy search over `answer_ids` using `probe_ids` with clue-ladder conditioning."""
    sc_after = int(profile_spec.get("super_class_unlock_after_guesses", 4 if profile_id == "ranked" else 1))
    scaf_after = int(profile_spec.get("scaffold_unlock_after_guesses", 5 if profile_id == "ranked" else 3))
    step_cap = int(max_steps if max_steps is not None else profile_spec.get("max_guesses", 12))

    def pick_probe(cands: List[str]) -> str:
        if len(cands) == 1:
            return cands[0]
        cand_set = set(cands)
        best_g = cands[0]
        best_key: Optional[Tuple[float, bool, str]] = None
        n_c = len(cands)
        for g_id in probe_ids:
            parts = Counter(sig_table[(g_id, t_id)] for t_id in cands)
            h = -sum((c / n_c) * math.log2(c / n_c) for c in parts.values())
            key = (round(h, 9), g_id in cand_set, g_id)
            if best_key is None or key > best_key:
                best_g, best_key = g_id, key
        return best_g

    results: Dict[str, Tuple[int, bool]] = {}
    for t_id in answer_ids:
        target = by_id[t_id]
        ladder = sorted((target.get("clue_ladder") or {}).get(profile_id) or [], key=lambda e: int(e.get("level", 0)))
        cands = list(answer_ids)
        if condition_on_clues:
            for clue in ladder:
                if int(clue.get("after_guesses", 99)) <= 0:
                    lvl = int(clue["level"])
                    txt = clue.get("text")
                    cands = [cid for cid in cands if _ladder_clue_text(by_id[cid], profile_id, lvl) == txt]

        steps = 0
        solved = False
        while steps < step_cap:
            g_id = pick_probe(cands)
            steps += 1
            if g_id == t_id:
                solved = True
                break
            obs_sig = sig_table[(g_id, t_id)]
            cands = [c for c in cands if sig_table[(g_id, c)] == obs_sig]
            if condition_on_clues:
                if steps >= sc_after:
                    cands = [c for c in cands if by_id[c].get("super_class") == target.get("super_class")]
                if steps >= scaf_after:
                    t_scaf = target.get("scaffold_smiles") or ""
                    cands = [c for c in cands if (by_id[c].get("scaffold_smiles") or "") == t_scaf]
                for clue in ladder:
                    if int(clue.get("after_guesses", 99)) <= steps:
                        lvl = int(clue["level"])
                        txt = clue.get("text")
                        cands = [c for c in cands if _ladder_clue_text(by_id[c], profile_id, lvl) == txt]
        results[t_id] = (steps, solved)
    return results


def _summarize_guess_counts(
    guess_results: Dict[str, Any],
    mode_cap: int = 6,
) -> Dict[str, Any]:
    """Summarizes solve step counts while keeping misses in an explicit `unsolved` bucket."""
    solved_vals: List[int] = []
    unsolved_targets: List[str] = []
    for key, val in guess_results.items():
        if isinstance(val, tuple):
            steps, is_solved = int(val[0]), bool(val[1])
        elif val is None:
            steps, is_solved = 0, False
        else:
            steps, is_solved = int(val), True
        if is_solved:
            solved_vals.append(steps)
        else:
            unsolved_targets.append(key)

    solved_vals.sort()
    total_count = len(guess_results)
    solved_count = len(solved_vals)
    unsolved_count = len(unsolved_targets)
    dist: Dict[str, int] = {str(k): v for k, v in sorted(Counter(solved_vals).items())}
    if unsolved_count > 0:
        dist["unsolved"] = unsolved_count

    if solved_vals:
        mean_v: Optional[float] = round(statistics.mean(solved_vals), 2)
        median_v: Optional[float] = round(float(statistics.median(solved_vals)), 2)
        p95_idx = min(len(solved_vals) - 1, int(math.ceil(0.95 * len(solved_vals))) - 1)
        p95_v: Optional[int] = solved_vals[p95_idx]
        max_v: Optional[int] = max(solved_vals)
    else:
        mean_v = None
        median_v = None
        p95_v = None
        max_v = None

    le2_pct = round(100.0 * sum(x <= 2 for x in solved_vals) / total_count, 2) if total_count else 0.0
    le6_pct = round(100.0 * sum(x <= 6 for x in solved_vals) / total_count, 2) if total_count else 0.0
    le8_pct = round(100.0 * sum(x <= 8 for x in solved_vals) / total_count, 2) if total_count else 0.0
    cap_pct = round(100.0 * sum(x <= mode_cap for x in solved_vals) / total_count, 2) if total_count else 0.0
    return {
        "distribution": dist,
        "target_count": total_count,
        "solved_count": solved_count,
        "unsolved_count": unsolved_count,
        "unsolved_targets": sorted(unsolved_targets),
        "mode_cap": int(mode_cap),
        "mean": mean_v,
        "median": median_v,
        "p95": p95_v,
        "max": max_v,
        "pct_solved_le2": le2_pct,
        "pct_solved_le6": le6_pct,
        "pct_solved_le8": le8_pct,
        "pct_solved_by_cap": cap_pct,
    }


def run_chemle_greedy_solver(
    records: Sequence[Dict[str, Any]],
    sim_matrix: Dict[str, Dict[str, int]],
    profile_spec: Dict[str, Any],
    profile_id: str = "ranked",
    condition_on_clues: bool = True,
    max_steps: Optional[int] = None,
    answer_ids: Optional[Sequence[str]] = None,
    probe_ids: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """Runs the greedy max-entropy solver over `answer_ids` and `probe_ids` under `profile_spec`."""
    by_id = {r["id"]: r for r in records}
    resolved_answers = (
        list(answer_ids)
        if answer_ids is not None
        else [
            r["id"]
            for r in records
            if profile_id in r.get("game_eligibility", ["learn", "ranked"])
        ]
    )
    resolved_probes = (
        list(probe_ids)
        if probe_ids is not None
        else [
            r["id"]
            for r in records
            if profile_id in r.get("guess_eligibility", r.get("game_eligibility", ["learn", "ranked"]))
        ]
    )
    if not resolved_answers:
        raise ValueError(f"[{profile_id}] Chemle solver answer_ids cannot be empty")
    if not set(resolved_answers).issubset(set(resolved_probes)):
        missing = sorted(set(resolved_answers) - set(resolved_probes))
        raise ValueError(f"[{profile_id}] Eligible answers must also be legal probes; missing={missing}")

    n_ans = len(resolved_answers)
    mode_cap = int(profile_spec.get("max_guesses", 6 if profile_id == "ranked" else 8))
    step_cap = int(max_steps if max_steps is not None else mode_cap)

    sig_table: Dict[Tuple[str, str], str] = {}
    for g_id in resolved_probes:
        for t_id in resolved_answers:
            sig_table[(g_id, t_id)] = evaluate_chemle_feedback_from_spec(
                by_id[g_id], by_id[t_id], sim_matrix[g_id][t_id], profile_spec
            )["signature"]

    # Opening-probe information metrics over the active answer universe (diagnostic)
    h_prior = math.log2(n_ans) if n_ans > 1 else 0.0
    probe_metrics: Dict[str, Dict[str, float]] = {}
    min_exp_rem = float("inf")
    for g_id in resolved_probes:
        cnts = Counter(sig_table[(g_id, t_id)] for t_id in resolved_answers)
        exp_post_h = sum((c / n_ans) * math.log2(c) for c in cnts.values() if c > 1)
        ig = round(h_prior - exp_post_h, 2)
        exp_rem = round(sum((c * c) / n_ans for c in cnts.values()), 2)
        min_exp_rem = min(min_exp_rem, exp_rem)
        probe_metrics[g_id] = {
            "information_gain_bits": ig,
            "expected_remaining_candidates": exp_rem,
        }

    guess_results = _solve_chemle_pool(
        resolved_answers, resolved_probes, by_id, sig_table, profile_spec, profile_id, condition_on_clues, max_steps=step_cap
    )
    summary = _summarize_guess_counts(guess_results, mode_cap=mode_cap)
    out: Dict[str, Any] = {
        **summary,
        "eligible_count": n_ans,
        "probe_count": len(resolved_probes),
        "min_opening_expected_remaining": round(min_exp_rem, 2),
        "probe_metrics": probe_metrics,
        "per_target_guesses": {k: (v[0] if v[1] else None) for k, v in guess_results.items()},
        "per_target_results": guess_results,
        "conditioned_on_clues": condition_on_clues,
    }
    if condition_on_clues:
        raw_results = _solve_chemle_pool(
            resolved_answers, resolved_probes, by_id, sig_table, profile_spec, profile_id, condition_on_clues=False, max_steps=12
        )
        out["unconditioned_summary"] = _summarize_guess_counts(raw_results, mode_cap=mode_cap)
    return out


def assert_solver_gate(solver_stats: Dict[str, Any], gate: Dict[str, Any], profile_id: str) -> None:
    """Asserts that `solver_stats` passes all thresholds in `gate` (P0-2)."""
    if solver_stats.get("mean") is None:
        raise AssertionError(f"[{profile_id}] Zero solved targets!")
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
    min_cap_pct = gate.get("min_solved_by_cap_pct")
    if min_cap_pct is not None and solver_stats["pct_solved_by_cap"] < min_cap_pct:
        raise AssertionError(
            f"[{profile_id}] Solver <={solver_stats['mode_cap']}-guess rate {solver_stats['pct_solved_by_cap']}% "
            f"below min {min_cap_pct}%"
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


def generate_legal_seqle_probes(
    aa_props: Dict[str, Dict[str, Any]],
    seed: int = 20261114,
    count: int = 10000,
) -> List[str]:
    """Generates a deterministic sample of legal 5-AA non-catalog probe sequences."""
    alphabet = sorted(aa_props.keys())
    rng = random.Random(seed)
    probes: set[str] = set()
    while len(probes) < count:
        probes.add("".join(rng.choices(alphabet, k=5)))
    return sorted(probes)


def evaluate_seqle_opening_groups(
    records: Sequence[Dict[str, Any]],
    aa_props: Dict[str, Dict[str, Any]],
    profile_spec: Dict[str, Any],
    profile_id: str = "ranked",
    answer_records: Optional[Sequence[Dict[str, Any]]] = None,
    probe_seqs: Optional[Sequence[str]] = None,
    tile_only: bool = False,
    precomputed_sigs: Optional[Dict[str, Dict[str, str]]] = None,
) -> List[Dict[str, Any]]:
    """Measures expected survivors, singleton coverage & two-step probability inside every opening group."""
    by_seq = {r["seq"]: dict(r) for r in records}
    answers = (
        list(answer_records)
        if answer_records is not None
        else [r for r in records if profile_id in r.get("game_eligibility", ["learn"])]
    )
    probes = list(dict.fromkeys(probe_seqs if probe_seqs is not None else [r["seq"] for r in records]))

    probe_bio: Dict[str, Dict[str, Any]] = {}
    if precomputed_sigs is None or tile_only:
        for g_seq in probes:
            if g_seq in by_seq:
                probe_bio[g_seq] = by_seq[g_seq]
            else:
                probe_bio[g_seq] = compute_peptide_biochem(g_seq, aa_props)
        for t in answers:
            probe_bio[t["seq"]] = t

    use_class = bool(profile_spec.get("show_positional_class_match", True))

    def feedback_sig(g_seq: str, t_seq: str) -> str:
        if tile_only:
            tiles = evaluate_seqle_guess(g_seq, t_seq, aa_props)
            if not use_class:
                tiles = ["absent" if t == "class" else t for t in tiles]
            return "".join(t[0].upper() for t in tiles)
        if precomputed_sigs is not None:
            row = precomputed_sigs.get(g_seq)
            if row is not None and t_seq in row:
                return row[t_seq]
        return evaluate_seqle_feedback_from_spec(
            probe_bio[g_seq], probe_bio[t_seq], g_seq, t_seq, aa_props, profile_spec
        )["signature"]

    groups: Dict[str, List[Dict[str, Any]]] = {}
    for t in answers:
        lvl0 = _ladder_clue_text(t, profile_id, 0) or "<UNGROUPED>"
        groups.setdefault(lvl0, []).append(t)

    reports: List[Dict[str, Any]] = []
    for label in sorted(groups.keys()):
        grp = groups[label]
        n_g = len(grp)
        grp_seqs = [t["seq"] for t in grp]
        grp_seq_set = set(grp_seqs)
        best_probe = grp_seqs[0]
        best_key: Optional[Tuple[float, int, str]] = None
        best_buckets: Counter[str] = Counter()

        for g_seq in probes:
            buckets = Counter(feedback_sig(g_seq, t_seq) for t_seq in grp_seqs)
            sum_sq = sum(c * c for c in buckets.values())
            # Exact audit ranking key: minimize (sum(n*n), -len(buckets), g_seq)
            key = (sum_sq, -len(buckets), g_seq)
            if best_key is None or key < best_key:
                best_key = key
                best_probe = g_seq
                best_buckets = buckets

        exp_surv = round(best_key[0] / n_g, 4) if n_g else 0.0
        singletons = sum(1 for c in best_buckets.values() if c == 1)
        singleton_pct = round(100.0 * singletons / n_g, 2) if n_g else 0.0

        # Post-guess-1 bucket including any clue unlocking at after_guesses <= 1 (e.g. Learn level 1)
        post1_buckets: Counter[Tuple[Any, ...]] = Counter()
        best_in_grp = best_probe in grp_seq_set
        for t in grp:
            t_seq = t["seq"]
            if t_seq == best_probe:
                continue
            ladder = sorted((t.get("clue_ladder") or {}).get(profile_id) or [], key=lambda e: int(e.get("level", 0)))
            post_clues = tuple(
                c.get("text") for c in ladder if 0 < int(c.get("after_guesses", 99)) <= 1
            )
            post1_buckets[(feedback_sig(best_probe, t_seq), post_clues)] += 1
        two_step_wins = (1 if best_in_grp else 0) + len(post1_buckets)
        two_step_prob = round(two_step_wins / n_g, 4) if n_g else 0.0

        reports.append(
            {
                "group_label": label,
                "group_size": n_g,
                "probe_count": len(probes),
                "best_probe": best_probe,
                "best_expected_survivors": exp_surv,
                "best_distinct_buckets": len(best_buckets),
                "singleton_count": singletons,
                "singleton_coverage_pct": singleton_pct,
                "two_step_success_prob": two_step_prob,
                "bucket_sizes": sorted(best_buckets.values(), reverse=True),
            }
        )
    return reports


def run_seqle_greedy_solver(
    records: Sequence[Dict[str, Any]],
    aa_props: Dict[str, Dict[str, Any]],
    profile_spec: Dict[str, Any],
    profile_id: str = "ranked",
    condition_on_clues: bool = True,
    answer_records: Optional[Sequence[Dict[str, Any]]] = None,
    probe_seqs: Optional[Sequence[str]] = None,
    strategy: str = "greedy_entropy",
    strategy_label: Optional[str] = None,
    max_steps: Optional[int] = None,
    include_unconditioned: bool = True,
) -> Dict[str, Any]:
    """Runs multi-step Seqle solver separating answers & probes with full clue-ladder parity."""
    by_seq = {r["seq"]: r for r in records}
    ans_list = (
        list(answer_records)
        if answer_records is not None
        else [r for r in records if profile_id in r.get("game_eligibility", ["learn"])]
    )
    target_seqs = [r["seq"] for r in ans_list]
    for r in ans_list:
        by_seq[r["seq"]] = r

    probes = list(dict.fromkeys(probe_seqs if probe_seqs is not None else [r["seq"] for r in records]))
    is_sampled = len(probes) > len(records)
    search_method = strategy_label or (
        f"{strategy}_sampled_upper_bound" if is_sampled else f"{strategy}_exact_catalog"
    )
    mode_cap = int(profile_spec.get("max_guesses", 6 if profile_id == "ranked" else 8))
    step_cap = int(max_steps if max_steps is not None else mode_cap)
    n_t = len(target_seqs)

    # Precompute biochem for any non-catalog probe
    probe_bio: Dict[str, Dict[str, Any]] = {}
    for g_seq in probes:
        probe_bio[g_seq] = by_seq[g_seq] if g_seq in by_seq else compute_peptide_biochem(g_seq, aa_props)
    for t_seq in target_seqs:
        probe_bio[t_seq] = by_seq[t_seq]

    sig_cache: Dict[str, Dict[str, str]] = {}

    def get_sig(g_seq: str, t_seq: str) -> str:
        row = sig_cache.get(g_seq)
        if row is None:
            g_b = probe_bio[g_seq]
            row = {
                ts: evaluate_seqle_feedback_from_spec(
                    g_b, probe_bio[ts], g_seq, ts, aa_props, profile_spec
                )["signature"]
                for ts in target_seqs
            }
            sig_cache[g_seq] = row
        return row[t_seq]

    # Precompute target signatures against catalog & opening diagnostics (whole-pool diagnostic semantics)
    min_exp_rem = float("inf")
    max_distinct_targets = 0
    yigsr_distinct = 0
    for g_seq in [r["seq"] for r in records]:
        if g_seq not in probe_bio:
            probe_bio[g_seq] = by_seq[g_seq]
        cnts = Counter(get_sig(g_seq, t_seq) for t_seq in target_seqs)
        n_dist = len(cnts)
        if n_dist > max_distinct_targets:
            max_distinct_targets = n_dist
        if g_seq == "YIGSR":
            yigsr_distinct = n_dist
        exp_rem = round(sum((c * c) / n_t for c in cnts.values()), 2) if n_t else 0.0
        if exp_rem < min_exp_rem:
            min_exp_rem = exp_rem

    def clue_key_at_step(t_seq: str, step_num: int) -> Tuple[Tuple[int, Optional[str]], ...]:
        ladder = sorted((by_seq[t_seq].get("clue_ladder") or {}).get(profile_id) or [], key=lambda e: int(e.get("level", 0)))
        return tuple(
            (int(c["level"]), c.get("text"))
            for c in ladder
            if int(c.get("after_guesses", 99)) <= step_num
        )

    probe_choice_cache: Dict[Tuple[Tuple[str, ...], int, bool, str], str] = {}

    def solve_pool(cond: bool, strat: str, solve_cap: int) -> Dict[str, Tuple[int, bool]]:
        def pick_probe(cands: List[str], current_steps: int) -> str:
            if len(cands) == 1:
                return cands[0]
            cache_k = (tuple(cands), current_steps + 1 if strat == "two_step" and cond else 0, cond, strat)
            cached = probe_choice_cache.get(cache_k)
            if cached is not None:
                return cached

            n_c = len(cands)
            cands_set = set(cands)
            best_g = cands[0]
            best_key: Optional[Tuple[Any, ...]] = None

            if strat == "two_step":
                post_clues = (
                    {ts: clue_key_at_step(ts, current_steps + 1) for ts in cands}
                    if cond
                    else {ts: () for ts in cands}
                )
                for g_seq in probes:
                    sig_row = sig_cache.get(g_seq)
                    if sig_row is None:
                        g_b = probe_bio[g_seq]
                        sig_row = {
                            ts: evaluate_seqle_feedback_from_spec(
                                g_b, probe_bio[ts], g_seq, ts, aa_props, profile_spec
                            )["signature"]
                            for ts in target_seqs
                        }
                        sig_cache[g_seq] = sig_row
                    parts: Counter[str] = Counter()
                    two_buckets: set[Tuple[str, Any]] = set()
                    in_cands = g_seq in cands_set
                    for t_seq in cands:
                        s = sig_row[t_seq]
                        parts[s] += 1
                        if t_seq != g_seq:
                            two_buckets.add((s, post_clues[t_seq]))
                    two_step_wins = (1 if in_cands else 0) + len(two_buckets)
                    sum_sq = sum(c * c for c in parts.values())
                    h = -sum((c / n_c) * math.log2(c / n_c) for c in parts.values())
                    key = (two_step_wins, -sum_sq, round(h, 9), in_cands, g_seq)
                    if best_key is None or key > best_key:
                        best_g, best_key = g_seq, key
            else:
                for g_seq in probes:
                    sig_row = sig_cache.get(g_seq)
                    if sig_row is None:
                        g_b = probe_bio[g_seq]
                        sig_row = {
                            ts: evaluate_seqle_feedback_from_spec(
                                g_b, probe_bio[ts], g_seq, ts, aa_props, profile_spec
                            )["signature"]
                            for ts in target_seqs
                        }
                        sig_cache[g_seq] = sig_row
                    parts = Counter(sig_row[t_seq] for t_seq in cands)
                    h = -sum((c / n_c) * math.log2(c / n_c) for c in parts.values())
                    key = (round(h, 9), g_seq in cands_set, g_seq)
                    if best_key is None or key > best_key:
                        best_g, best_key = g_seq, key

            probe_choice_cache[cache_k] = best_g
            return best_g

        results: Dict[str, Tuple[int, bool]] = {}
        for t_seq in target_seqs:
            cands = list(target_seqs)
            if cond:
                t_key0 = clue_key_at_step(t_seq, 0)
                cands = [s for s in cands if clue_key_at_step(s, 0) == t_key0]

            steps = 0
            solved = False
            while steps < solve_cap:
                g_seq = pick_probe(cands, steps)
                steps += 1
                if g_seq == t_seq:
                    solved = True
                    break
                obs_sig = get_sig(g_seq, t_seq)
                cands = [c for c in cands if get_sig(g_seq, c) == obs_sig]
                if cond:
                    t_ck = clue_key_at_step(t_seq, steps)
                    cands = [c for c in cands if clue_key_at_step(c, steps) == t_ck]
            results[t_seq] = (steps, solved)
        return results

    guess_results = solve_pool(condition_on_clues, strategy, step_cap)
    summary = _summarize_guess_counts(guess_results, mode_cap=mode_cap)
    opening_groups = evaluate_seqle_opening_groups(
        records,
        aa_props,
        profile_spec,
        profile_id=profile_id,
        answer_records=ans_list,
        probe_seqs=probes,
        tile_only=False,
        precomputed_sigs=sig_cache,
    )
    out: Dict[str, Any] = {
        **summary,
        "eligible_count": n_t,
        "probe_count": len(probes),
        "strategy": strategy,
        "search_method": search_method,
        "yigsr_distinct_signatures": yigsr_distinct,
        "max_single_probe_distinct_targets": max_distinct_targets,
        "min_opening_expected_remaining": round(min_exp_rem, 2),
        "whole_pool_opening_metrics": {
            "yigsr_distinct_signatures": yigsr_distinct,
            "max_single_probe_distinct_targets": max_distinct_targets,
            "min_opening_expected_remaining": round(min_exp_rem, 2),
        },
        "opening_group_metrics": opening_groups,
        "per_target_guesses": {k: (v[0] if v[1] else None) for k, v in guess_results.items()},
        "per_target_results": guess_results,
        "conditioned_on_clues": condition_on_clues,
    }
    if condition_on_clues and include_unconditioned:
        raw_results = solve_pool(False, strategy, 12)
        out["unconditioned_summary"] = _summarize_guess_counts(raw_results, mode_cap=mode_cap)
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


PROVENANCE_VERIFIED_DATE = "2026-11-14"
UNIPROT_VERIFIED_RELEASE = "2026_04"
VALID_PROVENANCE_STATUSES = ("verified_uniprot", "consensus_exemplar", "abstraction")
VALID_CITATION_ROLES = ("native_sequence", "exact_structure", "motif_family_context", "motif_definition")


def validate_copi_terminal_motif(seq: str, terminal: str) -> None:
    """Validates canonical COPI dilysine C-terminal retrieval motifs (-K-K-X-X or -K-X-K-X-X)."""
    s = (seq or "").strip().upper()
    if len(s) != 5 or terminal != "C":
        raise ValueError(f"[{s}] Canonical COPI dilysine motif requires a 5-mer at terminal='C'")
    is_kkxx = s[1] == "K" and s[2] == "K"  # positions -4/-3 from C-terminus
    is_kxkxx = s[0] == "K" and s[2] == "K"  # positions -5/-3 from C-terminus
    if not (is_kkxx or is_kxkxx):
        raise ValueError(
            f"[{s}] Invalid COPI retrieval placement: canonical motifs require Lys at (-4,-3) or (-5,-3), not -KKxxx (-5,-4)."
        )


def validate_nglycosylation_sequon(seq: str, consensus: Optional[str] = None) -> None:
    """Validates eukaryotic OST N-glycosylation sequon N-[^P]-[ST] without restricting downstream residues."""
    s = (seq or "").strip().upper()
    if len(s) < 3 or s[0] != "N" or s[1] == "P" or s[2] not in ("S", "T"):
        raise ValueError(
            f"[{s}] Invalid eukaryotic N-glycosylation sequon: requires N at +0, non-P at +1, and S/T at +2."
        )
    if consensus is not None:
        tokens = _consensus_tokens(consensus)
        if len(tokens) < 3 or tokens[0] != "N" or tokens[1] != "[^P]" or set(tokens[2].strip("[]")) != {"S", "T"}:
            raise ValueError(
                f"[{s}] Glycosylation consensus {consensus!r} must state N-[^P]-[ST] with explicit non-Pro exclusion."
            )


# Authoritative reference provenance registry verified against UniProtKB/Swiss-Prot
# FASTA (release 2026_04) and RCSB PDB polymer entities (BioWordle v0.2.2).
SEQLE_REFERENCE_REGISTRY: Dict[str, Dict[str, Any]] = {
    "GHSMG": {
        "status": "verified_uniprot",
        "accession": "G9BY57",
        "isoform": "canonical",
        "parent_residue_range": [160, 170],
        "motif_residue_range": [163, 167],
        "window_offset": 3,
        "verified_parent_sequence": "AVAGHSMGGGG",
        "functional_rule": "G-X-S-X-G",
        "pdb_ids": ["4EB0"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "G9BY57", "isoform": "canonical", "parent_residue_range": [160, 170], "motif_residue_range": [163, 167]},
            {"role": "exact_structure", "source": "PDB", "id": "4EB0", "entity_id": "1", "chains": ["A"], "residue_range": [128, 132], "deposited_sequence": "VAGHSMGGG", "mapped_uniprot": "G9BY57"},
        ],
        "notes": "Leaf-branch compost cutinase (LCC); PDB:4EB0 entity 1 contains mapped GHSMG at entity 128-132.",
    },
    "EKDEL": {
        "status": "verified_uniprot",
        "accession": "P11021",
        "isoform": "canonical",
        "parent_residue_range": [648, 654],
        "motif_residue_range": [650, 654],
        "window_offset": 2,
        "verified_parent_sequence": "TAEKDEL",
        "functional_rule": "[KRDEH]-[DENQ]-E-L",
        "pdb_ids": ["6I6H"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P11021", "isoform": "canonical", "parent_residue_range": [648, 654], "motif_residue_range": [650, 654]},
            {"role": "exact_structure", "source": "PDB", "id": "6I6H", "entity_id": "2", "chains": ["B"], "residue_range": [2, 6], "deposited_sequence": "AEKDEL", "mapped_uniprot": None},
        ],
        "notes": "Human HSPA5/BiP extreme C-terminal ER retention signal; PDB:6I6H entity 2 is synthetic bound AEKDEL peptide.",
    },
    "DEVDG": {
        "status": "verified_uniprot",
        "accession": "P09874",
        "isoform": "canonical",
        "parent_residue_range": [210, 219],
        "motif_residue_range": [211, 215],
        "window_offset": 1,
        "verified_parent_sequence": "GDEVDGVDEV",
        "functional_rule": "D-X-X-D-[GSAC]",
        "pdb_ids": [],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P09874", "isoform": "canonical", "parent_residue_range": [210, 219], "motif_residue_range": [211, 215]},
            {"role": "motif_definition", "source": "ELM", "id": "CLV_Caspase3-7"},
        ],
        "notes": "Human PARP1 Caspase-3/7 cleavage site (cleaved after Asp214).",
    },
    "GRGDS": {
        "status": "verified_uniprot",
        "accession": "P02751",
        "isoform": "canonical",
        "parent_residue_range": [1612, 1622],
        "motif_residue_range": [1614, 1618],
        "window_offset": 2,
        "verified_parent_sequence": "VTGRGDSPASS",
        "functional_rule": "R-G-D-[ST]",
        "pdb_ids": ["1FNF"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P02751", "isoform": "canonical", "parent_residue_range": [1612, 1622], "motif_residue_range": [1614, 1618]},
            {"role": "exact_structure", "source": "PDB", "id": "1FNF", "entity_id": "1", "chains": ["A"], "residue_range": [351, 355], "deposited_sequence": "VTGRGDSPA", "mapped_uniprot": "P02751"},
        ],
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
        "functional_rule": "L-P-X-T-G",
        "pdb_ids": ["1T2W"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P02976", "isoform": "canonical", "parent_residue_range": [480, 489], "motif_residue_range": [482, 486]},
            {"role": "exact_structure", "source": "PDB", "id": "1T2W", "entity_id": "2", "chains": ["D"], "residue_range": [1, 5], "deposited_sequence": "LPETG", "mapped_uniprot": None},
        ],
        "notes": "Staphylococcus aureus Protein A Sortase-A cell-wall sorting signal; PDB:1T2W entity 2 is substrate LPETG peptide.",
    },
    "HEIGH": {
        "status": "verified_uniprot",
        "accession": "P08254",
        "isoform": "canonical",
        "parent_residue_range": [215, 226],
        "motif_residue_range": [218, 222],
        "window_offset": 3,
        "verified_parent_sequence": "VAAHEIGHSLGL",
        "functional_rule": "H-E-X-X-H",
        "pdb_ids": ["1B3D"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P08254", "isoform": "canonical", "parent_residue_range": [215, 226], "motif_residue_range": [218, 222]},
            {"role": "exact_structure", "source": "PDB", "id": "1B3D", "entity_id": "1", "chains": ["A", "B"], "residue_range": [119, 123], "deposited_sequence": "AAHEIGHSL", "mapped_uniprot": "P08254"},
        ],
        "notes": "Human MMP3 Stromelysin-1 Zn-binding helix.",
    },
    "CCHHC": {
        "status": "abstraction",
        "accession": "PROSITE:PS50158",
        "isoform": None,
        "parent_residue_range": None,
        "motif_residue_range": None,
        "window_offset": 0,
        "verified_parent_sequence": "CCHHC",
        "functional_rule": "C-X2-C-X4-H-X4-C",
        "pdb_ids": [],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": "N/A",
        "citations": [
            {"role": "motif_definition", "source": "PROSITE", "id": "PS50158"},
            {"role": "motif_family_context", "source": "PDB", "id": "1MFS", "notes": "Retroviral gag zinc knuckle NMR structure (spaced C-X2-C-X4-H-X4-C ligands; not contiguous)."},
        ],
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
        "functional_rule": "Y-I-G-S-R",
        "pdb_ids": [],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P07942", "isoform": "canonical", "parent_residue_range": [946, 956], "motif_residue_range": [950, 954]},
        ],
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
        "functional_rule": "F-X-N-P-X-Y",
        "pdb_ids": [],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "motif_definition", "source": "ELM", "id": "TRG_ENDOCYTIC_2"},
            {"role": "motif_family_context", "source": "PDB", "id": "3SO6", "entity_id": "2", "chains": ["Q"], "mapped_uniprot": "P01130", "variant_sequence": "NPVYQ", "variant_residue_range": [825, 829], "notes": "Natural LDLR window carries NPVYQ (no universal Lys/Arg after Tyr established by LDLR/3SO6 evidence)."},
        ],
        "notes": "Synthetic teaching representative of F-X-N-P-X-Y core; native human LDLR P01130 (825-829) and PDB:3SO6 carry NPVYQ.",
    },
    "KCVIM": {
        "status": "verified_uniprot",
        "accession": "P01116-2",
        "isoform": "2B",
        "parent_residue_range": [175, 188],
        "motif_residue_range": [184, 188],
        "window_offset": 9,
        "verified_parent_sequence": "KKKKKKSKTKCVIM",
        "functional_rule": "C-a-a-X (C-terminal)",
        "pdb_ids": ["1D8D"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P01116-2", "isoform": "2B", "parent_residue_range": [175, 188], "motif_residue_range": [184, 188]},
            {"role": "exact_structure", "source": "PDB", "id": "1D8D", "entity_id": "3", "chains": ["P"], "residue_range": [7, 11], "deposited_sequence": "KTKCVIM", "mapped_uniprot": "P01116"},
        ],
        "notes": "Human KRAS Isoform 2B polybasic C-terminal CaaX farnesylation box; RCSB SIFTS maps 1D8D entity 3 to base P01116.",
    },
    "DEADE": {
        "status": "verified_uniprot",
        "accession": "P60842",
        "isoform": "canonical",
        "parent_residue_range": [180, 190],
        "motif_residue_range": [182, 186],
        "window_offset": 2,
        "verified_parent_sequence": "VLDEADEMLSR",
        "functional_rule": "D-E-A-D-[DE]",
        "pdb_ids": ["2ZU6"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P60842", "isoform": "canonical", "parent_residue_range": [180, 190], "motif_residue_range": [182, 186]},
            {"role": "exact_structure", "source": "PDB", "id": "2ZU6", "entity_id": "1", "chains": ["A", "C", "D", "F"], "residue_range": [164, 168], "deposited_sequence": "VLDEADEML", "mapped_uniprot": "P60842"},
        ],
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
        "functional_rule": "G-X4-G-K-[TS]",
        "pdb_ids": ["4OBE"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P01116", "isoform": "canonical", "parent_residue_range": [10, 17], "motif_residue_range": [13, 17]},
            {"role": "exact_structure", "source": "PDB", "id": "4OBE", "entity_id": "1", "chains": ["A", "B"], "residue_range": [14, 18], "deposited_sequence": "AGGVGKSAL", "mapped_uniprot": "P01116"},
        ],
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
        "functional_rule": "[VI]-D-F-G-[LM]",
        "pdb_ids": ["4Y83"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P41279", "isoform": "canonical", "parent_residue_range": [265, 277], "motif_residue_range": [269, 273]},
            {"role": "exact_structure", "source": "PDB", "id": "4Y83", "entity_id": "1", "chains": ["A", "B", "C"], "residue_range": [206, 210], "deposited_sequence": "VLVDFGLSV", "mapped_uniprot": "P41279"},
        ],
        "notes": "Human MAP3K8/COT activation-loop DFG switch.",
    },
    "HRDLA": {
        "status": "verified_uniprot",
        "accession": "P10721",
        "isoform": "canonical",
        "parent_residue_range": [789, 797],
        "motif_residue_range": [790, 794],
        "window_offset": 1,
        "verified_parent_sequence": "IHRDLAARN",
        "functional_rule": "H-R-D-L-A",
        "pdb_ids": ["1PKG"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P10721", "isoform": "canonical", "parent_residue_range": [789, 797], "motif_residue_range": [790, 794]},
            {"role": "exact_structure", "source": "PDB", "id": "1PKG", "entity_id": "1", "chains": ["A", "B"], "residue_range": [184, 188], "deposited_sequence": "CIHRDLAAR", "mapped_uniprot": "P10721"},
        ],
        "notes": "Human KIT catalytic loop.",
    },
    "WSEWS": {
        "status": "verified_uniprot",
        "accession": "P32927",
        "isoform": "canonical",
        "parent_residue_range": [421, 433],
        "motif_residue_range": [425, 429],
        "window_offset": 4,
        "verified_parent_sequence": "YNGIWSEWSEARS",
        "functional_rule": "W-S-X-W-S",
        "pdb_ids": ["1EGJ"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P32927", "isoform": "canonical", "parent_residue_range": [421, 433], "motif_residue_range": [425, 429]},
            {"role": "exact_structure", "source": "PDB", "id": "1EGJ", "entity_id": "1", "chains": ["A"], "residue_range": [88, 92], "deposited_sequence": "GIWSEWSEA", "mapped_uniprot": "P32927"},
        ],
        "notes": "Human CSF2RB cytokine receptor WSXWS box.",
    },
    "PPPVP": {
        "status": "verified_uniprot",
        "accession": "Q07889",
        "isoform": "canonical",
        "parent_residue_range": [1149, 1158],
        "motif_residue_range": [1150, 1154],
        "window_offset": 1,
        "verified_parent_sequence": "VPPPVPPRRR",
        "functional_rule": "P-P-X-[VIL]-P",
        "pdb_ids": [],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "Q07889", "isoform": "canonical", "parent_residue_range": [1149, 1158], "motif_residue_range": [1150, 1154]},
            {"role": "motif_family_context", "source": "PDB", "id": "1GBR", "entity_id": "2", "chains": ["B"], "mapped_uniprot": "Q02384", "variant_sequence": "SPLLPKLPPKTYKRE", "notes": "PDB:1GBR entity 2 is mouse Sos1 (Q02384) peptide SPLLPKLPPKTYKRE (motif-family SH3 context only; lacks human PPPVP)."},
            {"role": "motif_definition", "source": "ELM", "id": "LIG_SH3_2"},
        ],
        "notes": "Human SOS1 Class-II SH3-binding polyproline motif (1150-1154); PDB:1GBR is murine SH3-ligand context only.",
    },
    "INFEK": {
        "status": "verified_uniprot",
        "accession": "P01012",
        "isoform": "canonical",
        "parent_residue_range": [258, 265],
        "motif_residue_range": [260, 264],
        "window_offset": 2,
        "verified_parent_sequence": "SIINFEKL",
        "functional_rule": "S-I-I-N-F-E-K-L",
        "pdb_ids": ["1VAC"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P01012", "isoform": "canonical", "parent_residue_range": [258, 265], "motif_residue_range": [260, 264]},
            {"role": "exact_structure", "source": "PDB", "id": "1VAC", "entity_id": "3", "chains": ["P"], "residue_range": [3, 7], "deposited_sequence": "SIINFEKL", "mapped_uniprot": "P01012"},
        ],
        "notes": "Gallus gallus Ovalbumin H-2Kb epitope SIINFEKL (UniProt:P01012 residues 260-264 within 258-265).",
    },
    "YQRLN": {
        "status": "verified_uniprot",
        "accession": "P19814",
        "isoform": "canonical",
        "parent_residue_range": [348, 357],
        "motif_residue_range": [350, 354],
        "window_offset": 2,
        "verified_parent_sequence": "SDYQRLNLKL",
        "functional_rule": "Y-X-X-[LIVF]-N",
        "pdb_ids": ["1BXX"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P19814", "isoform": "canonical", "parent_residue_range": [348, 357], "motif_residue_range": [350, 354]},
            {"role": "exact_structure", "source": "PDB", "id": "1BXX", "entity_id": "2", "chains": ["P"], "residue_range": [2, 6], "deposited_sequence": "DYQRLN", "mapped_uniprot": None},
        ],
        "notes": "Rat TGN38 tyrosine-based AP2 internalization motif; PDB:1BXX entity 2 is synthetic DYQRLN peptide.",
    },
    "KEKSD": {
        "status": "verified_uniprot",
        "accession": "P39656",
        "isoform": "canonical",
        "parent_residue_range": [437, 456],
        "motif_residue_range": [452, 456],
        "window_offset": 15,
        "verified_parent_sequence": "LGLFIFSIVFLHMKEKEKSD",
        "functional_rule": "K-X-K-X-X (C-terminal)",
        "functional_rule_id": "COPI_DILYSINE_C_TERM",
        "pdb_ids": ["4J77"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P39656", "isoform": "canonical", "parent_residue_range": [437, 456], "motif_residue_range": [452, 456]},
            {"role": "exact_structure", "source": "PDB", "id": "4J77", "entity_id": "1", "chains": ["C", "D"], "residue_range": [2, 6], "deposited_sequence": "AKEKSD", "mapped_uniprot": None},
            {"role": "motif_definition", "source": "ELM", "id": "TRG_ER_diLys_1"},
        ],
        "notes": "Human OST48/DDOST extreme C-terminal -KxKxx COPI retrieval signal (452-456); PDB:4J77 entity 1 is co-crystallized Wbp1 C-terminal peptide AKEKSD.",
    },
    "CAACH": {
        "status": "verified_uniprot",
        "accession": "P07143",
        "isoform": "canonical",
        "parent_residue_range": [97, 109],
        "motif_residue_range": [101, 105],
        "window_offset": 4,
        "verified_parent_sequence": "YREVCAACHSLDR",
        "functional_rule": "C-X-X-C-H",
        "pdb_ids": ["3CX5"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P07143", "isoform": "canonical", "parent_residue_range": [97, 109], "motif_residue_range": [101, 105]},
            {"role": "exact_structure", "source": "PDB", "id": "3CX5", "entity_id": "4", "chains": ["D", "O"], "residue_range": [40, 44], "deposited_sequence": "EVCAACHSL", "mapped_uniprot": "P07143"},
        ],
        "notes": "Yeast mitochondrial Cytochrome c1 CXXCH heme site.",
    },
    "NLTGY": {
        "status": "consensus_exemplar",
        "accession": "ELM:MOD_N-GLC_1",
        "isoform": None,
        "parent_residue_range": None,
        "motif_residue_range": None,
        "window_offset": 0,
        "verified_parent_sequence": "NLTGY",
        "functional_rule": "N-[^P]-[ST]",
        "functional_rule_id": "N_GLC_SEQUON",
        "pdb_ids": [],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "motif_definition", "source": "ELM", "id": "MOD_N-GLC_1"},
            {"role": "motif_family_context", "source": "PDB", "id": "3WAK", "entity_id": "1", "chains": ["A"], "mapped_uniprot": "O29867", "notes": "PDB:3WAK is apo Archaeoglobus fulgidus AglB oligosaccharyltransferase (no bound acceptor peptide; UniProt O29867 lacks NLTGY)."},
        ],
        "notes": "Synthetic teaching example for the N-[^P]-[ST] eukaryotic N-glycosylation sequon (G-Y padding is pedagogical).",
    },
    "DSGLD": {
        "status": "verified_uniprot",
        "accession": "P25963",
        "isoform": "canonical",
        "parent_residue_range": [27, 39],
        "motif_residue_range": [31, 35],
        "window_offset": 4,
        "verified_parent_sequence": "DDRHDSGLDSMKD",
        "functional_rule": "D-S-G-X-X-S",
        "pdb_ids": ["9T9W"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P25963", "isoform": "canonical", "parent_residue_range": [27, 39], "motif_residue_range": [31, 35]},
            {"role": "exact_structure", "source": "PDB", "id": "9T9W", "entity_id": "2", "chains": ["C"], "residue_range": [4, 8], "deposited_sequence": "RHDSGLDSM", "mapped_uniprot": "P25963"},
            {"role": "motif_definition", "source": "ELM", "id": "DEG_SCF_TRCP1_1"},
        ],
        "notes": "Human NFKBIA / I-kappa-B-alpha beta-TrCP phosphodegron (PDB:9T9W entity 2 directly deposits phosphorylated I-kappa-B-alpha DRHDSGLDSMKD).",
    },
    "KENIP": {
        "status": "consensus_exemplar",
        "accession": "ELM:DEG_APCC_KENBOX_2",
        "isoform": None,
        "parent_residue_range": None,
        "motif_residue_range": None,
        "window_offset": 0,
        "verified_parent_sequence": "KENIP",
        "functional_rule": "K-E-N",
        "pdb_ids": [],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "motif_definition", "source": "ELM", "id": "DEG_APCC_KENBOX_2"},
            {"role": "motif_family_context", "source": "PDB", "id": "4BH6", "entity_id": "2", "chains": ["I", "J", "K", "L", "M", "N", "O", "P"], "mapped_uniprot": "Q08981", "variant_sequence": "KENLS", "variant_residue_range": [98, 102], "notes": "PDB:4BH6 Acm1 (Q08981) carries pseudosubstrate KENLS; Hsl1 (P34244 774-778) carries NKENE (neither requires downstream Pro)."},
        ],
        "notes": "Synthetic teaching example of the APC/C KEN-box core (K-E-N; I-P downstream padding is pedagogical, not a universal Pro requirement).",
    },
    "KKKRK": {
        "status": "verified_uniprot",
        "accession": "P03070",
        "isoform": "canonical",
        "parent_residue_range": [126, 132],
        "motif_residue_range": [127, 131],
        "window_offset": 1,
        "verified_parent_sequence": "PKKKRKV",
        "functional_rule": "[PKR]-K-[KR]-X-[KR]",
        "pdb_ids": ["1BK6"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P03070", "isoform": "canonical", "parent_residue_range": [126, 132], "motif_residue_range": [127, 131]},
            {"role": "exact_structure", "source": "PDB", "id": "1BK6", "entity_id": "2", "chains": ["C", "E"], "residue_range": [1, 5], "deposited_sequence": "KKKRKV", "mapped_uniprot": None},
        ],
        "notes": "SV40 Large T-antigen monopartite nuclear localization signal; PDB:1BK6 entity 2 is NLS peptide KKKRKV.",
    },
    "FGGFT": {
        "status": "verified_uniprot",
        "accession": "Q13519",
        "isoform": "canonical",
        "parent_residue_range": [130, 146],
        "motif_residue_range": [130, 134],
        "window_offset": 0,
        "verified_parent_sequence": "FGGFTGARKSARKLANQ",
        "functional_rule": "F-G-G-F-T (N-terminal message)",
        "pdb_ids": [],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "Q13519", "isoform": "canonical", "parent_residue_range": [130, 146], "motif_residue_range": [130, 134]},
        ],
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
        "functional_rule": "Y-G-G-F-M (N-terminal opioid message)",
        "pdb_ids": ["1PLW"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P01210", "isoform": "canonical", "parent_residue_range": [100, 104], "motif_residue_range": [100, 104]},
            {"role": "exact_structure", "source": "PDB", "id": "1PLW", "entity_id": "1", "chains": ["A"], "residue_range": [1, 5], "deposited_sequence": "YGGFM", "mapped_uniprot": "P01210"},
        ],
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
        "functional_rule": "[SG]-Y-[SG]-G-Y",
        "pdb_ids": ["5W3N"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P35637", "isoform": "canonical", "parent_residue_range": [34, 43], "motif_residue_range": [37, 41]},
            {"role": "exact_structure", "source": "PDB", "id": "5W3N", "entity_id": "1", "chains": ["A", "B", "C", "D", "E", "F", "G", "H", "I"], "residue_range": [64, 68], "deposited_sequence": "QQSYSGYSQ", "mapped_uniprot": "P35637"},
        ],
        "notes": "Human FUS low-complexity domain canonical repeat 37-41.",
    },
    "LDLAS": {
        "status": "abstraction",
        "accession": "ELM:TRG_NES_CRM1_1",
        "isoform": None,
        "parent_residue_range": None,
        "motif_residue_range": None,
        "window_offset": 0,
        "verified_parent_sequence": "LDLAS",
        "functional_rule": "L-X(2,3)-L-X(2,3)-L-X-L",
        "pdb_ids": [],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": "N/A",
        "citations": [
            {"role": "motif_definition", "source": "ELM", "id": "TRG_NES_CRM1_1"},
            {"role": "motif_family_context", "source": "PDB", "id": "3NBY", "notes": "CRM1-RanGTP nuclear export complex context (native NES spans 10-12 aa, not contiguous 5-mer)."},
        ],
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
        "functional_rule": "N-P-L-Y-K",
        "pdb_ids": ["1MK7"],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P05106", "isoform": "canonical", "parent_residue_range": [763, 788], "motif_residue_range": [770, 774]},
            {"role": "exact_structure", "source": "PDB", "id": "1MK7", "entity_id": "1", "chains": ["A", "C"], "residue_range": [10, 14], "deposited_sequence": "ANNPLYKE", "mapped_uniprot": "P05106"},
        ],
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
        "functional_rule": "G-H-H-P-H",
        "pdb_ids": [],
        "provenance_verified_date": PROVENANCE_VERIFIED_DATE,
        "uniprot_release": UNIPROT_VERIFIED_RELEASE,
        "citations": [
            {"role": "native_sequence", "source": "UniProt", "id": "P04196", "isoform": "canonical", "parent_residue_range": [383, 392], "motif_residue_range": [383, 387]},
        ],
        "notes": "Human HRG tandem HPR repeat 383-392 (motif 383-387; single-copy repeat also at 363-367).",
    },
}


def validate_seqle_registry_completeness(
    records: Sequence[Dict[str, Any]],
    registry: Optional[Dict[str, Dict[str, Any]]] = None,
) -> None:
    """Verifies a strict 1:1 bijective mapping between shipped Seqle records and the reference registry."""
    reg = SEQLE_REFERENCE_REGISTRY if registry is None else registry
    catalog_seqs = {str(r.get("seq", "")).strip().upper() for r in records}
    reg_seqs = set(reg.keys())
    missing = sorted(catalog_seqs - reg_seqs)
    stale = sorted(reg_seqs - catalog_seqs)
    if missing:
        raise ValueError(f"Shipped Seqle catalog has motifs missing from reference registry: {missing}")
    if stale:
        raise ValueError(f"Reference registry contains stale extra entries not in shipped catalog: {stale}")


def validate_seqle_motif_record(
    item: Dict[str, Any],
    aa_props: Dict[str, Dict[str, Any]],
    enforce_reference_registry: bool = True,
    registry: Optional[Dict[str, Dict[str, Any]]] = None,
) -> None:
    """Validates a Seqle motif record against P0-4 schema, biological rules, and typed provenance."""
    reg = SEQLE_REFERENCE_REGISTRY if registry is None else registry
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
        if len(tokens) == 5:
            regex = "^" + "".join(tokens) + "$"
            if not re.match(regex, seq):
                raise ValueError(f"[{seq}] Sequence does not match its 5-position fragment consensus {consensus!r}")

    elig = item.get("game_eligibility", [])
    if not elig or not set(elig).issubset({"learn", "ranked"}):
        raise ValueError(f"[{seq}] Invalid game_eligibility {elig!r}")
    if kind == "abstraction" and "ranked" in elig:
        raise ValueError(f"[{seq}] kind='abstraction' cannot be eligible for 'ranked' mode")

    if enforce_reference_registry:
        if seq not in reg:
            raise ValueError(f"[{seq}] Motif missing from SEQLE_REFERENCE_REGISTRY (enforce_reference_registry=True)")
        reg_ref = reg[seq]
        ref = item.get("reference_provenance") if isinstance(item.get("reference_provenance"), dict) else reg_ref
        status = item.get("provenance_status") if "provenance_status" in item else ref.get("status")
        if status not in VALID_PROVENANCE_STATUSES or status != reg_ref.get("status"):
            raise ValueError(f"[{seq}] Unknown or mismatched provenance status {status!r}")
        if parent_seq != ref.get("verified_parent_sequence"):
            raise ValueError(
                f"[{seq}] parent_sequence {parent_seq!r} does not match verified reference "
                f"{ref.get('verified_parent_sequence')!r} ({ref.get('accession')})"
            )
        if offset != ref.get("window_offset"):
            raise ValueError(
                f"[{seq}] window_offset {offset} != verified offset {ref.get('window_offset')} ({ref.get('accession')})"
            )

        citations = item.get("citations") if "citations" in item else ref.get("citations")
        if not isinstance(citations, list) or not citations:
            raise ValueError(f"[{seq}] Provenance requires non-empty typed citations list")
        exact_pdb_ids: set[str] = set()
        context_pdb_ids: set[str] = set()
        for cit in citations:
            if not isinstance(cit, dict):
                raise ValueError(f"[{seq}] Malformed citation entry: {cit!r}")
            role = cit.get("role")
            if role not in VALID_CITATION_ROLES:
                raise ValueError(f"[{seq}] Invalid citation role {role!r}")
            if not cit.get("source") or not cit.get("id"):
                raise ValueError(f"[{seq}] Citation missing source/id: {cit!r}")
            if role == "native_sequence" and cit.get("source") == "UniProt":
                if cit.get("id") != ref.get("accession"):
                    raise ValueError(f"[{seq}] native_sequence UniProt id {cit.get('id')!r} != {ref.get('accession')!r}")
                reg_native = next(
                    (rc for rc in (reg_ref.get("citations") or []) if isinstance(rc, dict) and rc.get("role") == "native_sequence"),
                    None,
                )
                expected_iso = (reg_native.get("isoform") if reg_native else None) or "canonical"
                if cit.get("isoform") != expected_iso:
                    raise ValueError(f"[{seq}] native_sequence isoform {cit.get('isoform')!r} != {expected_iso!r}")
                if cit.get("motif_residue_range") != ref.get("motif_residue_range"):
                    raise ValueError(
                        f"[{seq}] native_sequence motif_residue_range {cit.get('motif_residue_range')!r} != {ref.get('motif_residue_range')!r}"
                    )
                if cit.get("parent_residue_range") != ref.get("parent_residue_range"):
                    raise ValueError(
                        f"[{seq}] native_sequence parent_residue_range {cit.get('parent_residue_range')!r} != {ref.get('parent_residue_range')!r}"
                    )
            if role == "exact_structure":
                if status in ("consensus_exemplar", "abstraction"):
                    raise ValueError(f"[{seq}] {status} cannot claim exact_structure citation {cit['id']!r}")
                ent = cit.get("entity_id")
                chains = cit.get("chains")
                r_rng = cit.get("residue_range")
                dep_seq = str(cit.get("deposited_sequence") or "")
                if not ent or not isinstance(chains, list) or not chains:
                    raise ValueError(f"[{seq}] exact_structure {cit['id']!r} must store entity_id and chains")
                if not (isinstance(r_rng, list) and len(r_rng) == 2 and all(isinstance(v, int) for v in r_rng) and r_rng[0] >= 1 and r_rng[1] - r_rng[0] + 1 == 5):
                    raise ValueError(f"[{seq}] exact_structure {cit['id']!r} invalid residue_range {r_rng!r}")
                if seq not in dep_seq:
                    raise ValueError(f"[{seq}] exact_structure {cit['id']!r} deposited_sequence {dep_seq!r} lacks motif")
                mapped = cit.get("mapped_uniprot")
                if status == "verified_uniprot" and mapped is not None:
                    base_acc = str(ref["accession"]).split("-")[0]
                    if mapped != ref["accession"] and mapped != base_acc:
                        raise ValueError(
                            f"[{seq}] Conflicting mapped_uniprot {mapped!r} in exact_structure {cit['id']!r} vs accession {ref['accession']!r}"
                        )
                exact_pdb_ids.add(str(cit["id"]))
            elif role == "motif_family_context" and cit.get("source") == "PDB":
                context_pdb_ids.add(str(cit["id"]))

        if set(ref.get("pdb_ids") or []) != exact_pdb_ids:
            raise ValueError(
                f"[{seq}] Registry pdb_ids {ref.get('pdb_ids')!r} disagree with exact_structure citations {sorted(exact_pdb_ids)!r}"
            )
        if exact_pdb_ids & context_pdb_ids:
            raise ValueError(f"[{seq}] PDB ID(s) {sorted(exact_pdb_ids & context_pdb_ids)} marked both exact and context")

        if status == "verified_uniprot":
            acc = str(ref["accession"])
            if f"UniProt:{acc}" not in exemplar:
                raise ValueError(f"[{seq}] exemplar_source {exemplar!r} must cite verified UniProt:{acc}")
            p_rng = ref["parent_residue_range"]
            m_rng = ref["motif_residue_range"]
            if not (isinstance(p_rng, list) and len(p_rng) == 2 and p_rng[1] - p_rng[0] + 1 == len(parent_seq)):
                raise ValueError(f"[{seq}] Invalid parent_residue_range {p_rng!r} for len={len(parent_seq)}")
            expected_motif_start = p_rng[0] + offset
            if m_rng != [expected_motif_start, expected_motif_start + 4]:
                raise ValueError(
                    f"[{seq}] motif_residue_range {m_rng!r} inconsistent with parent range {p_rng!r} + offset {offset}"
                )
            cited_coords = f"{m_rng[0]}-{m_rng[1]}"
            # Reject stale wrong range mentions when explicit coordinates appear in exemplar_source
            for m_coord in re.findall(r"\b(\d+-\d+)\b", exemplar):
                if m_coord != cited_coords and m_coord != f"{p_rng[0]}-{p_rng[1]}":
                    # allow documented legacy mature coordinate parentheticals ONLY if canonical range is present
                    if cited_coords not in exemplar:
                        raise ValueError(f"[{seq}] exemplar_source {exemplar!r} cites {m_coord} instead of canonical {cited_coords}")
        elif status in ("consensus_exemplar", "abstraction"):
            if parent_seq != seq or offset != 0:
                raise ValueError(
                    f"[{seq}] {status} must not fabricate flanking residues (got parent={parent_seq!r}, offset={offset})"
                )
            if re.search(r"UniProt:[A-Z0-9]+", exemplar) and "context" not in exemplar.lower() and "carries" not in exemplar.lower() and "consensus" not in exemplar.lower():
                raise ValueError(f"[{seq}] {status} exemplar_source {exemplar!r} cannot imply direct native UniProt ownership")
            if status == "abstraction" and kind != "abstraction":
                raise ValueError(f"[{seq}] abstraction status requires kind='abstraction'")

        rule_id = ref.get("functional_rule_id")
        if rule_id == "COPI_DILYSINE_C_TERM":
            validate_copi_terminal_motif(seq, str(terminal))
        elif rule_id == "N_GLC_SEQUON":
            validate_nglycosylation_sequon(seq, consensus)


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


# ---------------------------------------------------------------------------
# 10. Offline Chemle Compound Manifest Validator (Commit 3)
# ---------------------------------------------------------------------------
ALLOWED_CHEMLE_REVIEW_STATUSES = frozenset(
    {"verified", "baseline_legacy", "pending", "draft", "rejected"}
)
PLAYABLE_CHEMLE_REVIEW_STATUSES = frozenset({"verified", "baseline_legacy"})
ALLOWED_CHEMLE_STEREO_STATUSES = frozenset(
    {"achiral", "defined_stereo", "racemate", "unspecified_legacy"}
)
ALLOWED_CHEMLE_PROTONATION_POLICIES = frozenset(
    {"neutral_canonical", "neutral_major_tautomer", "zwitterion_neutral_net"}
)
_CHEMLE_VALID_PROFILES = frozenset({"learn", "ranked"})
_CHEMLE_UNCHARGER = rdMolStandardize.Uncharger()
_CHEMLE_LARGEST_FRAG = rdMolStandardize.LargestFragmentChooser()
_NITRO_GROUP_SMARTS = Chem.MolFromSmarts("[N+](=O)[O-]")


def validate_chemle_manifest_entry(item: Dict[str, Any]) -> Tuple[Chem.Mol, Dict[str, Any], str]:
    """Validates one Chemle manifest compound strictly offline via RDKit.

    Returns `(std_mol, std_info, parent_uncharged_inchikey)`.
    Raises `ValueError` on schema, stereochemical, salt/parent, formula,
    label/graph, or unreviewed playable eligibility violations.
    """
    if not isinstance(item, dict):
        raise ValueError("Manifest entry must be a JSON object")

    mid = str(item.get("id") or "").strip()
    if not mid or not re.match(r"^[a-z0-9_]+$", mid):
        raise ValueError(f"Invalid or missing manifest id: {item.get('id')!r}")

    required_str_fields = (
        "id",
        "name",
        "Category",
        "primary_class",
        "input_isomeric_smiles",
        "expected_standardized_smiles",
        "expected_inchikey",
        "expected_formula",
        "stereo_status",
        "protonation_policy",
        "review_status",
        "hint_role",
        "hint_origin",
        "structural_note",
    )
    for field in required_str_fields:
        if not isinstance(item.get(field), str) or not item[field].strip():
            raise ValueError(f"[{mid}] Missing or empty required string field {field!r}")

    review_status = item["review_status"].strip()
    if review_status not in ALLOWED_CHEMLE_REVIEW_STATUSES:
        raise ValueError(f"[{mid}] Unsupported review_status: {review_status!r}")

    stereo_status = item["stereo_status"].strip()
    if stereo_status not in ALLOWED_CHEMLE_STEREO_STATUSES:
        raise ValueError(f"[{mid}] Unsupported stereo_status: {stereo_status!r}")

    protonation_policy = item["protonation_policy"].strip()
    if protonation_policy not in ALLOWED_CHEMLE_PROTONATION_POLICIES:
        raise ValueError(f"[{mid}] Unsupported protonation_policy: {protonation_policy!r}")

    role_tags = item.get("role_tags")
    if not isinstance(role_tags, list) or not role_tags or not all(
        isinstance(t, str) and t.strip() for t in role_tags
    ):
        raise ValueError(f"[{mid}] role_tags must be a non-empty list of strings")

    aliases = item.get("aliases")
    if aliases is None or not isinstance(aliases, list) or not all(
        isinstance(a, str) and a.strip() for a in aliases
    ):
        raise ValueError(f"[{mid}] aliases must be a list of non-empty strings")

    g_elig = item.get("game_eligibility")
    q_elig = item.get("guess_eligibility")
    if not isinstance(g_elig, list) or not isinstance(q_elig, list):
        raise ValueError(f"[{mid}] game_eligibility and guess_eligibility must be lists")
    g_set = set(g_elig)
    q_set = set(q_elig)
    if not g_set.issubset(_CHEMLE_VALID_PROFILES) or not q_set.issubset(_CHEMLE_VALID_PROFILES):
        raise ValueError(f"[{mid}] Unknown profile in eligibility: game={g_elig!r}, guess={q_elig!r}")
    if not g_set.issubset(q_set):
        raise ValueError(f"[{mid}] game_eligibility {g_elig!r} must be a subset of guess_eligibility {q_elig!r}")

    if (g_set or q_set) and review_status not in PLAYABLE_CHEMLE_REVIEW_STATUSES:
        raise ValueError(
            f"[{mid}] Unreviewed record (review_status={review_status!r}) cannot have playable eligibility"
        )

    s_refs = item.get("structure_refs")
    if not isinstance(s_refs, dict):
        raise ValueError(f"[{mid}] structure_refs must be an object")
    cid = s_refs.get("pubchem_cid")
    if isinstance(cid, bool) or not isinstance(cid, int) or cid <= 0:
        raise ValueError(f"[{mid}] structure_refs.pubchem_cid must be a positive integer, got {cid!r}")
    chebi_id = s_refs.get("chebi_id")
    if not isinstance(chebi_id, str) or not re.match(r"^CHEBI:\d+$", chebi_id):
        raise ValueError(f"[{mid}] structure_refs.chebi_id must match 'CHEBI:<digits>', got {chebi_id!r}")

    e_refs = item.get("evidence_refs")
    if not isinstance(e_refs, list) or len(e_refs) == 0:
        raise ValueError(f"[{mid}] evidence_refs must be a non-empty list")
    for idx, ref in enumerate(e_refs):
        if not isinstance(ref, dict) or not all(
            isinstance(ref.get(k), str) and ref[k].strip() for k in ("source", "id", "url", "role")
        ):
            raise ValueError(f"[{mid}] evidence_refs[{idx}] missing source/id/url/role")

    raw_smi = item["input_isomeric_smiles"].strip()
    if "." in raw_smi:
        raise ValueError(
            f"[{mid}] Disconnected multi-fragment / salt input SMILES disallowed in manifest: {raw_smi!r}"
        )

    raw_mol = Chem.MolFromSmiles(raw_smi)
    if raw_mol is None:
        raise ValueError(f"[{mid}] Invalid input_isomeric_smiles: {raw_smi!r}")

    std_mol, std_info = standardize_and_audit(raw_smi)
    if std_info["canonical_smiles"] != item["expected_standardized_smiles"]:
        raise ValueError(
            f"[{mid}] expected_standardized_smiles mismatch: expected {item['expected_standardized_smiles']!r}, got {std_info['canonical_smiles']!r}"
        )
    if std_info["inchikey"] != item["expected_inchikey"]:
        raise ValueError(
            f"[{mid}] expected_inchikey mismatch: expected {item['expected_inchikey']!r}, got {std_info['inchikey']!r}"
        )

    computed_formula = rdMolDescriptors.CalcMolFormula(std_mol)
    if computed_formula != item["expected_formula"]:
        raise ValueError(
            f"[{mid}] Label/formula conflict: expected {item['expected_formula']!r}, got {computed_formula!r}"
        )

    net_q = Chem.GetFormalCharge(std_mol)
    nitro_atom_idxs = {
        idx for match in std_mol.GetSubstructMatches(_NITRO_GROUP_SMARTS) for idx in match
    }
    nonzero_atoms = [
        a.GetFormalCharge()
        for a in std_mol.GetAtoms()
        if a.GetFormalCharge() != 0 and a.GetIdx() not in nitro_atom_idxs
    ]
    if protonation_policy in ("neutral_canonical", "neutral_major_tautomer"):
        if net_q != 0 or nonzero_atoms:
            raise ValueError(
                f"[{mid}] Protonation policy {protonation_policy!r} expects uncharged neutral graph, got charge={net_q}"
            )
    elif protonation_policy == "zwitterion_neutral_net":
        if net_q != 0:
            raise ValueError(f"[{mid}] Zwitterion policy expects net charge 0, got {net_q}")

    centers = Chem.FindMolChiralCenters(std_mol, includeUnassigned=True)
    pot_stereo = list(Chem.FindPotentialStereo(std_mol))
    unspec_stereo = [p for p in pot_stereo if p.specified == Chem.StereoSpecified.Unspecified]

    if stereo_status == "achiral":
        if centers or pot_stereo:
            raise ValueError(
                f"[{mid}] stereo_status='achiral' conflicts with stereogenic elements (centers={centers}, pot={len(pot_stereo)})"
            )
    elif stereo_status == "defined_stereo":
        if not (centers or pot_stereo):
            raise ValueError(f"[{mid}] stereo_status='defined_stereo' set on achiral molecule")
        if unspec_stereo:
            raise ValueError(
                f"[{mid}] stereo_status='defined_stereo' has {len(unspec_stereo)} unspecified stereogenic element(s)"
            )
    elif stereo_status == "racemate":
        if not unspec_stereo:
            raise ValueError(
                f"[{mid}] Misleading stereo_status='racemate' on molecule with zero unspecified stereocenters"
            )
        if not isinstance(item.get("review_notes"), str) or not item["review_notes"].strip():
            raise ValueError(f"[{mid}] stereo_status='racemate' requires explicit review_notes rationale")
    elif stereo_status == "unspecified_legacy":
        if review_status != "baseline_legacy":
            raise ValueError(f"[{mid}] stereo_status='unspecified_legacy' requires review_status='baseline_legacy'")

    if review_status == "baseline_legacy":
        if not isinstance(item.get("review_notes"), str) or not item["review_notes"].strip():
            raise ValueError(f"[{mid}] review_status='baseline_legacy' requires explicit review_notes")

    parent_mol = _CHEMLE_UNCHARGER.uncharge(_CHEMLE_LARGEST_FRAG.choose(std_mol))
    Chem.SanitizeMol(parent_mol)
    parent_ik = Chem.MolToInchiKey(parent_mol)
    return std_mol, std_info, parent_ik


def validate_chemle_manifest(items: Sequence[Dict[str, Any]]) -> None:
    """Validates an entire Chemle compound manifest for structural integrity and uniqueness."""
    if not isinstance(items, (list, tuple)) or len(items) == 0:
        raise ValueError("Chemle manifest must be a non-empty list of compounds")

    seen_ids: Dict[str, str] = {}
    seen_names: Dict[str, str] = {}
    seen_inchikeys: Dict[str, str] = {}
    seen_parent_iks: Dict[str, str] = {}

    for item in items:
        _, std_info, parent_ik = validate_chemle_manifest_entry(item)
        mid = item["id"].strip()
        norm_name = item["name"].strip().lower()
        ik = std_info["inchikey"]

        if mid in seen_ids:
            raise ValueError(f"Duplicate compound id {mid!r}")
        if norm_name in seen_names:
            raise ValueError(f"Duplicate compound name {item['name']!r} ({mid} vs {seen_names[norm_name]})")
        if ik in seen_inchikeys:
            raise ValueError(f"Duplicate standardized InChIKey {ik!r} ({mid} vs {seen_inchikeys[ik]})")
        if parent_ik in seen_parent_iks:
            raise ValueError(
                f"Duplicate desalted/uncharged parent InChIKey {parent_ik!r} ({mid} vs {seen_parent_iks[parent_ik]})"
            )

        seen_ids[mid] = mid
        seen_names[norm_name] = mid
        seen_inchikeys[ik] = mid
        seen_parent_iks[parent_ik] = mid

