"""Prepare a non-playable research intake from the saved source responses.

This is an intake, not a scientific approval. Missing source matches and
claim-level evidence remain visible in each row's v04_review field.
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from chem_engine import standardize_and_audit

CATEGORIES = {
    "Neuroactive metabolites": "Neurotransmitter",
    "Primary metabolism": "Primary Metabolite",
    "Microbial and related natural products": "Microbial Natural Product",
    "Purines and nucleobases": "Primary Metabolite",
}
# Narrow study leads, not approved gameplay claims. Other rows intentionally
# carry no invented biological mechanism or pathway.
LEADS = {
    "2-phenylethylamine": ("39596410", "Selected Staphylococcus skin isolates produced phenylethylamine in the study; SadA expression reproduced the reported product profile."),
    "agmatine": ("37878710", "In Synechocystis sp. PCC 6803, isotope tracing, genetics and biochemistry support agmatine as a precursor in one spermidine pathway."),
    "spermidine": ("37878710", "A Synechocystis sp. PCC 6803 pathway to spermidine via carboxyaminopropylagmatine was investigated using isotope tracing and enzyme assays."),
    "indole": ("7487051", "A mixed pig-faecal bacterial culture converted added tryptophan to indole under the reported conditions."),
    "3-methylindole": ("30310076", "An indoleacetate decarboxylase was biochemically characterized as producing skatole in certain anaerobic bacteria."),
    "kojic acid": ("20849972", "Aspergillus oryzae gene-disruption experiments identified two genes required for kojic acid production."),
    "itaconic acid": ("32162091", "In Aspergillus pseudoterreus, deletion experiments identified cluster genes needed for itaconic acid production."),
}

def main():
    source = json.loads((ROOT / "tools/v04_candidate_b_source_records.json").read_text(encoding="utf-8"))
    existing = json.loads((ROOT / "data/chemle_compounds.json").read_text(encoding="utf-8"))["compounds"]
    existing_ids = {x["id"] for x in existing}
    existing_keys = {x["expected_inchikey"] for x in existing}
    out = []
    for src in source["records"]:
        name = src["name"]
        props = src.get("pubchem", [])
        if len(props) != 1:
            raise ValueError(f"{name}: expected exactly one PubChem property record, got {len(props)}")
        pc = props[0]
        smiles = pc.get("SMILES") or pc.get("IsomericSMILES")
        if not smiles:
            raise ValueError(f"{name}: no PubChem isomeric SMILES")
        mol, audit = standardize_and_audit(smiles)
        key = audit["inchikey"]
        formula = rdMolDescriptors.CalcMolFormula(mol)
        identifier = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
        if identifier in existing_ids or key in existing_keys:
            raise ValueError(f"duplicate baseline ID or InChIKey: {name}")
        chebi = src.get("chebi", {})
        chebi_id = chebi.get("chebi_accession")
        if not chebi_id:
            raise ValueError(f"{name}: no ChEBI identifier found")
        ckey = (chebi.get("default_structure") or {}).get("standard_inchi_key")
        chebi_smiles = (chebi.get("default_structure") or {}).get("smiles")
        chebi_mol = Chem.MolFromSmiles(chebi_smiles) if chebi_smiles else None
        source_graph = Chem.MolFromSmiles(smiles)
        same_graph = bool(chebi_mol and source_graph and
                          Chem.MolToSmiles(chebi_mol, isomericSmiles=True) ==
                          Chem.MolToSmiles(source_graph, isomericSmiles=True))
        potential = list(Chem.FindPotentialStereo(mol))
        unspec = [x for x in potential if x.specified == Chem.StereoSpecified.Unspecified]
        if unspec:
            raise ValueError(f"{name}: unspecified potential stereochemistry; select a supported form before preview")
        status = "defined_stereo" if potential else "achiral"
        issues = []
        if key != pc.get("InChIKey") or formula != pc.get("MolecularFormula"):
            issues.append("Project normalization differs from PubChem key or formula; investigate.")
        if not ckey:
            issues.append("ChEBI source graph/key not yet retrieved and checked.")
        elif ckey != key:
            issues.append("ChEBI and selected PubChem InChIKeys differ; classify graph, charge, stereo or tautomer before approval.")
        if chebi_smiles and chebi_mol is None:
            issues.append("ChEBI displayed SMILES does not parse in RDKit; independently resolve its graph.")
        elif chebi_mol and not same_graph:
            issues.append("ChEBI and PubChem isomeric graphs differ despite any key match; review tautomer and bond states.")
        issues.append("No claim-bound primary biological mechanism or biosynthesis reference reviewed.")
        cid = pc["CID"]
        row = {
            "id": identifier, "name": name, "aliases": [], "Category": CATEGORIES[src["family"]],
            "primary_class": "v04_research_candidate", "role_tags": ["research_intake"],
            "input_isomeric_smiles": smiles, "expected_standardized_smiles": audit["canonical_smiles"],
            "expected_inchikey": key, "expected_formula": formula, "stereo_status": status,
            "protonation_policy": "neutral_canonical", "game_eligibility": [], "guess_eligibility": [],
            "structure_refs": {"pubchem_cid": cid, "chebi_id": chebi_id,
                               "pubchem_matched_inchikey": key == pc.get("InChIKey")},
            "evidence_refs": [
                {"source": "PubChem", "id": f"CID:{cid}", "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}", "role": "structure_reference"},
                {"source": "ChEBI", "id": chebi_id, "url": f"https://www.ebi.ac.uk/chebi/{chebi_id}", "role": "independent_structure_reference"}],
            "hint_role": "Compare the molecular properties of this example.",
            "hint_origin": "Consult reviewed literature for its biological context.",
            "structural_note": "Structure and biological context require review before use.",
            "review_status": "pending", "review_notes": "Research intake only; no biological or gameplay approval.",
            "v04_family": src["family"],
            "v04_review": {"status": "pending", "pubchem_key": pc.get("InChIKey"),
                           "chebi_key": ckey, "pubchem_smiles": smiles, "chebi_smiles": chebi_smiles,
                           "chebi_rdkit_graph_matches_pubchem": same_graph,
                           "issues": issues},
        }
        if name in LEADS:
            pmid, claim = LEADS[name]
            url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
            row["evidence_refs"].append({"source": "PubMed", "id": f"PMID:{pmid}", "url": url, "role": "candidate_biological_context"})
            row["v04_review"]["biological_evidence"] = {"claim": claim, "source": url, "status": "abstract-level lead; independent full-study review required"}
        out.append(row)
    if len(out) != 24:
        raise ValueError(f"Expected 24 source records, got {len(out)}")
    (ROOT / "tools/v04_candidates_b.json").write_text(json.dumps({"schema_version": "v04-research-intake", "compounds": out}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Prepared {len(out)} pending rows; none is authorized for play.")

if __name__ == "__main__":
    main()
