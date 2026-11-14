"""Audit provisional intake against baseline without modifying playable data."""
import collections
import json
import sys
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from chem_engine import validate_chemle_manifest

def main():
    source = json.loads((ROOT / "data/chemle_compounds.json").read_text(encoding="utf-8"))["compounds"]
    candidates = json.loads((ROOT / "BIOWORDLE_v0.4_CANDIDATE_MANIFEST.json").read_text(encoding="utf-8"))["compounds"]
    validate_chemle_manifest(source + candidates)
    assert len(source) == 78 and len(candidates) == 42
    assert all(x["review_status"] == "pending" and not x["game_eligibility"] and not x["guess_eligibility"] for x in candidates)
    rows = source + candidates
    new = {r["id"] for r in candidates}
    group = collections.defaultdict(list)
    taut = collections.defaultdict(list)
    enum = rdMolStandardize.TautomerEnumerator()
    truncated = []
    mols = []
    for row in rows:
        mol = Chem.MolFromSmiles(row["expected_standardized_smiles"])
        if mol is None:
            raise ValueError(row["id"])
        mols.append(mol)
        group[row["expected_inchikey"].split("-")[0]].append(row["id"])
        tautomer_result = enum.Enumerate(mol)
        if str(tautomer_result.status) != "Completed":
            truncated.append({"id": row["id"], "status": str(tautomer_result.status)})
        canonical_tautomer = enum.PickCanonical(tautomer_result)
        taut[Chem.MolToSmiles(canonical_tautomer, isomericSmiles=False)].append(row["id"])
    def relevant_groups(groups):
        return [sorted(g) for g in groups.values() if len(g) > 1 and any(x in new for x in g)]
    gen = AllChem.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    fps = [gen.GetFingerprint(mol) for mol in mols]
    close = []
    for i in range(len(rows)):
        for j in range(i):
            if rows[i]["id"] not in new and rows[j]["id"] not in new:
                continue
            similarity = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            if similarity >= .8:
                close.append({"a":rows[i]["id"],"b":rows[j]["id"],"similarity":round(similarity, 4)})
    graph_review = [r["id"] for r in candidates if r.get("v04_review", {}).get("chebi_rdkit_graph_matches_pubchem") is False]
    result = {"validator":"PASS — IDs, full keys and uncharged parents under current validator",
              "base_count":len(source),"candidate_count":len(candidates),
              "candidate_review_status":"pending and not playable",
              "candidate_families":dict(collections.Counter(r["v04_family"] for r in candidates)),
              "same_inchi_connectivity_block_groups":relevant_groups(group),
              "same_rdkit_canonical_tautomer_ignoring_stereo_groups":relevant_groups(taut),
              "tautomer_enumeration_incomplete": truncated,
              "similarity_threshold":.8,"near_pairs":sorted(close,key=lambda x:(-x["similarity"],x["a"],x["b"])),
              "second_batch_graph_reviews":graph_review,
              "limitation":"Tautomer enumeration and fingerprint comparisons are screening heuristics, not exhaustive structure/biological validation."}
    (ROOT / "tools/v04_manifest_audit_results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))

if __name__ == "__main__":
    main()
