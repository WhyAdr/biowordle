"""Compare represented molecular space in the checked-in and research builds.

Run after the preview: python tools/v04_space_audit.py --preview tools/v04_compiled_preview.json
"""
import argparse
import collections
import json
import statistics
from pathlib import Path

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[1]

def audit(path):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    records = data["chemle"]
    mols = [Chem.MolFromSmiles(r["canonical_smiles"]) for r in records]
    if any(m is None for m in mols):
        raise ValueError("invalid compiled SMILES")
    gen = AllChem.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    fps = [gen.GetFingerprint(m) for m in mols]
    nearest = []
    for i, fp in enumerate(fps):
        sims = DataStructs.BulkTanimotoSimilarity(fp, fps)
        sims[i] = -1
        j = max(range(len(sims)), key=sims.__getitem__)
        nearest.append({"id": records[i]["id"], "nearest_id": records[j]["id"], "similarity": round(sims[j], 5)})
    counts = collections.Counter()
    scaffolds = set()
    for m in mols:
        symbols = {a.GetSymbol() for a in m.GetAtoms()}
        for element in ("N", "O", "S", "P", "F", "Cl", "Br", "I"):
            counts[element] += int(element in symbols)
        counts["carbon_only"] += int(symbols <= {"C"})
        counts["has_stereochemistry"] += int(bool(Chem.FindMolChiralCenters(m, includeUnassigned=False)))
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=m)
        scaffolds.add(scaffold)
        counts["acyclic_scaffold"] += int(not scaffold)
    def bounds(v):
        return {"min": round(min(v), 3), "median": round(statistics.median(v), 3), "max": round(max(v), 3)}
    return {"n": len(records), "elements_or_features": dict(counts),
            "murcko_scaffolds_including_empty": len(scaffolds),
            "mw": bounds([Descriptors.MolWt(m) for m in mols]),
            "fraction_csp3": bounds([rdMolDescriptors.CalcFractionCSP3(m) for m in mols]),
            "nearest_neighbor_mean": round(statistics.mean(x["similarity"] for x in nearest), 5),
            "nearest_neighbor_median": round(statistics.median(x["similarity"] for x in nearest), 5),
            "nearest_neighbor_max": max(x["similarity"] for x in nearest),
            "nearest_neighbors": nearest}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--preview", required=True)
    p.add_argument("--output", default="tools/v04_space_results.json")
    a = p.parse_args()
    result = {"method": "RDKit Morgan radius 2, 2048 bits, chirality; nearest neighbor within each catalog; scaffold is canonical Murcko including empty", "base_78": audit(ROOT / "data/puzzles.json"), "research_120": audit(a.preview)}
    Path(a.output).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for k in ("base_78", "research_120"):
        print(k, {key: val for key, val in result[k].items() if key != "nearest_neighbors"})

if __name__ == "__main__":
    main()
