"""Fetch independent public structure records for the provisional v0.4 intake.

Run online with: python tools/v04_fetch_candidates_b.py
Results remain research data and do not confer biological or gameplay approval.
"""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NAMES = {
    "Neuroactive metabolites": ["2-phenylethylamine", "agmatine", "putrescine", "cadaverine", "spermidine", "taurine"],
    "Primary metabolism": ["L-alanine", "L-leucine", "L-serine", "L-glutamic acid", "L-lactic acid", "2-oxoglutaric acid"],
    "Microbial and related natural products": ["indole", "3-methylindole", "kojic acid", "indigo", "itaconic acid", "tropolone"],
    "Purines and nucleobases": ["guanine", "cytosine", "thymine", "uric acid", "uridine", "inosine"],
}
# Explicit source-specific forms: general-name hits need not identify a graph.
PUBCHEM_FORMS = {"indigo": 5318432}
CHEBI_FORMS = {"indole": "CHEBI:16881", "uric acid": "CHEBI:46811"}

def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "BioWordle-v04-research/1.0"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.load(r)

def main():
    results = []
    for family, names in NAMES.items():
        for name in names:
            row = {"name": name, "family": family}
            try:
                identity = ("cid/" + str(PUBCHEM_FORMS[name])) if name in PUBCHEM_FORMS else ("name/" + urllib.parse.quote(name, safe=""))
                url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/" + identity + "/property/MolecularFormula,IsomericSMILES,InChIKey,IUPACName/JSON"
                row["pubchem_request"] = url
                row["pubchem"] = get(url)["PropertyTable"]["Properties"]
                time.sleep(0.3)
                url = "https://www.ebi.ac.uk/ols4/api/search?" + urllib.parse.urlencode({"q": name, "ontology": "chebi", "exact": "on", "rows": 30})
                docs = get(url)["response"]["docs"]
                row["chebi_search"] = [{k: d.get(k) for k in ("obo_id", "label", "exact_synonyms")} for d in docs if d.get("obo_id", "").startswith("CHEBI:")]
                exact = [d for d in docs if str(d.get("label", "")).casefold() == name.casefold() and str(d.get("obo_id", "")).startswith("CHEBI:")]
                if not exact and len(docs) == 1:
                    exact = docs  # single exact-synonym result; retain the returned label for review
                if name in CHEBI_FORMS:
                    row["chebi"] = get("https://www.ebi.ac.uk/chebi/backend/api/public/compound/" + CHEBI_FORMS[name] + "/")
                elif len(exact) == 1:
                    row["chebi"] = get("https://www.ebi.ac.uk/chebi/backend/api/public/compound/" + exact[0]["obo_id"] + "/")
                else:
                    row["review_issue"] = "ChEBI exact preferred label missing or ambiguous; review candidates manually"
            except Exception as exc:
                row["error"] = str(exc)
            results.append(row)
            print(f"{len(results):02d} {name}: {'ERROR' if row.get('error') else 'review' if row.get('review_issue') else 'retrieved'}", flush=True)
            (ROOT / "tools" / "v04_candidate_b_source_records.json").write_text(json.dumps({"schema": 1, "records": results}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            time.sleep(0.2)

if __name__ == "__main__":
    main()
