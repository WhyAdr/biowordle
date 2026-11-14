"""Standalone Online Manifest Verifier for BioWordle Chemle Manifest.

Opt-in network audit (not run during offline build or unit tests): queries
PubChem PUG REST and EBI OLS4 ChEBI for every compound in
`data/chemle_compounds.json`, checking that recorded PubChem CIDs, ChEBI IDs,
and InChIKeys resolve consistently.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Dict, List

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from chem_engine import validate_chemle_manifest


def _fetch_json(url: str, timeout: float = 15.0) -> Dict[str, Any]:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "BioWordle-Manifest-Verifier/0.3"},
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def verify_manifest_online(manifest_path: Path, delay_s: float = 0.22) -> int:
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    compounds: List[Dict[str, Any]] = payload["compounds"]
    validate_chemle_manifest(compounds)
    print(f"[Offline Validator PASS] {len(compounds)} manifest compounds verified locally.")

    failures: List[str] = []
    for item in compounds:
        mid = item["id"]
        cid = int(item["structure_refs"]["pubchem_cid"])
        chebi_id = item["structure_refs"]["chebi_id"]
        exp_ik = item["expected_inchikey"]
        review_status = item["review_status"]

        pc_url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
            "/property/Title,InChIKey,MolecularFormula/JSON"
        )
        try:
            pc_data = _fetch_json(pc_url)
            props = pc_data.get("PropertyTable", {}).get("Properties", [{}])[0]
            pc_ik = props.get("InChIKey")
        except Exception as exc:
            failures.append(f"{mid}: PubChem CID {cid} lookup failed: {exc}")
            continue

        chebi_url = (
            "https://www.ebi.ac.uk/ols4/api/ontologies/chebi/terms?"
            + urllib.parse.urlencode({"obo_id": chebi_id})
        )
        try:
            chebi_data = _fetch_json(chebi_url)
            terms = chebi_data.get("_embedded", {}).get("terms", [])
            chebi_ok = any(t.get("obo_id") == chebi_id for t in terms)
        except Exception as exc:
            failures.append(f"{mid}: ChEBI {chebi_id} lookup failed: {exc}")
            continue

        ik_match = pc_ik == exp_ik
        if review_status == "verified" and not ik_match:
            failures.append(
                f"{mid}: verified compound InChIKey {exp_ik} != PubChem CID {cid} InChIKey {pc_ik}"
            )
        if not chebi_ok:
            failures.append(f"{mid}: ChEBI {chebi_id} not resolved in OLS4")

        tag = "OK" if ik_match else f"LEGACY_DIFF (PubChem={pc_ik})"
        print(f"  [{review_status:15s}] {mid:22s} | CID:{cid:<10d} | {chebi_id:13s} | {tag}")
        if delay_s > 0:
            time.sleep(delay_s)

    if failures:
        print("\n[Online Manifest Failures]")
        for f in failures:
            print("  -", f)
        return 1

    print(f"[Online Manifest PASS] All {len(compounds)} PubChem/ChEBI references resolved.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify data/chemle_compounds.json against live PubChem & ChEBI.")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "data" / "chemle_compounds.json",
        help="Path to chemle_compounds.json",
    )
    parser.add_argument("--delay", type=float, default=0.2, help="Inter-request delay in seconds")
    args = parser.parse_args()
    return verify_manifest_online(args.manifest, delay_s=args.delay)


if __name__ == "__main__":
    raise SystemExit(main())

