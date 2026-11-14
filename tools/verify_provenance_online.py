#!/usr/bin/env python3
"""Explicit online provenance refresh & source-comparison tool for BioWordle Seqle.

Intentionally lives outside the normal offline build path (`build_puzzles.py`).
Fetches UniProtKB FASTA records (`rest.uniprot.org`) and RCSB PDB core/entity
records (`data.rcsb.org`) to cross-check every entry in `SEQLE_REFERENCE_REGISTRY`.
"""

from __future__ import annotations

import json
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Tuple

ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

from chem_engine import SEQLE_REFERENCE_REGISTRY  # noqa: E402


def _fetch_text(url: str) -> Tuple[str, Dict[str, str]]:
    req = urllib.request.Request(url, headers={"User-Agent": "BioWordle-Provenance-Verifier/0.2.2"})
    with urllib.request.urlopen(req, timeout=20) as resp:
        headers = {k.lower(): v for k, v in resp.headers.items()}
        return resp.read().decode("utf-8"), headers


def _fetch_uniprot_fasta(accession: str) -> Tuple[str, str]:
    fasta, hdrs = _fetch_text(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta")
    seq = "".join(
        line.strip() for line in fasta.splitlines() if line and not line.startswith(">")
    )
    release = hdrs.get("x-uniprot-release") or hdrs.get("last-modified") or "unknown"
    return seq, release


def _fetch_rcsb_entity(pdb_id: str, entity_id: str) -> Dict[str, Any]:
    raw, _ = _fetch_text(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}")
    return json.loads(raw)


def verify_all(verbose: bool = True) -> int:
    releases: set[str] = set()
    failures: List[str] = []

    for motif, ref in sorted(SEQLE_REFERENCE_REGISTRY.items()):
        status = ref.get("status")
        if status == "verified_uniprot":
            acc = str(ref["accession"])
            try:
                full_seq, rel = _fetch_uniprot_fasta(acc)
                releases.add(rel)
                p0, p1 = ref["parent_residue_range"]
                m0, m1 = ref["motif_residue_range"]
                obs_parent = full_seq[p0 - 1 : p1]
                obs_motif = full_seq[m0 - 1 : m1]
                if obs_parent != ref["verified_parent_sequence"]:
                    failures.append(
                        f"[{motif}] UniProt {acc} [{p0}..{p1}]={obs_parent!r} != {ref['verified_parent_sequence']!r}"
                    )
                if obs_motif != motif:
                    failures.append(
                        f"[{motif}] UniProt {acc} motif [{m0}..{m1}]={obs_motif!r} != {motif!r}"
                    )
                if verbose:
                    print(f"OK UniProt | {motif} | {acc} | parent={p0}-{p1} motif={m0}-{m1} | release={rel}")
            except (urllib.error.URLError, KeyError, IndexError, ValueError) as exc:
                failures.append(f"[{motif}] UniProt {acc} fetch/check error: {exc}")

        for cit in ref.get("citations") or []:
            if cit.get("source") != "PDB":
                continue
            pdb_id = str(cit["id"])
            role = cit.get("role")
            ent_id = cit.get("entity_id")
            if not ent_id:
                if role == "exact_structure":
                    failures.append(f"[{motif}] exact_structure PDB:{pdb_id} missing entity_id")
                elif verbose:
                    print(f"OK Context | {motif} | PDB:{pdb_id} ({role})")
                continue
            try:
                ent = _fetch_rcsb_entity(pdb_id, str(ent_id))
                poly = ent.get("entity_poly") or {}
                can_seq = str(poly.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "").strip()
                chains = list(poly.get("pdbx_strand_id").split(",")) if isinstance(poly.get("pdbx_strand_id"), str) else []
                chains = [c.strip() for c in chains if c.strip()]
                ups = list((ent.get("rcsb_polymer_entity_container_identifiers") or {}).get("uniprot_ids") or [])

                if role == "exact_structure":
                    r0, r1 = cit["residue_range"]
                    obs = can_seq[r0 - 1 : r1]
                    if obs != motif:
                        failures.append(
                            f"[{motif}] PDB:{pdb_id}:{ent_id} [{r0}..{r1}]={obs!r} != {motif!r}"
                        )
                    dep = str(cit.get("deposited_sequence") or "")
                    if dep not in can_seq:
                        failures.append(
                            f"[{motif}] PDB:{pdb_id}:{ent_id} deposited_sequence {dep!r} not in entity sequence"
                        )
                    if sorted(cit.get("chains") or []) != sorted(chains):
                        failures.append(
                            f"[{motif}] PDB:{pdb_id}:{ent_id} chains {cit.get('chains')} != RCSB {chains}"
                        )
                    expected_up = cit.get("mapped_uniprot")
                    if expected_up is not None and expected_up not in ups:
                        failures.append(
                            f"[{motif}] PDB:{pdb_id}:{ent_id} expected mapped_uniprot {expected_up!r} not in {ups}"
                        )
                    if verbose:
                        print(f"OK PDB     | {motif} | {pdb_id}:{ent_id} chains={chains} [{r0}-{r1}]={obs} ups={ups}")
                elif role == "motif_family_context":
                    var_seq = cit.get("variant_sequence")
                    if var_seq and var_seq not in can_seq:
                        failures.append(
                            f"[{motif}] PDB:{pdb_id}:{ent_id} context variant_sequence {var_seq!r} not found in {can_seq[:40]!r}..."
                        )
                    if verbose:
                        print(f"OK Context | {motif} | {pdb_id}:{ent_id} chains={chains} ups={ups}")
            except (urllib.error.URLError, KeyError, ValueError) as exc:
                failures.append(f"[{motif}] PDB:{pdb_id}:{ent_id} fetch/check error: {exc}")

    print(f"\n[Online Provenance Verification] Releases observed: {sorted(releases)}")
    if failures:
        print("FAILURES:", file=sys.stderr)
        for f in failures:
            print("  -", f, file=sys.stderr)
        return 1
    print(f"[Online Provenance Verification] All {len(SEQLE_REFERENCE_REGISTRY)} registry entries verified cleanly.")
    return 0


if __name__ == "__main__":
    sys.exit(verify_all())

