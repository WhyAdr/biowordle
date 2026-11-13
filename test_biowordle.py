"""Automated Scientific, Game-Theoretic & Runtime Audit Suite for BioWordle v0.2 ("Restraint").

Covers all requirements from Sol's technical review and `BIOWORDLE_PATCH.md`:
1. Compiler independence & schema v0.2.0 metadata
2. Molecular identity, isomer separation, chiral Morgan FP, and ACS-1996 normalized SVG hashes
3. Chiral Morgan (ECFP4, r=2, 2048b) Tanimoto matrix symmetry & self-identity
4. Descriptor validity, honest `CNS_MPO_4prop` naming, and acyclic scaffold fallback
5. P0-1: Opening-clue k-anonymity gate (k_min >= 4 in Ranked, >= 3 in Learn) & non-leakage
6. P0-2: Greedy max-entropy solver-bot calibration gate (Ranked & Learn profiles)
7. P0-4: Seqle motif schema validator, scientific corrections, and negative tests (KDELR, CVIMK)
8. P1-3: Order-independent canonical dataset_version_hash invariance & sensitivity
9. P1-5: Non-overlapping split carbonyl SMARTS verification
10. P2-3: Pre-reveal hint text leakage lint & negative stem/3-letter AA tests
11. P1-1 / P0-3 / P1-2 / P1-4: Cross-language Python <-> JS golden vectors (4,392 pairs) & DOM/state suite
"""

from __future__ import annotations

import copy
import json
import math
import os
import random
import subprocess
import unittest
from pathlib import Path

from rdkit import Chem, DataStructs

import build_puzzles
import chem_engine
from tools import solver_bot


BASE_DIR = os.path.dirname(os.path.abspath(__file__))


class TestBioWordleAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        json_path = os.path.join(BASE_DIR, "data", "puzzles.json")
        with open(json_path, "r", encoding="utf-8") as f:
            cls.data = json.load(f)
        spec_path = os.path.join(BASE_DIR, "data", "feedback_spec.json")
        with open(spec_path, "r", encoding="utf-8") as f:
            cls.feedback_spec = json.load(f)

    def test_01_compiler_independence_and_schema(self) -> None:
        """Verify neither chem_engine.py nor build_puzzles.py references external .gemini paths."""
        for fname in ("chem_engine.py", "build_puzzles.py"):
            with open(os.path.join(BASE_DIR, fname), "r", encoding="utf-8") as f:
                src = f.read()
            self.assertNotIn(".gemini", src, f"{fname} must not reference external .gemini paths")
            self.assertNotIn("skill_paths", src, f"{fname} must not import external skill_paths")
        self.assertEqual(self.data["metadata"]["schema_version"], "0.2.0")

    def test_02_molecular_identity_stereochem_and_acs1996_svg(self) -> None:
        """Verify unique IDs, InChIKeys, canonical SMILES, preserved stereocenters, and normalized SVG hashes."""
        chemle = self.data["chemle"]
        self.assertEqual(len(chemle), 36)

        ids = [c["id"] for c in chemle]
        inchikeys = [c["inchikey"] for c in chemle]
        can_smiles = [c["canonical_smiles"] for c in chemle]

        self.assertEqual(len(set(ids)), 36, "Duplicate compound IDs found")
        self.assertEqual(len(set(inchikeys)), 36, "Duplicate InChIKeys after standardization")
        self.assertEqual(len(set(can_smiles)), 36, "Duplicate canonical SMILES after standardization")
        self.assertEqual(ids, sorted(ids), "Chemle catalog must be canonically sorted by id")

        by_id = {c["id"]: c for c in chemle}
        pca = by_id["protocatechuic_acid"]
        gen = by_id["gentisic_acid"]
        self.assertEqual(pca["Formula"], "C7H6O4")
        self.assertEqual(gen["Formula"], "C7H6O4")
        self.assertNotEqual(pca["inchikey"], gen["inchikey"])
        self.assertLess(self.data["chemle_similarity"]["protocatechuic_acid"]["gentisic_acid"], 100)

        for chiral_id in ("adrenaline", "tryptophan", "ascorbic_acid", "biotin", "quinine", "artemisinin"):
            self.assertGreater(by_id[chiral_id]["ChiralCenters"], 0)
            self.assertIn("@", by_id[chiral_id]["canonical_smiles"])

        for c in chemle:
            self.assertEqual(len(c["svg_sha256"]), 12)
            self.assertIn('width="100%" height="100%"', c["svg_full"])

        r_adr = Chem.MolFromSmiles("CNC[C@H](O)c1ccc(O)c(O)c1")
        s_adr = Chem.MolFromSmiles("CNC[C@@H](O)c1ccc(O)c(O)c1")
        fp_r = chem_engine.compute_chiral_morgan_fp(r_adr)
        fp_s = chem_engine.compute_chiral_morgan_fp(s_adr)
        self.assertLess(DataStructs.TanimotoSimilarity(fp_r, fp_s), 1.0)

    def test_03_tanimoto_matrix_symmetry_and_spec(self) -> None:
        """Verify fingerprint specification and 36x36 Tanimoto matrix properties."""
        fp_spec = self.data["metadata"]["fingerprint_spec"]
        self.assertEqual(fp_spec["radius"], 2)
        self.assertEqual(fp_spec["fpSize"], 2048)
        self.assertTrue(fp_spec["includeChirality"])
        self.assertTrue(fp_spec["useBondTypes"])

        sim = self.data["chemle_similarity"]
        ids = [c["id"] for c in self.data["chemle"]]
        for i in ids:
            self.assertEqual(sim[i][i], 100)
            for j in ids:
                self.assertEqual(sim[i][j], sim[j][i])
                self.assertTrue(0 <= sim[i][j] <= 100)
                if i != j:
                    self.assertLess(sim[i][j], 100)

    def test_04_descriptors_and_acyclic_scaffold(self) -> None:
        """Verify no NaN/Inf descriptors, honest CNS_MPO_4prop naming, and acyclic scaffold fallback."""
        numeric_keys = (
            "MW",
            "exact_mass",
            "cLogP",
            "TPSA",
            "HBD",
            "HBA",
            "RotB",
            "Fsp3",
            "HeavyAtoms",
            "Rings",
            "AromaticRings",
            "AliphaticRings",
            "ChiralCenters",
            "QED",
            "SA_Score",
            "CNS_MPO_4prop",
        )
        for c in self.data["chemle"]:
            for k in numeric_keys:
                val = c[k]
                self.assertIsNotNone(val, f"{c['id']}.{k} is None")
                self.assertFalse(math.isnan(float(val)) or math.isinf(float(val)))
            self.assertTrue(0.0 <= c["CNS_MPO_4prop"] <= 4.0)
            self.assertIn("4-parameter", c["CNS_MPO_Method"])

        by_id = {c["id"]: c for c in self.data["chemle"]}
        citric = by_id["citric_acid"]
        self.assertEqual(citric["Rings"], 0)
        self.assertEqual(citric["scaffold_mode"], "carbon_skeleton")
        self.assertNotEqual(citric["svg_scaffold"], citric["svg_full"])

    def test_05_p0_1_opening_clue_k_anonymity(self) -> None:
        """Verify Level-0 opening clues satisfy k_min >= 4 (Ranked) and k_min >= 3 (Learn) without leaking."""
        k_audit = self.data["metadata"]["k_anonymity_audit"]
        self.assertGreaterEqual(k_audit["chemle_ranked"]["k_min"], 4)
        self.assertGreaterEqual(k_audit["chemle_learn"]["k_min"], 3)
        self.assertGreaterEqual(k_audit["seqle_ranked"]["k_min"], 4)
        self.assertGreaterEqual(k_audit["seqle_learn"]["k_min"], 3)

        for c in self.data["chemle"]:
            for prof in ("ranked", "learn"):
                lvl0 = c["clue_ladder"][prof][0]["text"]
                self.assertNotIn(str(c["MW"]), lvl0)
                self.assertNotIn(c["Category"], lvl0)
                self.assertNotIn("bits", lvl0.lower())

        for m in self.data["seqle"]:
            for prof in ("ranked", "learn"):
                lvl0 = m["clue_ladder"][prof][0]["text"]
                self.assertNotIn(m["family"], lvl0)
                self.assertNotIn(m["seq"], lvl0)

    def test_06_p0_2_solver_bot_calibration_gate(self) -> None:
        """Run the greedy max-entropy solver across Ranked and Learn profiles and enforce difficulty gates."""
        report = solver_bot.run_all_gates(
            Path(BASE_DIR) / "data" / "puzzles.json",
            Path(BASE_DIR) / "data" / "feedback_spec.json",
        )
        r_stats = report["chemle"]["ranked"]
        self.assertTrue(3.5 <= r_stats["mean"] <= 4.5)
        self.assertLess(r_stats["pct_solved_le2"], 5.0)
        self.assertGreaterEqual(r_stats["pct_solved_le6"], 98.0)
        self.assertGreaterEqual(r_stats["min_opening_expected_remaining"], 2.5)

        l_stats = report["chemle"]["learn"]
        self.assertLessEqual(l_stats["mean"], 3.5)
        self.assertLess(l_stats["pct_solved_le2"], 25.0)
        self.assertGreaterEqual(l_stats["pct_solved_le6"], 98.0)
        self.assertGreaterEqual(l_stats["min_opening_expected_remaining"], 2.5)

    def test_07_p0_4_seqle_motif_validator_and_negative_cases(self) -> None:
        """Verify all 30 curated motifs pass P0-4 schema validation and test negative rejection of KDELR/CVIMK."""
        aa = self.data["aa_properties"]
        seqle = self.data["seqle"]
        self.assertEqual(len(seqle), 30)
        by_seq = {m["seq"]: m for m in seqle}

        # Verify scientific replacements from P0-4 table
        for expected_seq in (
            "EKDEL",
            "KCVIM",
            "GRGDS",
            "DEADE",
            "GVGKS",
            "DSGLD",
            "KENIP",
            "KKKRK",
            "FGGFT",
            "SYSGY",
            "GHHPH",
            "INFEK",
            "CCHHC",
            "LDLAS",
        ):
            self.assertIn(expected_seq, by_seq, f"Missing audited motif {expected_seq}")

        # Verify abstractions CCHHC and LDLAS are excluded from Ranked mode
        for abs_seq in ("CCHHC", "LDLAS"):
            self.assertEqual(by_seq[abs_seq]["kind"], "abstraction")
            self.assertNotIn("ranked", by_seq[abs_seq]["game_eligibility"])

        # Negative test 1: KDELR with terminal='C' but non-C-terminal offset in parent sequence must fail
        bad_kdelr = {
            "seq": "KDELR",
            "kind": "exact",
            "consensus": "K-D-E-L-R",
            "parent_sequence": "TAEKDELRGG",
            "window_offset": 3,
            "terminal": "C",
            "game_eligibility": ["learn", "ranked"],
            "exemplar_source": "UniProt:P11021",
        }
        with self.assertRaises(ValueError):
            chem_engine.validate_seqle_motif_record(bad_kdelr, aa)

        # Negative test 2: CVIMK with terminal='C' when C-terminal CAAX is -CVIM must fail
        bad_cvimk = {
            "seq": "CVIMK",
            "kind": "fragment",
            "consensus": "C-V-I-M-K",
            "parent_sequence": "KKKKKKSKTKCVIM",
            "window_offset": 10,
            "terminal": "C",
            "game_eligibility": ["learn", "ranked"],
            "exemplar_source": "UniProt:P01116",
        }
        with self.assertRaises(ValueError):
            chem_engine.validate_seqle_motif_record(bad_cvimk, aa)

        # Negative test 3: abstraction in ranked pool must fail
        bad_abstraction = copy.deepcopy(by_seq["CCHHC"])
        bad_abstraction["game_eligibility"] = ["learn", "ranked"]
        with self.assertRaises(ValueError):
            chem_engine.validate_seqle_motif_record(bad_abstraction, aa)

    def test_08_p1_3_order_independent_canonical_hash(self) -> None:
        """Verify shuffling catalog order leaves dataset_version_hash invariant while descriptor edits change it."""
        chem_copy = copy.deepcopy(self.data["chemle"])
        seq_copy = copy.deepcopy(self.data["seqle"])
        rng = random.Random(42)
        rng.shuffle(chem_copy)
        rng.shuffle(seq_copy)

        h_orig = self.data["metadata"]["dataset_version_hash"]
        h_shuffled = build_puzzles.compute_canonical_dataset_hash(
            chem_copy, seq_copy, self.feedback_spec
        )
        self.assertEqual(h_orig, h_shuffled, "dataset_version_hash must be invariant to catalog input ordering")

        chem_mutated = copy.deepcopy(self.data["chemle"])
        chem_mutated[0]["MW"] = round(chem_mutated[0]["MW"] + 1.0, 1)
        h_mutated = build_puzzles.compute_canonical_dataset_hash(
            chem_mutated, self.data["seqle"], self.feedback_spec
        )
        self.assertNotEqual(h_orig, h_mutated, "dataset_version_hash must change when a gameplay descriptor changes")

    def test_09_p1_5_split_carbonyl_smarts(self) -> None:
        """Verify non-overlapping split carbonyl SMARTS distinguish ketone, aldehyde, acid, ester, and amide."""
        cases = [
            ("CC(=O)C", {"Ketone (C-C(=O)-C)"}, {"Aldehyde (-CHO)", "Carboxylic Acid (-COOH)", "Ester / Lactone", "Amide / Lactam"}),
            ("CC=O", {"Aldehyde (-CHO)"}, {"Ketone (C-C(=O)-C)", "Carboxylic Acid (-COOH)", "Ester / Lactone", "Amide / Lactam"}),
            ("CC(=O)O", {"Carboxylic Acid (-COOH)"}, {"Ketone (C-C(=O)-C)", "Aldehyde (-CHO)", "Ester / Lactone", "Amide / Lactam"}),
            ("CC(=O)OC", {"Ester / Lactone"}, {"Ketone (C-C(=O)-C)", "Aldehyde (-CHO)", "Carboxylic Acid (-COOH)", "Amide / Lactam"}),
            ("CC(=O)NC", {"Amide / Lactam"}, {"Ketone (C-C(=O)-C)", "Aldehyde (-CHO)", "Carboxylic Acid (-COOH)", "Ester / Lactone"}),
            ("CC(=O)Oc1ccccc1C(=O)O", {"Ester / Lactone", "Carboxylic Acid (-COOH)"}, {"Ketone (C-C(=O)-C)", "Aldehyde (-CHO)", "Amide / Lactam"}),
        ]
        for smi, must_have, must_not_have in cases:
            mol = Chem.MolFromSmiles(smi)
            fgs = set(chem_engine.detect_functional_groups(mol))
            self.assertNotIn("Ketone / Carbonyl (C=O)", fgs)
            for req in must_have:
                self.assertIn(req, fgs, f"{smi} failed to match {req}; got {fgs}")
            for excl in must_not_have:
                self.assertNotIn(excl, fgs, f"{smi} falsely matched {excl}; got {fgs}")

    def test_10_p2_3_hint_leakage_audit(self) -> None:
        """Verify all catalog items pass hint leakage lint and negative cases raise AssertionError."""
        for c in self.data["chemle"]:
            chem_engine.audit_hint_text_leakage(c, ("hint_role", "hint_origin"), "id")
        for m in self.data["seqle"]:
            chem_engine.audit_hint_text_leakage(m, ("hint_role", "hint_detail"), "seq")

        with self.assertRaises(AssertionError):
            chem_engine.audit_hint_text_leakage(
                {"id": "penicillin_g", "name": "Penicillin G", "hint_role": "Classic penicillin antibiotic."},
                ("hint_role",),
                "id",
            )
        with self.assertRaises(AssertionError):
            chem_engine.audit_hint_text_leakage(
                {"seq": "KKKRK", "title": "SV40 NLS", "hint_detail": "Contains Lys-Lys-Lys-Arg-Lys."},
                ("hint_detail",),
                "seq",
            )

    def test_11_cross_language_golden_vectors_and_dom_suite(self) -> None:
        """Execute Node.js cross-language golden vector & DOM/persistence/lockout test suite."""
        proc = subprocess.run(
            ["node", os.path.join(BASE_DIR, "tests", "verify_golden_and_dom.mjs")],
            cwd=BASE_DIR,
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(
            proc.returncode,
            0,
            f"Node golden/DOM suite failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}",
        )
        self.assertIn("2592 Chemle vectors + 1800 Seqle vectors", proc.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
