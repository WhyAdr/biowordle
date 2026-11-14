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
        self.assertEqual(self.data["metadata"]["schema_version"], "0.2.2")

    def test_02_molecular_identity_stereochem_and_acs1996_svg(self) -> None:
        """Verify unique IDs, InChIKeys, canonical SMILES, preserved stereocenters, and normalized SVG hashes."""
        chemle = self.data["chemle"]
        self.assertEqual(len(chemle), 78)
        self.assertEqual(self.data["metadata"]["chemle_ranked_count"], 36)
        self.assertEqual(self.data["metadata"]["chemle_ranked_guess_count"], 36)
        self.assertEqual(self.data["metadata"]["chemle_learn_count"], 78)
        self.assertEqual(self.data["metadata"]["chemle_learn_guess_count"], 78)

        ids = [c["id"] for c in chemle]
        inchikeys = [c["inchikey"] for c in chemle]
        can_smiles = [c["canonical_smiles"] for c in chemle]

        self.assertEqual(len(set(ids)), 78, "Duplicate compound IDs found")
        self.assertEqual(len(set(inchikeys)), 78, "Duplicate InChIKeys after standardization")
        self.assertEqual(len(set(can_smiles)), 78, "Duplicate canonical SMILES after standardization")
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
                if prof not in c.get("game_eligibility", []):
                    continue
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
        """Run clue-conditioned Chemle & Seqle greedy solvers, enforce gates, and verify pre-tuning YIGSR 28/28 collapse."""
        report = solver_bot.run_all_gates(
            Path(BASE_DIR) / "data" / "puzzles.json",
            Path(BASE_DIR) / "data" / "feedback_spec.json",
        )
        r_stats = report["chemle"]["ranked"]
        self.assertTrue(r_stats["conditioned_on_clues"])
        self.assertTrue(3.0 <= r_stats["mean"] <= 4.3)
        self.assertLess(r_stats["pct_solved_le2"], self.feedback_spec["profiles"]["ranked"]["solver_gate"]["max_le2_pct"])
        self.assertGreaterEqual(r_stats["pct_solved_le6"], 98.0)
        self.assertGreaterEqual(r_stats["min_opening_expected_remaining"], 2.5)

        l_stats = report["chemle"]["learn"]
        self.assertTrue(l_stats["conditioned_on_clues"])
        self.assertLessEqual(l_stats["mean"], 3.5)
        self.assertLess(l_stats["pct_solved_le2"], self.feedback_spec["profiles"]["learn"]["solver_gate"]["max_le2_pct"])
        self.assertGreaterEqual(l_stats["pct_solved_le6"], 98.0)
        self.assertGreaterEqual(l_stats["min_opening_expected_remaining"], 1.0)

        sr_stats = report["seqle"]["ranked"]
        sl_stats = report["seqle"]["learn"]
        self.assertEqual(sr_stats["eligible_count"], 28)
        self.assertEqual(sl_stats["eligible_count"], 30)
        self.assertLessEqual(
            sr_stats["max_single_probe_distinct_targets"],
            self.feedback_spec["profiles"]["ranked"]["seqle_solver_gate"]["max_single_probe_distinct_targets"],
        )
        self.assertLess(sr_stats["yigsr_distinct_signatures"], 28)

        # Reproduce pre-tuning v0.2.0 Ranked Seqle feedback against corrected NPLYK and legacy NPLYR catalogs
        raw_ranked_spec = copy.deepcopy(self.feedback_spec["profiles"]["ranked"])
        raw_ranked_spec["show_positional_class_match"] = True
        raw_ranked_spec["seqle_strip_fields"] = [
            {
                "name": "net_charge_ph74",
                "label": "Q (pH 7.4)",
                "decimals": 2,
                "kind": "delta",
                "exact_tol": 0.4,
                "close_tol": 1.8,
                "close_directional": False,
                "far_directional": True,
                "shown": True,
            },
            {
                "name": "gravy",
                "label": "GRAVY",
                "decimals": 2,
                "kind": "delta",
                "exact_tol": 0.6,
                "close_tol": 1.8,
                "close_directional": False,
                "far_directional": True,
                "shown": True,
            },
        ]
        raw_seq_stats = chem_engine.run_seqle_greedy_solver(
            self.data["seqle"],
            self.data["aa_properties"],
            raw_ranked_spec,
            profile_id="ranked",
            condition_on_clues=False,
        )
        self.assertEqual(raw_seq_stats["yigsr_distinct_signatures"], 27)
        self.assertEqual(raw_seq_stats["max_single_probe_distinct_targets"], 28)
        self.assertEqual(raw_seq_stats["max"], 2)

        # Confirm exact v0.2.0 NPLYR baseline reproduction: YIGSR splits all 28/28 targets
        legacy_seqle = copy.deepcopy(self.data["seqle"])
        for rec in legacy_seqle:
            if rec["seq"] == "NPLYK":
                bio = chem_engine.compute_peptide_biochem("NPLYR", self.data["aa_properties"])
                rec["seq"] = "NPLYR"
                rec.update(bio)
                break
        legacy_raw_stats = chem_engine.run_seqle_greedy_solver(
            legacy_seqle,
            self.data["aa_properties"],
            raw_ranked_spec,
            profile_id="ranked",
            condition_on_clues=False,
        )
        self.assertEqual(legacy_raw_stats["yigsr_distinct_signatures"], 28)
        self.assertEqual(legacy_raw_stats["max"], 2)

    def test_07_p0_4_seqle_motif_validator_and_negative_cases(self) -> None:
        """Verify all 30 curated motifs pass P0-4 & UniProt provenance validation and test negative rejections."""
        aa = self.data["aa_properties"]
        seqle = self.data["seqle"]
        self.assertEqual(len(seqle), 30)
        by_seq = {m["seq"]: m for m in seqle}

        # Verify scientific replacements & v0.2.1 NPLYK (ITGB3 P05106) correction
        self.assertIn("NPLYK", by_seq)
        self.assertNotIn("NPLYR", by_seq)
        self.assertIn("UniProt:P05106", by_seq["NPLYK"]["exemplar_source"])
        self.assertEqual(by_seq["NPLYK"]["parent_sequence"], "AKWDTANNPLYKEATSTFTNITYRGT")
        self.assertEqual(by_seq["NPLYK"]["reference_provenance"]["status"], "verified_uniprot")

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
            "NPLYK",
            "GHHPH",
            "INFEK",
            "CCHHC",
            "LDLAS",
        ):
            self.assertIn(expected_seq, by_seq, f"Missing audited motif {expected_seq}")

        for m in seqle:
            chem_engine.validate_seqle_motif_record(m, aa, enforce_reference_registry=True)

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

        # Negative test 4: artificial NPLYR parent override or fabricated NPLYK flank must fail reference validation
        bad_nplyr_override = copy.deepcopy(by_seq["NPLYK"])
        bad_nplyr_override["parent_sequence"] = "AKWDTANNPLYREATST"
        with self.assertRaises(ValueError):
            chem_engine.validate_seqle_motif_record(bad_nplyr_override, aa, enforce_reference_registry=True)

        bad_nplyk_fabricated = copy.deepcopy(by_seq["NPLYK"])
        bad_nplyk_fabricated["parent_sequence"] = "ZZZZZZZNPLYKZZZZZZ"
        with self.assertRaises(ValueError):
            chem_engine.validate_seqle_motif_record(bad_nplyk_fabricated, aa, enforce_reference_registry=True)

    def test_08_p1_3_order_independent_canonical_hash(self) -> None:
        """Verify dataset_version_hash and schedule_version_hash order-invariance and decoupled sensitivity."""
        chem_copy = copy.deepcopy(self.data["chemle"])
        seq_copy = copy.deepcopy(self.data["seqle"])
        rng = random.Random(42)
        rng.shuffle(chem_copy)
        rng.shuffle(seq_copy)

        h_orig = self.data["metadata"]["dataset_version_hash"]
        s_orig = self.data["metadata"]["schedule_version_hash"]
        self.assertEqual(
            h_orig,
            build_puzzles.compute_canonical_dataset_hash(chem_copy, seq_copy, self.feedback_spec),
            "dataset_version_hash must be invariant to catalog input ordering",
        )
        self.assertEqual(
            s_orig,
            build_puzzles.compute_schedule_version_hash(chem_copy, seq_copy),
            "schedule_version_hash must be invariant to catalog input ordering",
        )

        # Changing descriptors/clues/spec updates dataset_version_hash WITHOUT shifting daily schedule_version_hash
        chem_mutated = copy.deepcopy(self.data["chemle"])
        ranked_rec = next(r for r in chem_mutated if "ranked" in r.get("game_eligibility", []))
        ranked_rec["MW"] = round(ranked_rec["MW"] + 1.0, 1)
        ranked_rec["clue_ladder"]["ranked"][0]["text"] = "Updated band text"
        self.assertNotEqual(
            h_orig,
            build_puzzles.compute_canonical_dataset_hash(chem_mutated, self.data["seqle"], self.feedback_spec),
            "dataset_version_hash must change when a gameplay descriptor changes",
        )
        self.assertEqual(
            s_orig,
            build_puzzles.compute_schedule_version_hash(chem_mutated, self.data["seqle"]),
            "schedule_version_hash must remain stable across descriptor/feedback tuning changes",
        )

        # Changing pool membership / game_eligibility updates schedule_version_hash
        seq_elig_mut = copy.deepcopy(self.data["seqle"])
        seq_elig_mut[0]["game_eligibility"] = ["learn"]
        self.assertNotEqual(
            s_orig,
            build_puzzles.compute_schedule_version_hash(self.data["chemle"], seq_elig_mut),
            "schedule_version_hash must change when game_eligibility changes",
        )

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
        self.assertIn("12168 Chemle vectors + 1800 Seqle vectors", proc.stdout)


    def test_12_p0_5_provenance_registry_and_offline_fixtures(self):
        """Verify 1:1 registry completeness, typed citation roles, COPI/N-glc validators, and offline negative fixtures."""
        import copy
        from chem_engine import (
            SEQLE_REFERENCE_REGISTRY,
            validate_copi_terminal_motif,
            validate_nglycosylation_sequon,
            validate_seqle_motif_record,
            validate_seqle_registry_completeness,
        )

        seqle_records = self.data["seqle"]
        validate_seqle_registry_completeness(seqle_records)

        # Retired KKMPQ must not be in live pool; verified KEKSD, INFEK (260-264), DSGLD (9T9W), PPPVP (1GBR context-only) must exist
        by_seq = {r["seq"]: r for r in seqle_records}
        self.assertNotIn("KKMPQ", by_seq)
        self.assertIn("KEKSD", by_seq)
        self.assertEqual(by_seq["INFEK"]["window_offset"], 2)
        self.assertEqual(by_seq["INFEK"]["citations"][0]["motif_residue_range"], [260, 264])
        dsgld_cits = by_seq["DSGLD"]["citations"]
        self.assertTrue(any(c.get("id") == "9T9W" and c.get("entity_id") == "2" for c in dsgld_cits))
        self.assertFalse(any(c.get("id") == "1P22" for c in dsgld_cits))
        pppvp_cits = by_seq["PPPVP"]["citations"]
        self.assertTrue(any(c.get("id") == "1GBR" and c.get("role") == "motif_family_context" for c in pppvp_cits))

        # 1. Missing registry key / stale extra registry key
        with self.assertRaises(ValueError):
            validate_seqle_registry_completeness(seqle_records[:-1])
        bad_rec = copy.deepcopy(seqle_records[0])
        bad_rec["seq"] = "ZZZZZ"
        with self.assertRaises(ValueError):
            validate_seqle_motif_record(bad_rec, self.data["aa_properties"], enforce_reference_registry=True)

        extra_reg = dict(SEQLE_REFERENCE_REGISTRY)
        extra_reg["ZZZZZ"] = copy.deepcopy(SEQLE_REFERENCE_REGISTRY["CAACH"])
        with self.assertRaises(ValueError):
            validate_seqle_registry_completeness(seqle_records, registry=extra_reg)

        # 2. Unknown provenance status, wrong isoform, off-by-one residue coordinates, copied flanks
        aa_props = self.data["aa_properties"]
        for mutator in [
            lambda r: r.__setitem__("provenance_status", "invented_status"),
            lambda r: r["citations"][0].__setitem__("isoform", "Q99999-2"),
            lambda r: r["citations"][0].__setitem__("motif_residue_range", [999, 1003]),
            lambda r: r.__setitem__("parent_sequence", "AAAAAAAAAAAAA"),
        ]:
            mut = copy.deepcopy(by_seq["CAACH"])
            mutator(mut)
            with self.assertRaises(ValueError):
                validate_seqle_motif_record(mut, aa_props, enforce_reference_registry=True)

        # 3. 1GBR cannot masquerade as exact_structure on PPPVP, nor wrong UniProt on DSGLD
        bad_dsgld = copy.deepcopy(by_seq["DSGLD"])
        bad_dsgld["citations"][0]["id"] = "WRONG_UNIPROT"
        with self.assertRaises(ValueError):
            validate_seqle_motif_record(bad_dsgld, aa_props, enforce_reference_registry=True)

        # 4. COPI C-terminal and N-glycosylation validators
        validate_copi_terminal_motif("KEKSD", "C")
        validate_copi_terminal_motif("AKKDA", "C")
        for invalid_copi, term in [("KKMPQ", "C"), ("KEKSD", "none"), ("AALVS", "C")]:
            with self.assertRaises(ValueError):
                validate_copi_terminal_motif(invalid_copi, term)

        validate_nglycosylation_sequon("NLTGY", "N-[^P]-[ST]")
        with self.assertRaises(ValueError):
            validate_nglycosylation_sequon("NPTGY", "N-[^P]-[ST]")
        with self.assertRaises(ValueError):
            validate_nglycosylation_sequon("NLAGY", "N-[^P]-[ST]")

    def test_13_p3_solver_separation_clue_levels_and_unsolved_bucket(self):
        """Verify legal-probe upper bound, opening groups, Level 2/3 clue conditioning, and explicit unsolved bucket."""
        import copy
        from chem_engine import (
            compute_peptide_biochem,
            generate_legal_seqle_probes,
            run_seqle_greedy_solver,
        )

        import json
        from pathlib import Path

        aa_props = self.data["aa_properties"]
        ranked_spec = json.loads((Path(__file__).resolve().parent / "data" / "feedback_spec.json").read_text(encoding="utf-8"))["profiles"]["ranked"]
        seqle_records = self.data["seqle"]
        ranked_targets = [r for r in seqle_records if "ranked" in r.get("game_eligibility", [])]

        # Confirm exact Ranked solver distribution & legal-probe opening splits
        cat_stats = run_seqle_greedy_solver(
            seqle_records, aa_props, ranked_spec, profile_id="ranked", condition_on_clues=True
        )
        self.assertEqual(cat_stats["target_count"], 28)
        self.assertEqual(cat_stats["unsolved_count"], 0)
        self.assertIn("whole_pool_opening_metrics", cat_stats)
        self.assertEqual(len(cat_stats["opening_group_metrics"]), 3)

        probes = sorted(set(r["seq"] for r in seqle_records) | set(generate_legal_seqle_probes(aa_props, seed=20261114, count=10000)))
        legal_greedy = run_seqle_greedy_solver(
            seqle_records,
            aa_props,
            ranked_spec,
            profile_id="ranked",
            condition_on_clues=True,
            answer_records=ranked_targets,
            probe_seqs=probes,
            strategy="greedy_entropy",
            include_unconditioned=False,
        )
        self.assertEqual(legal_greedy["pct_solved_le2"], 100.0)
        self.assertTrue(all(g["singleton_count"] == g["group_size"] for g in legal_greedy["opening_group_metrics"]))

        # Synthetic Level 2 & Level 3 clue unlocking fixture (all 8 targets give identical tile feedback CCCCA to wrong probes)
        def mk_target(seq, fam, ctx, rule):
            rec = copy.deepcopy(ranked_targets[0])
            rec["seq"] = seq
            rec["game_eligibility"] = ["ranked", "learn"]
            rec.update(compute_peptide_biochem(seq, aa_props))
            rec["clue_ladder"] = {
                "ranked": [
                    {"level": 0, "after_guesses": 0, "text": "Charge Band (pH 7.4): Near-Neutral"},
                    {"level": 1, "after_guesses": 2, "text": fam},
                    {"level": 2, "after_guesses": 4, "text": ctx},
                    {"level": 3, "after_guesses": 5, "text": rule},
                ]
            }
            return rec

        synth_spec = {
            **ranked_spec,
            "max_guesses": 6,
            "show_biochem_deltas": False,
            "show_positional_class_match": False,
        }
        synth = [
            mk_target("AAAAC", "FamA", "Ctx1", "Rule1"),
            mk_target("AAAAD", "FamA", "Ctx1", "Rule2"),
            mk_target("AAAAE", "FamA", "Ctx2", "Rule1"),
            mk_target("AAAAF", "FamA", "Ctx2", "Rule2"),
            mk_target("AAAAG", "FamB", "Ctx1", "Rule1"),
            mk_target("AAAAH", "FamB", "Ctx1", "Rule2"),
            mk_target("AAAAI", "FamB", "Ctx2", "Rule1"),
            mk_target("AAAAK", "FamB", "Ctx2", "Rule2"),
        ]
        # With max_steps=1, 7 targets fail to solve -> explicit unsolved bucket without omitting them from target_count!
        capped = run_seqle_greedy_solver(
            synth,
            aa_props,
            synth_spec,
            profile_id="ranked",
            condition_on_clues=True,
            answer_records=synth,
            probe_seqs=[r["seq"] for r in synth],
            max_steps=1,
            include_unconditioned=False,
        )
        self.assertEqual(capped["target_count"], 8)
        self.assertEqual(capped["solved_count"], 1)
        self.assertEqual(capped["unsolved_count"], 7)
        self.assertEqual(capped["solved_count"] + capped["unsolved_count"], capped["target_count"])

        # Confirm Level 2 (guess 4) & Level 3 (guess 5) clues solve all 8 within 6 guesses while unconditioned leaves unsolved=2
        cond_res = run_seqle_greedy_solver(
            synth,
            aa_props,
            synth_spec,
            profile_id="ranked",
            condition_on_clues=True,
            answer_records=synth,
            probe_seqs=[r["seq"] for r in synth],
            max_steps=6,
            include_unconditioned=True,
        )
        self.assertEqual(cond_res["unsolved_count"], 0)
        self.assertEqual(cond_res["pct_solved_by_cap"], 100.0)
        self.assertLess(cond_res["unconditioned_summary"]["pct_solved_by_cap"], 100.0)
        self.assertGreater(cond_res["unconditioned_summary"]["max"], 6)

    def test_14_p4_per_game_schedule_hash_isolation(self):
        """Verify Chemle & Seqle daily schedules are isolated across provenance/Learn/opposite-game edits."""
        import copy
        from build_puzzles import (
            compute_chemle_schedule_hash,
            compute_schedule_version_hash,
            compute_seqle_schedule_hash,
        )

        chemle = copy.deepcopy(self.data["chemle"])
        seqle = copy.deepcopy(self.data["seqle"])
        self.assertEqual(compute_chemle_schedule_hash(chemle), "9257bd168fa6")
        base_seqle_h = compute_seqle_schedule_hash(seqle)

        # Adding/changing Learn-only text or provenance on Chemle must NOT rotate Chemle Ranked schedule
        chemle_learn_edit = copy.deepcopy(chemle)
        chemle_learn_edit[0]["hint_summary"] = "Updated educational summary"
        chemle_learn_edit[0]["clue_ladder"]["learn"][0]["text"] = "Edited learn text"
        self.assertEqual(compute_chemle_schedule_hash(chemle_learn_edit), "9257bd168fa6")

        # Rotating a Ranked Seqle target rotates seqle_schedule_hash & combined schedule_version_hash, NOT chemle
        seqle_edit = copy.deepcopy(seqle)
        for r in seqle_edit:
            if "ranked" in r.get("game_eligibility", []):
                r["seq"] = "KKKKK"
                break
        self.assertEqual(compute_chemle_schedule_hash(chemle), "9257bd168fa6")
        self.assertNotEqual(compute_seqle_schedule_hash(seqle_edit), base_seqle_h)
        self.assertNotEqual(
            compute_schedule_version_hash(chemle, seqle_edit),
            compute_schedule_version_hash(chemle, seqle),
        )

    def test_15_chemle_pool_isolation_with_learn_only_fixture(self):
        """Learn-only compounds (even descriptor extremes) must not mutate Ranked quantiles, Level-0 clues, k-anonymity, or solver/schedule."""
        import copy
        from build_puzzles import (
            CHEMLE_RAW_COMPOUNDS,
            build_all,
            compute_chemle_schedule_hash,
            compute_seqle_schedule_hash,
        )

        from rdkit.Chem import rdMolDescriptors
        from chem_engine import standardize_and_audit

        ext_smi = "CC(=CCC/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/CCC=C(C)C)\\C)\\C)\\C)/C)/C)/C)C"
        ext_mol, ext_std = standardize_and_audit(ext_smi)
        learn_extreme = {
            "id": "synthetic_lycopene_probe",
            "name": "Synthetic Polyene Extreme",
            "aliases": [],
            "Category": "Plant Natural Product",
            "primary_class": "botanical_secondary_metabolite",
            "role_tags": ["pigment", "polyene_fixture"],
            "input_isomeric_smiles": ext_smi,
            "expected_standardized_smiles": ext_std["canonical_smiles"],
            "expected_inchikey": ext_std["inchikey"],
            "expected_formula": rdMolDescriptors.CalcMolFormula(ext_mol),
            "stereo_status": "defined_stereo",
            "protonation_policy": "neutral_canonical",
            "structure_refs": {"pubchem_cid": 446925, "chebi_id": "CHEBI:15948", "pubchem_matched_inchikey": True},
            "evidence_refs": [{"source": "PubChem", "id": "CID:446925", "url": "https://pubchem.ncbi.nlm.nih.gov/compound/446925", "role": "canonical_structure"}],
            "review_status": "verified",
            "review_notes": "Synthetic extreme-LogP test fixture.",
            "hint_role": "Lipophilic conjugated hydrocarbon pigment used to test pool separation.",
            "hint_origin": "Test fixture with elevated molecular weight and zero heteroatoms.",
            "structural_note": "All-trans conjugated polyene chain with no heteroatoms.",
            "game_eligibility": ["learn"],
            "guess_eligibility": ["learn"],
        }
        feedback_spec = copy.deepcopy(self.data["feedback_spec"])
        augmented_raw = copy.deepcopy(CHEMLE_RAW_COMPOUNDS) + [learn_extreme]
        ext_db = build_all(feedback_spec, chemle_raw_compounds=augmented_raw)

        self.assertEqual(ext_db["metadata"]["chemle_ranked_count"], 36)
        self.assertEqual(ext_db["metadata"]["chemle_ranked_guess_count"], 36)
        self.assertEqual(ext_db["metadata"]["chemle_learn_count"], 79)
        self.assertEqual(ext_db["metadata"]["chemle_learn_guess_count"], 79)
        self.assertEqual(ext_db["metadata"]["chemle_schedule_hash"], "9257bd168fa6")
        self.assertEqual(ext_db["metadata"]["seqle_schedule_hash"], "623f02623c18")
        self.assertEqual(compute_chemle_schedule_hash(ext_db["chemle"]), "9257bd168fa6")
        self.assertEqual(compute_seqle_schedule_hash(ext_db["seqle"]), "623f02623c18")

        base_by_id = {m["id"]: m for m in self.data["chemle"]}
        ext_by_id = {m["id"]: m for m in ext_db["chemle"]}
        self.assertNotIn("ranked", ext_by_id["synthetic_lycopene_probe"]["game_eligibility"])
        self.assertNotIn("ranked", ext_by_id["synthetic_lycopene_probe"]["guess_eligibility"])

        for mid, base_m in base_by_id.items():
            ext_m = ext_by_id[mid]
            self.assertEqual(ext_m["game_eligibility"], base_m["game_eligibility"])
            self.assertEqual(ext_m["guess_eligibility"], base_m["guess_eligibility"])
            if "ranked" in base_m.get("game_eligibility", []):
                self.assertEqual(ext_m["clue_ladder"]["ranked"], base_m["clue_ladder"]["ranked"])
                self.assertEqual(
                    ext_m["opening_probe_quality"]["ranked"],
                    base_m["opening_probe_quality"]["ranked"],
                )

        self.assertEqual(
            ext_db["metadata"]["k_anonymity_audit"]["chemle_ranked"],
            self.data["metadata"]["k_anonymity_audit"]["chemle_ranked"],
        )
        self.assertEqual(
            ext_db["metadata"]["solver_gate_audit"]["ranked"],
            self.data["metadata"]["solver_gate_audit"]["ranked"],
        )

        # Invariant: every profile where a compound is an answer must also allow it as a guess
        bad_raw = copy.deepcopy(CHEMLE_RAW_COMPOUNDS)
        bad_raw[0]["game_eligibility"] = ["ranked"]
        bad_raw[0]["guess_eligibility"] = ["learn"]
        with self.assertRaises(ValueError):
            build_all(feedback_spec, chemle_raw_compounds=bad_raw)

    def test_16_chemle_manifest_validation_and_negative_cases(self):
        """Verify offline manifest validation, provenance checks, and all required negative/positive fixtures."""
        import copy
        from rdkit.Chem import rdMolDescriptors
        from build_puzzles import load_chemle_manifest
        from chem_engine import (
            standardize_and_audit,
            validate_chemle_manifest,
            validate_chemle_manifest_entry,
        )

        manifest = load_chemle_manifest()
        self.assertEqual(len(manifest), 78)
        validate_chemle_manifest(manifest)

        # Every verified/baseline_legacy manifest entry surfaces provenance fields in compiled puzzles.json
        for rec in self.data["chemle"]:
            self.assertIn(rec.get("review_status"), ("verified", "baseline_legacy"))
            self.assertTrue(rec.get("structure_refs", {}).get("pubchem_cid"))
            self.assertTrue(rec.get("structure_refs", {}).get("chebi_id"))
            self.assertTrue(rec.get("evidence_refs"))

        by_id = {m["id"]: copy.deepcopy(m) for m in manifest}

        # 1. Invalid SMILES
        bad_smi = copy.deepcopy(by_id["caffeine"])
        bad_smi["input_isomeric_smiles"] = "NOT_A_VALID_SMILES((("
        with self.assertRaises(ValueError):
            validate_chemle_manifest_entry(bad_smi)

        # 2. Disconnected multi-fragment salt input
        bad_salt = copy.deepcopy(by_id["ibuprofen"])
        bad_salt["input_isomeric_smiles"] = "CC(C)Cc1ccc(C(C)C(=O)O)cc1.[Na+]"
        with self.assertRaises(ValueError):
            validate_chemle_manifest_entry(bad_salt)

        # 3. Duplicate single-ion salt/conjugate base resolving to the same uncharged parent InChIKey
        ibu_mol, ibu_std = standardize_and_audit("CC(C)Cc1ccc(C(C)C(=O)[O-])cc1")
        dup_parent = copy.deepcopy(by_id["ibuprofen"])
        dup_parent["id"] = "ibuprofen_anion"
        dup_parent["name"] = "Ibuprofen Anion"
        dup_parent["input_isomeric_smiles"] = "CC(C)Cc1ccc(C(C)C(=O)[O-])cc1"
        dup_parent["expected_standardized_smiles"] = ibu_std["canonical_smiles"]
        dup_parent["expected_inchikey"] = ibu_std["inchikey"]
        dup_parent["expected_formula"] = rdMolDescriptors.CalcMolFormula(ibu_mol)
        dup_parent["protonation_policy"] = "zwitterion_neutral_net"
        # Wait: formally charged anion fails either net_q != 0 OR parent_ik collision!
        with self.assertRaises(ValueError):
            validate_chemle_manifest([by_id["ibuprofen"], dup_parent])

        # 4. Wrong E/Z double-bond or stereoisomer against pinned expected_inchikey (cis vs trans-resveratrol)
        wrong_ez = copy.deepcopy(by_id["resveratrol"])
        wrong_ez["input_isomeric_smiles"] = "Oc1ccc(/C=C\\c2cc(O)cc(O)c2)cc1"
        with self.assertRaises(ValueError):
            validate_chemle_manifest_entry(wrong_ez)

        # 5. Label/graph conflict (wrong formula declared for caffeine graph)
        bad_formula = copy.deepcopy(by_id["caffeine"])
        bad_formula["expected_formula"] = "C7H8N4O2"
        with self.assertRaises(ValueError):
            validate_chemle_manifest_entry(bad_formula)

        # 6. Misleading stereo_status='racemate' on achiral molecule
        bad_rac = copy.deepcopy(by_id["caffeine"])
        bad_rac["stereo_status"] = "racemate"
        bad_rac["review_notes"] = "Invalid racemate claim on achiral purine."
        with self.assertRaises(ValueError):
            validate_chemle_manifest_entry(bad_rac)

        # 7. Duplicate ID and duplicate expected_inchikey
        dup_id = copy.deepcopy(by_id["caffeine"])
        with self.assertRaises(ValueError):
            validate_chemle_manifest([by_id["caffeine"], dup_id])

        dup_ik = copy.deepcopy(by_id["caffeine"])
        dup_ik["id"] = "caffeine_clone"
        dup_ik["name"] = "1,3,7-Trimethylxanthine"
        with self.assertRaises(ValueError):
            validate_chemle_manifest([by_id["caffeine"], dup_ik])

        # 8. Unreviewed playable record
        unreviewed = copy.deepcopy(by_id["caffeine"])
        unreviewed["review_status"] = "pending"
        with self.assertRaises(ValueError):
            validate_chemle_manifest_entry(unreviewed)

        # 9. Positive case: legitimate enantiomer pair (R- vs S-2-butanol) + regioisomer pair (o- vs p-xylene) pass!
        iso_specs = [
            ("r_2butanol", "(R)-2-Butanol", "CC[C@@H](C)O", "BTANRVKWQNVYAZ-SCSAIBSYSA-N", "defined_stereo", 6568, "CHEBI:41236"),
            ("s_2butanol", "(S)-2-Butanol", "CC[C@H](C)O", "BTANRVKWQNVYAZ-BYPYZUCNSA-N", "defined_stereo", 6568, "CHEBI:41237"),
            ("o_xylene", "o-Xylene", "Cc1ccccc1C", "CTQNGGLPUBDAKN-UHFFFAOYSA-N", "achiral", 7237, "CHEBI:28413"),
            ("p_xylene", "p-Xylene", "Cc1ccc(C)cc1", "URLKBWYHVLBVBO-UHFFFAOYSA-N", "achiral", 7809, "CHEBI:28633"),
        ]
        iso_entries = []
        for mid, mname, smi, _, st, pcid, chid in iso_specs:
            mol, std = standardize_and_audit(smi)
            iso_entries.append({
                "id": mid,
                "name": mname,
                "aliases": [],
                "Category": "Plant Natural Product",
                "primary_class": "isomer_pair",
                "role_tags": ["isomer_test"],
                "input_isomeric_smiles": smi,
                "expected_standardized_smiles": std["canonical_smiles"],
                "expected_inchikey": std["inchikey"],
                "expected_formula": rdMolDescriptors.CalcMolFormula(mol),
                "stereo_status": st,
                "protonation_policy": "neutral_canonical",
                "game_eligibility": ["learn"],
                "guess_eligibility": ["learn"],
                "structure_refs": {"pubchem_cid": pcid, "chebi_id": chid, "pubchem_matched_inchikey": True},
                "evidence_refs": [{"source": "PubChem", "id": f"CID:{pcid}", "url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{pcid}", "role": "canonical_structure"}],
                "hint_role": "Isomer test fixture verifying legitimate isomer pairs coexist cleanly.",
                "hint_origin": "Synthetic pair test.",
                "structural_note": "Distinct InChIKeys allow both isomers in the manifest.",
                "review_status": "verified",
                "review_notes": "Positive isomer pair check.",
            })
        validate_chemle_manifest(iso_entries)


if __name__ == "__main__":
    unittest.main(verbosity=2)
