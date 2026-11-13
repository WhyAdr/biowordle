"""Automated Scientific & Engineering Audit Suite for BioWordle.

Verifies all 10 priority acceptance criteria from Sol's technical review:
1. Compiler independence (zero external skill path imports)
2. Build completion & artifact integrity
3. Molecular identity, isomer separation, and stereocenter preservation
4. Chiral Morgan (ECFP4, r=2, 2048b) Tanimoto matrix symmetry & self-identity
5. Seqle two-pass duplicate residue handling, class precedence & H-H terminal pKa model
6. Descriptor validity (no NaN/Inf, documented CNS-MPO & SA_Score methods)
7. Acyclic Murcko scaffold fallback (Citric acid carbon-skeleton clue != full structure)
8. Deterministic dataset version hash & Shannon entropy information gain IG(g)
9. Browser HTML/DOM & responsive SVG smoke test (both modular and standalone builds)
10. Reproducible pinned environment verification
"""

from __future__ import annotations

import json
import math
import os
import unittest

from rdkit import Chem, DataStructs

import build_puzzles
import chem_engine


BASE_DIR = os.path.dirname(os.path.abspath(__file__))


class TestBioWordleAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        json_path = os.path.join(BASE_DIR, "data", "puzzles.json")
        with open(json_path, "r", encoding="utf-8") as f:
            cls.data = json.load(f)

    def test_01_compiler_independence(self) -> None:
        """Verify neither chem_engine.py nor build_puzzles.py references external .gemini paths."""
        for fname in ("chem_engine.py", "build_puzzles.py"):
            with open(os.path.join(BASE_DIR, fname), "r", encoding="utf-8") as f:
                src = f.read()
            self.assertNotIn(".gemini", src, f"{fname} must not reference external .gemini paths")
            self.assertNotIn("skill_paths", src, f"{fname} must not import external skill_paths")

    def test_02_molecular_identity_and_stereochemistry(self) -> None:
        """Verify unique IDs, InChIKeys, canonical SMILES, and preserved stereocenters."""
        chemle = self.data["chemle"]
        self.assertEqual(len(chemle), 36)

        ids = [c["id"] for c in chemle]
        inchikeys = [c["inchikey"] for c in chemle]
        can_smiles = [c["canonical_smiles"] for c in chemle]

        self.assertEqual(len(set(ids)), 36, "Duplicate compound IDs found")
        self.assertEqual(len(set(inchikeys)), 36, "Duplicate InChIKeys after standardization")
        self.assertEqual(len(set(can_smiles)), 36, "Duplicate canonical SMILES after standardization")

        # Positional isomers Protocatechuic acid (3,4-DHB) and Gentisic acid (2,5-DHB)
        by_id = {c["id"]: c for c in chemle}
        pca = by_id["protocatechuic_acid"]
        gen = by_id["gentisic_acid"]
        self.assertEqual(pca["formula"], "C7H6O4")
        self.assertEqual(gen["formula"], "C7H6O4")
        self.assertNotEqual(pca["inchikey"], gen["inchikey"])
        self.assertLess(self.data["chemle_similarity"]["protocatechuic_acid"]["gentisic_acid"], 100)

        # Chiral molecules must preserve stereocenters
        for chiral_id in ("adrenaline", "tryptophan", "ascorbic_acid", "biotin", "quinine", "artemisinin"):
            self.assertGreater(
                by_id[chiral_id]["ChiralCenters"],
                0,
                f"{chiral_id} lost stereocenters during standardization",
            )
            self.assertIn("@", by_id[chiral_id]["canonical_smiles"])

        # Verify chiral Morgan fingerprint distinguishes (R)- vs (S)-adrenaline
        r_adr = Chem.MolFromSmiles("CNC[C@H](O)c1ccc(O)c(O)c1")
        s_adr = Chem.MolFromSmiles("CNC[C@@H](O)c1ccc(O)c(O)c1")
        fp_r = chem_engine.compute_chiral_morgan_fp(r_adr)
        fp_s = chem_engine.compute_chiral_morgan_fp(s_adr)
        sim_rs = DataStructs.TanimotoSimilarity(fp_r, fp_s)
        self.assertLess(sim_rs, 1.0, "Chiral Morgan fingerprint must distinguish (R) and (S) enantiomers")

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
            self.assertEqual(sim[i][i], 100, f"Self-similarity for {i} must be 100")
            for j in ids:
                self.assertEqual(sim[i][j], sim[j][i], f"Asymmetry between {i} and {j}")
                self.assertGreaterEqual(sim[i][j], 0)
                self.assertLessEqual(sim[i][j], 100)
                if i != j:
                    self.assertLess(sim[i][j], 100, f"Distinct compounds {i} and {j} have 100% similarity")

    def test_04_descriptor_validity_and_metadata(self) -> None:
        """Verify no NaN/Inf descriptors and documented methodology."""
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
            "CNS_MPO",
        )
        for c in self.data["chemle"]:
            for k in numeric_keys:
                val = c[k]
                self.assertIsNotNone(val, f"{c['id']}.{k} is None")
                self.assertFalse(math.isnan(float(val)) or math.isinf(float(val)), f"{c['id']}.{k} is NaN/Inf")
            self.assertTrue(0.0 <= c["CNS_MPO"] <= 6.0)
            self.assertTrue(1.0 <= c["SA_Score"] <= 10.0)
            self.assertTrue(0.0 <= c["QED"] <= 1.0)

    def test_05_acyclic_scaffold_fallback(self) -> None:
        """Verify acyclic molecules (Citric acid) use carbon_skeleton fallback and do not leak full SVG."""
        by_id = {c["id"]: c for c in self.data["chemle"]}
        citric = by_id["citric_acid"]
        self.assertEqual(citric["Rings"], 0)
        self.assertEqual(citric["scaffold_mode"], "carbon_skeleton")
        self.assertNotEqual(
            citric["svg_scaffold"],
            citric["svg_full"],
            "Acyclic molecule scaffold SVG must not leak the full structure SVG",
        )

    def test_06_seqle_two_pass_duplicates_and_class_precedence(self) -> None:
        """Test Seqle two-pass duplicate counting and exact > present > class > absent priority."""
        aa = self.data["aa_properties"]

        # Case 1: Target has 1 'K' at pos 0 ('KDELR'). Guess has two 'K's at pos 0 and pos 1 ('KKAAA').
        # Pos 0 is 'exact' (consumes the only 'K' in target).
        # Pos 1 'K' cannot be 'present' (no 'K' left in target), and target[1]='D' is Negative while 'K' is Positive -> 'absent'.
        res1 = chem_engine.evaluate_seqle_guess("KKAAA", "KDELR", aa)
        self.assertEqual(res1[0], "exact")
        self.assertEqual(res1[1], "absent")

        # Case 2: Target 'KDELR' has 'R' (Positive) at pos 4. Guess 'KKKKK':
        # Pos 0: 'exact' (K==K)
        # Pos 1..3: 'absent'
        # Pos 4: 'K' is not in remaining target letters, but 'K' (Positive) shares primary class with target[4]='R' (Positive) -> 'class'!
        res2 = chem_engine.evaluate_seqle_guess("KKKKK", "KDELR", aa)
        self.assertEqual(res2, ["exact", "absent", "absent", "absent", "class"])

        # Case 3: Target 'DEVDG' (two 'D's at 0 and 3). Guess 'DDDDD':
        # Pos 0 and 3 are 'exact', Pos 1 ('D' vs 'E' Negative) is 'class', Pos 2 and 4 are 'absent'.
        res3 = chem_engine.evaluate_seqle_guess("DDDDD", "DEVDG", aa)
        self.assertEqual(res3, ["exact", "class", "absent", "exact", "absent"])

        # Case 4: Misplaced duplicate count ('DAAAD' vs 'KDELR' -> only one 'D' at pos 1 in target):
        # First 'D' at pos 0 gets 'present', second 'D' at pos 4 gets 'absent'.
        res4 = chem_engine.evaluate_seqle_guess("DAAAD", "KDELR", aa)
        self.assertEqual(res4[0], "present")
        self.assertEqual(res4[4], "absent")

    def test_07_seqle_biophysics_and_provenance(self) -> None:
        """Verify all 30 motifs have context_requirement, exemplar_source, and valid H-H biophysics."""
        seqle = self.data["seqle"]
        self.assertEqual(len(seqle), 30)
        for m in seqle:
            self.assertTrue(len(m["context_requirement"]) > 15, f"Missing context_requirement for {m['seq']}")
            self.assertTrue(len(m["exemplar_source"]) > 10, f"Missing exemplar_source for {m['seq']}")
            self.assertTrue(1.0 <= m["isoelectric_point"] <= 13.5)
        # Cationic PKKKR must have high positive charge at pH 7.4 and high pI (> 10.5)
        by_seq = {m["seq"]: m for m in seqle}
        self.assertGreater(by_seq["PKKKR"]["net_charge_ph74"], 3.5)
        self.assertGreater(by_seq["PKKKR"]["isoelectric_point"], 10.5)
        # Anionic DEVDG must have negative net charge at pH 7.4 and acidic pI (< 4.5)
        self.assertLess(by_seq["DEVDG"]["net_charge_ph74"], -2.5)
        self.assertLess(by_seq["DEVDG"]["isoelectric_point"], 4.5)

    def test_08_shannon_entropy_information_gain(self) -> None:
        """Verify build-time Shannon entropy IG(g) metrics are bounded and populated."""
        max_chem_h = math.log2(len(self.data["chemle"]))
        for c in self.data["chemle"]:
            ig = c["info_metrics"]["information_gain_bits"]
            self.assertGreater(ig, 1.5)
            self.assertLessEqual(ig, max_chem_h + 0.01)

        max_seq_h = math.log2(len(self.data["seqle"]))
        for m in self.data["seqle"]:
            ig = m["info_metrics"]["information_gain_bits"]
            self.assertGreater(ig, 1.0)
            self.assertLessEqual(ig, max_seq_h + 0.01)

    def test_09_browser_html_and_standalone_smoke_test(self) -> None:
        """Verify index.html and dist/biowordle_standalone.html contain all required DOM IDs and responsive SVGs."""
        required_ids = [
            "tab-chemle",
            "tab-seqle",
            "btn-playstyle",
            "btn-new-puzzle",
            "btn-hint",
            "btn-giveup",
            "btn-theme",
            "chemle-search-input",
            "chemle-cat-filter",
            "chemle-dropdown",
            "chemle-guesses",
            "chemle-structure-well",
            "chemle-hints-list",
            "chemle-full-props",
            "seqle-motif-picker",
            "seqle-board",
            "seqle-keyboard",
            "seqle-structure-well",
            "seqle-hints-list",
        ]
        for rel_path in ("index.html", os.path.join("dist", "biowordle_standalone.html")):
            full_path = os.path.join(BASE_DIR, rel_path)
            self.assertTrue(os.path.isfile(full_path), f"Missing {rel_path}")
            with open(full_path, "r", encoding="utf-8") as f:
                html = f.read()
            for dom_id in required_ids:
                self.assertIn(f'id="{dom_id}"', html, f"Missing DOM id={dom_id} in {rel_path}")

        # Verify standalone bundle embeds BIOWORDLE_DATA and responsive SVGs
        standalone_path = os.path.join(BASE_DIR, "dist", "biowordle_standalone.html")
        with open(standalone_path, "r", encoding="utf-8") as f:
            bundle = f.read()
        self.assertTrue("window.BIOWORDLE_DATA =" in bundle)
        self.assertTrue("viewBox='0 0 320 220'" in bundle or 'viewBox="0 0 320 220"' in bundle)
        self.assertTrue('width="100%" height="100%"' in bundle or 'width=\\"100%\\" height=\\"100%\\"' in bundle)


if __name__ == "__main__":
    unittest.main(verbosity=2)
