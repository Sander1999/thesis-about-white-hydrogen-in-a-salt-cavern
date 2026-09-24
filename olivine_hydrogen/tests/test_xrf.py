"""Data/chemistry regressions against independently stated input expectations."""
import copy
import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np

from olivine_hydrogen.data.process_xrf import (
    RAW_FILE, FE_KG_MOL, H2_KG_MOL, FORSTERITE_SIO2_MASS_FRACTION,
    load_readings, element_statistics, summarize, parse_concentration,
    iron_capacity, legacy_assumed_mg_number, process,
)


class XrfDataTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.readings, cls.channels = load_readings()

    def test_dataset_counts_and_groups(self):
        self.assertEqual(len(self.readings), 56)
        summary, elements = summarize(self.readings, self.channels)
        self.assertEqual(len(summary), 11)
        self.assertEqual(len(elements), 11 * 33)
        self.assertEqual(sum(s["n_readings"] for s in summary), 56)
        self.assertEqual([s["n_readings"] for s in summary], [5, 5, 6, 5, 5, 5, 5, 5, 5, 5, 5])
        self.assertEqual(summary[3]["reading_numbers"], "4474 4475 4476 4477 4478")
        self.assertEqual(summary[3]["subgroup"], "big grains")

    def test_means_do_not_mix_concentrations_and_sigma(self):
        sample1 = [r for r in self.readings if r["sample"] == "sample1"]
        stats = element_statistics(sample1, "Fe")
        self.assertAlmostEqual(stats["mean_ppm"], 58452.34844, places=7)
        self.assertAlmostEqual(stats["sample_sd_ppm"], 1274.96393307, places=7)
        self.assertNotAlmostEqual(stats["mean_ppm"], 29311.94285, places=3)
        self.assertAlmostEqual(stats["sem_ppm"], 1274.96393307 / np.sqrt(5), places=7)
        instrument = np.array([173.4674, 171.4910, 168.4424, 175.0654, 169.2201])
        self.assertAlmostEqual(stats["nominal_instrument_1sigma_mean_ppm"], np.linalg.norm(instrument / 2) / 5, places=10)

    def test_uncertainty_cannot_create_an_element_detection(self):
        stats = element_statistics(self.readings, "Bi")
        self.assertEqual(stats["n_detected"], 0)
        self.assertEqual(stats["n_below_detection"], 56)
        self.assertIsNone(stats["mean_ppm"])
        self.assertIsNone(stats["nominal_instrument_1sigma_mean_ppm"])
        self.assertTrue(any(r["reported_2sigma_ppm"]["Bi"] > 0 for r in self.readings))

    def test_censored_values_keep_status_and_optional_limit(self):
        self.assertEqual(parse_concentration("<LOD"), (None, "below_detection", None))
        self.assertEqual(parse_concentration("< 3.5"), (None, "below_detection", 3.5))
        self.assertEqual(parse_concentration(""), (None, "missing", None))
        self.assertEqual(parse_concentration("0"), (0.0, "detected", None))
        for invalid in ["NaN", "inf", "-1", "<-3"]:
            with self.assertRaises(ValueError):
                parse_concentration(invalid)

    def test_theoretical_capacity_uses_elemental_fe_moles_and_hydrogen_molar_mass(self):
        # 55,845 ppm = 55.845 g Fe per kg rock = exactly one mole Fe/kg.
        capacity = iron_capacity(55845)
        self.assertAlmostEqual(float(capacity["magnetite_H2_mol_per_kg_rock"]), 1 / 3, places=14)
        self.assertAlmostEqual(float(capacity["ferric_ceiling_H2_mol_per_kg_rock"]), 1 / 2, places=14)
        self.assertAlmostEqual(float(capacity["magnetite_H2_kg_per_tonne_rock"]), 0.67196, places=12)
        self.assertAlmostEqual(float(capacity["ferric_ceiling_H2_kg_per_tonne_rock"]), 1.00794, places=12)
        scaled = iron_capacity(55845, ferrous_fraction=0.6, accessible_fraction=0.2)
        self.assertAlmostEqual(float(scaled["magnetite_H2_mol_per_kg_rock"]), 0.04, places=14)
        self.assertEqual(float(iron_capacity(55845, accessible_fraction=0)["magnetite_H2_mol_per_kg_rock"]), 0)
        with self.assertRaises(ValueError):
            iron_capacity(55845, ferrous_fraction=1.1)

    def test_legacy_mg_is_only_a_declared_forsterite_assumption(self):
        self.assertAlmostEqual(FORSTERITE_SIO2_MASS_FRACTION, 60.0843 / 140.6931, places=12)
        self.assertGreater(FORSTERITE_SIO2_MASS_FRACTION, 0.427)
        self.assertNotIn("Mg", self.channels)
        # Under the counterfactual assumption, one kg pure forsterite has
        # 2 / .1406931 mol Mg; choose the same Fe moles to obtain Mg#=.5.
        bal_ppm = 1e5
        fe_ppm = 1e6 * (0.1 * 2 / 0.1406931) * FE_KG_MOL
        self.assertAlmostEqual(float(legacy_assumed_mg_number(fe_ppm, bal_ppm)), 0.5, places=12)

    def test_nominal_scatter_propagates_linearly_to_capacity_not_to_a_rate(self):
        summaries, _ = summarize(self.readings, self.channels)
        row = summaries[0]
        expected = row["Fe_sample_sd_ppm"] * 1e-6 / FE_KG_MOL / 3 * H2_KG_MOL * 1000
        self.assertAlmostEqual(row["magnetite_H2_sample_sd_kg_per_tonne_rock"], expected, places=14)
        self.assertFalse(any("per_day" in key or "rate" in key for key in row))

    def test_reference_units_counts_and_duplicate_ids_fail_explicitly(self):
        with RAW_FILE.open() as stream:
            original = list(csv.DictReader(stream))
        def write(tmp, rows):
            path = Path(tmp) / "fixture.csv"
            with path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader(); writer.writerows(rows)
            return path
        with tempfile.TemporaryDirectory(prefix="olivine_xrf_test_", dir="/private/tmp") as tmp:
            wrong_units = copy.deepcopy(original)
            wrong_units[0]["Units"] = "wt%"
            with self.assertRaisesRegex(ValueError, "ppm"):
                load_readings(write(tmp, wrong_units))
            with self.assertRaisesRegex(ValueError, "56"):
                load_readings(write(tmp, original[:-1]))
            duplicate = copy.deepcopy(original)
            duplicate[-1]["Reading No"] = duplicate[0]["Reading No"]
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                load_readings(write(tmp, duplicate))

    def test_read_order_does_not_change_original_subgroup_mapping(self):
        with RAW_FILE.open() as stream:
            original = list(csv.DictReader(stream))
        with tempfile.TemporaryDirectory(prefix="olivine_xrf_test_", dir="/private/tmp") as tmp:
            path = Path(tmp) / "reversed.csv"
            with path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(original[0]))
                writer.writeheader(); writer.writerows(reversed(original))
            reordered, channels = load_readings(path)
            reference, _ = summarize(self.readings, self.channels)
            result, _ = summarize(reordered, channels)
            for expected, actual in zip(reference, result):
                self.assertEqual(expected["reading_numbers"], actual["reading_numbers"])
                self.assertEqual(expected["subgroup"], actual["subgroup"])
                self.assertAlmostEqual(expected["Fe_mean_ppm"], actual["Fe_mean_ppm"], places=7)

    def test_export_contains_no_imputed_measured_magnesium(self):
        with tempfile.TemporaryDirectory(prefix="olivine_xrf_test_", dir="/private/tmp") as tmp:
            audit, summary = process(output=tmp, make_plots=False)
            with (Path(tmp) / "processed/readings.csv").open() as stream:
                headers = next(csv.reader(stream))
            self.assertNotIn("Mg_ppm", headers)
            self.assertEqual(audit["rows"], 56)
            self.assertIn("legacy_assumption_Mg_number_mean", summary[0])
            self.assertTrue((Path(tmp) / "audit.json").is_file())


if __name__ == "__main__":
    unittest.main()
