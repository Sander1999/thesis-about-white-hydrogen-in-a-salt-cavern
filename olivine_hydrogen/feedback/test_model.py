"""Independent balances, null cases and numerical sensitivity for feedback ODE."""
from dataclasses import replace
import unittest
import numpy as np
from olivine_hydrogen.batch import M_H2, M_WATER, capacity_mol
from .model import Config, MODES, simulate, feedback_factors


class FeedbackTests(unittest.TestCase):
    def test_closed_and_open_component_ledgers(self):
        for mode in MODES:
            h,s=simulate(mode,replace(Config(),film_damkohler=1.,product_inhibition_coefficient=1.))
            for key in ("h2_balance_error_kg","hydrogen_recovery_balance_error_kg"):
                self.assertLess(np.max(np.abs(h[key])),1e-9)
            for key in ("water_balance_error_kg","salt_balance_error_kg"):
                self.assertLess(np.max(np.abs(h[key])),1e-6)
            self.assertTrue(np.all(h["generated_h2_kg"]>=-1e-12))
            self.assertLessEqual(h["generated_h2_kg"][-1],s["hydrogen_capacity_kg"]*(1+1e-10))
            np.testing.assert_allclose(h["solid_mass_gain_kg"]+h["generated_h2_kg"],h["water_consumed_kg"],atol=1e-12)
            if mode!="once_through":
                np.testing.assert_array_equal(h["export_h2_mol"],0.)
                np.testing.assert_array_equal(h["feed_water_kg"],0.)

    def test_analytic_intrinsic_null(self):
        c=replace(Config(),solid_expansion_fraction=0.,maximum_area_factor=1.,
                  passivation_coefficient=0.,porosity_access_exponent=0.,
                  film_damkohler=0.,product_inhibition_coefficient=0.)
        expected=capacity_mol(c)*M_H2*(1-np.exp(-np.log(2)*np.array([0.,3.,365.])/365))
        for mode in MODES:
            h,_=simulate(mode,c,[0.,3.,365.])
            np.testing.assert_allclose(h["generated_h2_kg"],expected,rtol=1e-9,atol=1e-10)

    def test_zero_source_and_no_iron(self):
        for c in (replace(Config(),rate_multiplier=0.,duration_days=3.),replace(Config(),fe_mass_fraction=0.,duration_days=3.)):
            for mode in MODES:
                h,_=simulate(mode,c,[0.,3.]);self.assertEqual(h["generated_h2_kg"][-1],0.)
                self.assertEqual(h["total_deliverable_h2_kg"][-1],0.)

    def test_geometric_and_source_bounds(self):
        for expansion,relief in ((0.,0.),(.4,.5),(.6,0.),(.6,1.)):
            c=replace(Config(),solid_expansion_fraction=expansion,crack_relief_fraction=relief)
            g=feedback_factors(np.linspace(0,1,1000),c)
            self.assertGreaterEqual(g["porosity"].min(),c.porosity_floor)
            self.assertLessEqual(g["porosity"].max(),c.initial_porosity)
            self.assertTrue(np.all((g["area_factor"]>=1)&(g["area_factor"]<=5)))
            self.assertTrue(np.all((g["permeability_ratio"]>0)&(g["permeability_ratio"]<=1)))
        h,_=simulate("stationary",replace(Config(),duration_days=1095.,half_time_days=10.,solid_expansion_fraction=.6,crack_relief_fraction=0.))
        self.assertLessEqual(h["generated_h2_kg"][-1],capacity_mol(Config())*M_H2)
        self.assertGreaterEqual(h["water_kg"][-1],0.)

    def test_hydration_and_salt_accounting(self):
        c=replace(Config(),solid_expansion_fraction=0.,maximum_area_factor=1.,passivation_coefficient=0.,porosity_access_exponent=0.)
        h,_=simulate("stationary",c,[0.,365.])
        expected=.5*c.hydration_water_kg_per_kg_rock*c.rock_kg+.5*capacity_mol(c)*M_WATER
        self.assertAlmostEqual(h["water_consumed_kg"][-1],expected,places=7)
        self.assertAlmostEqual(h["salt_molality"][-1],c.water_kg*c.salt_molality/h["water_kg"][-1],places=12)

    def test_parcel_recovery_not_pooling(self):
        h,_=simulate("once_through",Config(),[0.,3.,365.])
        self.assertGreater(h["exported_h2_kg"][-1],h["retained_h2_kg"][-1])
        # A year of dilute effluent must not be pooled to manufacture a rich gas stream.
        self.assertEqual(h["outlet_collected_h2_kg"][-1],0.)
        self.assertAlmostEqual(h["feed_water_kg"][-1],Config().water_flow_kg_day*365,places=5)

    def test_time_refinement(self):
        c=replace(Config(),film_damkohler=1.,product_inhibition_coefficient=1.)
        for mode in MODES:
            h,_=simulate(mode,c,[0.,3.,365.])
            fine,_=simulate(mode,replace(c,max_step_days=.5),[0.,3.,365.])
            for key in ("generated_h2_kg","retained_h2_kg","total_deliverable_h2_kg","flow_pump_kwh"):
                np.testing.assert_allclose(h[key],fine[key],rtol=2e-6,atol=2e-7)

    def test_domain_guard(self):
        for c in (replace(Config(),maximum_area_factor=6.),replace(Config(),crack_relief_fraction=1.1),replace(Config(),initial_porosity=.005)):
            with self.assertRaises(ValueError):simulate("stationary",c)
        with self.assertRaises(ValueError):simulate("stationary",replace(Config(),hydration_water_kg_per_kg_rock=2.,half_time_days=1.))


if __name__=="__main__":unittest.main()
