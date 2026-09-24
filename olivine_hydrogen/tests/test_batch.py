import unittest
from dataclasses import replace
import numpy as np
from olivine_hydrogen.batch import *


class BatchTests(unittest.TestCase):
    def test_primary_solubility_benchmarks(self):
        points = [(323.19,121.706,0,.001544), (323.21,100.354,1,.000972),
                  (323.20,100.832,3,.000631), (323.19,66.385,5,.000293)]
        for t,p,m,x in points:
            actual = solubility_mol_kg(t,p,m)
            expected = x/(M_WATER*(1-x))
            self.assertLess(abs(actual/expected-1), .04)

    def test_domain_guard(self):
        for t,p,m in [(323.15,1,4.5),(300,120,4.5),(323.15,120,5.1)]:
            with self.assertRaises(ValueError): solubility_mol_kg(t,p,m)
        for c in [replace(Setup(),warm_system_extra_dp_bar=-1),
                  replace(Setup(),feed_pressure_bar=0),replace(Setup(),feed_temperature_k=0)]:
            with self.assertRaises(ValueError): c.validate()

    def test_balanced_redox(self):
        # 3FeO+H2O -> Fe3O4+H2: Fe,H,O atom counts.
        np.testing.assert_equal(np.array([3,0,3])+[0,2,1], np.array([3,0,4])+[0,2,0])

    def test_half_time_and_exhaustion(self):
        c=Setup()
        n=reaction_extent([0,365,1e9],365,c)
        np.testing.assert_allclose(n,capacity_mol(c)*np.array([0,.5,1]))
        self.assertAlmostEqual(capacity_mol(c)*M_H2,7.1218,delta=.002)

    def test_all_ledgers(self):
        c=Setup(); x=evaluate([0,3,365,3650,1e8],365,c)
        np.testing.assert_allclose(x['generated_h2_kg'], x['free_reactor_h2_kg']+x['dissolved_reactor_h2_kg'])
        np.testing.assert_allclose(x['h2_ledger_error_kg'],0,atol=1e-13)
        np.testing.assert_allclose(x['water_remaining_kg']+x['generated_h2_kg']/M_H2*M_WATER,c.water_kg)
        np.testing.assert_allclose(x['salt_kg'],c.water_kg*c.salt_molality*M_SALT)
        self.assertTrue(np.all(x['remaining_reactive_fe_kg']>=0))

    def test_zero_iron_zero_accessibility_zero_rate(self):
        for c,h in [(replace(Setup(),fe_mass_fraction=0),3),
                    (replace(Setup(),accessible_fraction=0),3),(Setup(),np.inf)]:
            x=evaluate([0,3,1000],h,c)
            np.testing.assert_array_equal(x['generated_h2_kg'],0)
            np.testing.assert_array_equal(x['gross_electricity_kwh'],0)

    def test_no_collection_no_electricity(self):
        x=evaluate([3,10000],3,replace(Setup(),collection_efficiency=0))
        np.testing.assert_array_equal(x['gross_electricity_kwh'],0)

    def test_no_initial_hydrogen_and_repeatability(self):
        x=evaluate([0,3],365)
        self.assertEqual(x['generated_h2_kg'][0],0)
        self.assertEqual(evaluate([3],365)['generated_h2_kg'][0],x['generated_h2_kg'][1])

    def test_pumping_units_and_heat(self):
        c=Setup();x=evaluate([3],365,c)
        v=(c.water_kg+c.water_kg*c.salt_molality*M_SALT)/c.brine_density_kg_m3
        self.assertAlmostEqual(x['initial_pump_electricity_kwh'][0],119*1e5*v/(.7*3.6e6))
        self.assertLess(x['standalone_balance_kwh'][0],0)

    def test_water_bound_even_if_redox_limit_is_large(self):
        c=replace(Setup(), water_kg=.01, salt_molality=0)
        self.assertLessEqual(reaction_extent([1e8],3,c)[0],.01/M_WATER)


if __name__=='__main__': unittest.main()
