"""Properties and native conservation/refinement checks for the packed bed.

Run: ./run_darts.sh olivine_hydrogen/darts/test_model.py
Native diagnostic results are retained in output/validation.json; temporary
simulation files are removed. Fast property tests are also unittest-compatible.
"""
import contextlib,io,json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
if str(Path(__file__).resolve().parents[2]) not in sys.path:sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from olivine_hydrogen.darts.properties import load_config,BrineFlash,solubility_x,MW,FE_MW
from olivine_hydrogen.darts.model import Model
from olivine_hydrogen.darts.run import run_case
from darts.models.darts_model import DartsModel

BENCHMARKS=[(323.19,121.706,0.,.001544),(323.21,100.354,1.,.000972),(323.20,100.832,3.,.000631),(323.19,66.385,5.,.000293)]


class PropertiesTest(unittest.TestCase):
    def test_chabab_measured_reference(self):
        for t,p,m,observed in BENCHMARKS:
            self.assertLess(abs(float(solubility_x(p,t,m))/observed-1),.035)
    def test_analytical_flash_partition_conserves_salt(self):
        system=BrineFlash(load_config())
        for p in [10.,120.,230.]:
            for h in [1e-10,1e-5,.001,.05]:
                z=np.array([h,1.,4.5*MW[1]/1000]);z/=z.sum()
                q=system.state(p,z)
                np.testing.assert_allclose(q['nu']@q['x'],z,rtol=1e-12,atol=1e-14)
                self.assertEqual(q['x'][0,2],0.)
                self.assertAlmostEqual(q['molality'],4.5)
                if q['nu'][0]>0:
                    aq=q['x'][1,0]/(q['x'][1,0]+q['x'][1,1])
                    self.assertAlmostEqual(aq,solubility_x(p,323.15,4.5),places=12)
    def test_rejected_attempt_does_not_consume_iron(self):
        model=Model();model.init();before=model.fe_remaining.copy()
        with patch.object(DartsModel,'run_timestep',return_value=False):
            self.assertFalse(model.run_timestep(.1,0.))
        np.testing.assert_array_equal(model.fe_remaining,before)
        self.assertEqual(model.generated_kmol,0.)


def validate():
    result=unittest.TextTestRunner(stream=io.StringIO()).run(unittest.defaultTestLoader.loadTestsFromTestCase(PropertiesTest))
    if not result.wasSuccessful():raise AssertionError(str(result.failures)+str(result.errors))
    cases={};checks={}
    with tempfile.TemporaryDirectory(prefix='packed_bed_tests_',dir='/private/tmp') as tmp:
        for name,overrides in [('baseline',{}),('zero_source',{'reaction_enabled':False}),
                ('half_timestep',{'max_timestep_days':.01}),('spatial_5x5x5',{'nx':5,'ny':5,'nz':5}),
                ('fine_interpolation',{'obl_pressure_step_bar':.1,'obl_hydrogen_step':1e-6,'obl_water_step':1e-6}),
                ('rapid_reaction_stress',{'half_time_days':.5,'max_timestep_days':.01})]:
            with (Path(tmp)/(name+'.stdout.log')).open('w') as f,contextlib.redirect_stdout(f):
                summary,history=run_case(overrides,Path(tmp)/name,plots=False)
            cases[name]=summary
            checks[name+'_native_balance']=max(summary['max_native_relative_balance_error'].values())<2e-5
            checks[name+'_finite_iron']=0<=summary['generated_H2_kg']<=summary['redox_capacity_H2_kg']*(1+1e-12)
            initial_fe=summary['parameters']['rock_mass_kg']*summary['parameters']['fe_mass_fraction']
            fe_remaining=summary['remaining_accessible_Fe_kg'];generated=summary['generated_H2_kg']/MW[0]
            checks[name+'_Fe_ledger']=abs(fe_remaining+3*generated*FE_MW-initial_fe)<1e-8
            checks[name+'_full_fluid_solid_mass']=abs(sum(history[-1][n+'_balance_kg'] for n in ['H2','H2O','NaCl']))<1e-5
        checks['zero_source_no_generation']=cases['zero_source']['generated_H2_kg']==0
        checks['stress_exsolves_H2']=cases['rapid_reaction_stress']['gas_H2_exported_kg']>0
        checks['default_outlet_separator_no_gas']=cases['baseline']['separator_H2_available_kg']==0.
        checks['stress_outlet_separator_produces_gas']=cases['rapid_reaction_stress']['separator_H2_available_kg']>0.
        spatial=cases['spatial_5x5x5'];base=cases['baseline']
        metrics={'spatial_aqueous_H2_export_relative_difference':abs(base['aqueous_H2_exported_kg']-spatial['aqueous_H2_exported_kg'])/spatial['aqueous_H2_exported_kg'],
                 'spatial_remaining_dissolved_H2_difference_kg':spatial['remaining_dissolved_H2_kg']-base['remaining_dissolved_H2_kg'],
                 'spatial_pressure_range_difference_bar':[spatial['final_pressure_range_bar'][j]-base['final_pressure_range_bar'][j] for j in range(2)]}
        for name in ['half_timestep','fine_interpolation']:
            a=cases['baseline']['aqueous_H2_exported_kg'];b=cases[name]['aqueous_H2_exported_kg']
            metrics[name+'_aqueous_H2_export_relative_difference']=abs(a-b)/b
            checks[name+'_export_convergence']=abs(a-b)/b<.01
    report={'passed':all(checks.values()),'property_tests':result.testsRun,'checks':{k:bool(v) for k,v in checks.items()},'metrics':metrics,
        'Chabab_measured_reference':[{'temperature_k':t,'pressure_bar':p,'molality':m,'observed_x':x,'predicted_x':float(solubility_x(p,t,m))} for t,p,m,x in BENCHMARKS],
        'results':cases,'limitations':'Time and OBL refinement are tested. A5x5x5 spatial sensitivity changes corner-cell well locations and completion lengths as well as cell sizes; it does not establish spatial convergence. The assumed spatial mesh/permeability and reaction half-time are not calibrated. The rapid-reaction case is solely a numerical phase-appearance/finite-inventory stress test.'}
    output=Path(__file__).with_name('output')/'validation.json';output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'passed':report['passed'],'checks':report['checks'],'metrics':metrics},indent=2))
    if not report['passed']:raise AssertionError('Packed-bed validation failed')
    return report

if __name__=='__main__':validate()
