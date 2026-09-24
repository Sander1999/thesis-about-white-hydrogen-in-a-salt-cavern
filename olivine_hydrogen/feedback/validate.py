"""Retain actual unit-test and timestep-refinement results as validation.json."""
from dataclasses import replace
from pathlib import Path
import json
import unittest
import numpy as np
from .model import Config, MODES, simulate
from .test_model import FeedbackTests


def validate(output=None):
    tests=unittest.defaultTestLoader.loadTestsFromTestCase(FeedbackTests)
    result=unittest.TextTestRunner(verbosity=2).run(tests)
    metrics={}
    cfg=replace(Config(),film_damkohler=1.,product_inhibition_coefficient=1.)
    for mode in MODES:
        h,s=simulate(mode,cfg,[0.,3.,365.])
        fine,_=simulate(mode,replace(cfg,max_step_days=.5),[0.,3.,365.])
        metrics[mode]={"max_abs_balance_errors":s["max_abs_balance_errors"],
                       "half_step_relative_endpoint_changes":{
                           key:abs(float(fine[key][-1]-h[key][-1]))/max(abs(float(fine[key][-1])),1e-12)
                           for key in ("generated_h2_kg","retained_h2_kg","total_deliverable_h2_kg","flow_pump_kwh")}}
    report={"passed":result.wasSuccessful(),"tests_run":result.testsRun,
            "failures":len(result.failures),"errors":len(result.errors),
            "refinement_coarse_max_step_days":1.,"refinement_fine_max_step_days":.5,
            "metrics":metrics,
            "scope":"Numerical and accounting checks; not validation of reaction or fracture hypotheses."}
    output=Path(output or Path(__file__).with_name("output")/"validation.json")
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    if not result.wasSuccessful():raise RuntimeError("Feedback validation failed")
    return report


if __name__=="__main__":validate()
