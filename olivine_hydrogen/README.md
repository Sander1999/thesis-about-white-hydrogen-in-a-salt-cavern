# Olivine hydrogen and electricity model

This project connects 56 XRF readings to a finite iron inventory, batch hydrogen
partition/recovery, an explicit energy budget, a comparison of stationary,
recirculating and once-through operation, and native 3D DARTS transport through
an assumed olivine bed. The positive reaction rates are scenarios, not
measurements. The zero-production reference is essential at 50°C.

Open `Olivine_data_transformation.ipynb` for the full walkthrough. The notebook
uses the same simple explanatory style as the data analysis and executes both
the batch model, feedback comparison and native flow case. The thesis source is
`thesis/main.tex`.

## Run

From this directory or its parent:

```bash
bash /absolute/path/to/olivine_hydrogen/run_model.sh all
bash /absolute/path/to/olivine_hydrogen/run_model.sh check
bash /absolute/path/to/olivine_hydrogen/run_model.sh notebook
```

The launcher uses the existing shared interpreter
`/Users/sanderbertdacosta/.local/share/python-envs/darts-py311/bin/python`.
Set `OLIVINE_PYTHON` only to an independently verified compatible DARTS runtime.
It creates and cleans temporary plotting/Jupyter caches. No new environment is
created. The notebook executor selects this exact interpreter with a temporary
kernel; it does not use an unrelated default `python3` kernel. To work interactively
in Jupyter, select/register the same interpreter. Recorded versions are in
`environment.json` and `requirements.txt`; native DARTS requires the matching
installed engine, not merely a similarly named pip package.

`data`, `batch`, `feedback`, and `darts` can be run separately. `all` runs these
four steps, and `check` includes the feedback tests. `darts/config.json` controls
the spatial case. Batch settings are the `Setup` dataclass in `batch.py`; the
notebook demonstrates changing them with `dataclasses.replace`. The reduced
feedback settings are `Config` in `feedback/model.py`. Its full default suite
compares eight scenarios across three modes (24 cases) over 365 days; this duration differs
from the 72-hour native DARTS example.

## Main outcomes and boundaries

For 10 t fresh Eifel material, all Fe assumed ferrous and accessible:

- Magnetite-pathway capacity: 7.12 kg H₂; complete-ferric electron ceiling: 10.68 kg.
- Complete magnetite conversion: 6.31 kg collected, approximately 105 kWh gross
  electricity at 90% gas collection and 50% LHV electrical efficiency.
- Heating the baseline rock/brine from 20 to 50°C: approximately 765 kWh thermal.
  Electrical heating plus surface pressurization makes the selected standalone
  balance approximately −803 kWh even at complete conversion.
- Native flow case, assumed 365-day half-time over 72 h: 40.46 g H₂ generated,
  38.59 g exported dissolved and 1.87 g remaining dissolved. The 10-bar outlet
  separator gives no collected gas and no gross electricity in this case.

The batch water inventory is 20 t water plus 5.26 t salt. The DARTS alternative
instead injects that quantity over 72 h, in addition to approximately 1.51 t initial
porewater. These are different physical operating cases. Batch endpoints are
independent one-time recoveries, not cumulative withdrawals. The conditional warm
energy balance excludes heating and initial pressurization; recompressing fluid
after separation is an additional cost if it is recycled.

The Chabab correlation is restricted to 323.15–373.15 K, 10–230 bar, 0–5 molal
NaCl. Three-dimensional fields represent a porous bed, not an open salt cavern.
Permeability, accessible Fe, conversion time, density, viscosity and energy
efficiencies are declared assumptions. There is no pH/mineral-equilibrium,
coupled thermal or cavern-geomechanics prediction. The feedback comparison uses
an assumed passivation factor; this is not a calibrated mineral-reaction law.

## What the feedback comparison means

`feedback/` is a well-mixed finite-inventory model. It compares three operations:

- **Stationary:** retain the initial brine and hydrogen until a terminal flash.
- **Recirculating:** circulate the same fluid without renewal or gas stripping.
- **Once-through:** add fresh brine and withdraw an equal water mass. Each outlet
  parcel is flashed separately; dissolved exports are not automatically gas recovery.

All modes start with the same 10 t rock and 20 t water inventory. Flowing cases
circulate or replace 20 t water per three days, so their resource demands differ.
The comparison reports extra water, pumping and heating explicitly. The endpoint
recovery adds gas already collected from the outlet to one possible final flash;
terminal endpoints must not be summed as repeated withdrawals.

The no-feedback null isolates the analytic finite-inventory reaction law. The
intrinsic-limited case retains geometry feedbacks but is a null comparison for a flow-induced rate benefit:
mixing and product removal do not increase generation in this case. It still has
an assumed positive intrinsic reaction rate. The separate zero-source case gives
no hydrogen. Film-only and inhibition-only cases isolate the two proposed flow
effects. Additional scenarios combine them with pore closure, passivation or
greater accessible area. These test consequences of
hypotheses and do not show that circulation necessarily speeds hydrogen formation.

Expansion, cracking and porosity are scalar feedback assumptions. This component
does not simulate stress or propagating fractures and is not the native DARTS
model. The native three-dimensional transport calculation remains in `darts/`.
See `feedback/README.md` for its equations and model limits.

## Files and verification

- `data/raw/`: original measurement input; `data/processed/`: corrected summaries.
- `batch.py`, `run_batch.py`, `results/`: finite inventory, recovery and energy.
- `feedback/`: reduced operating-mode comparison, scenarios, ledgers and figures.
- `darts/`: native model, conservation checks, CSV fields, 3D figures and outputs.
- `tests/`: XRF and batch numerical checks.
- `thesis/`: complete editable LaTeX source and original sample photographs.

All 20 XRF/batch tests, 8 reduced-feedback tests and 3 DARTS property/rollback
tests pass, together with 30 native conservation/refinement checks. The notebook
has 11 executed code cells, no error outputs and 9 embedded figures. All 48 pages
of the thesis PDF were visually checked; mathematical references resolve.
Feedback checks include an analytic intrinsic limit, zero source, finite inventory,
hydrogen/water/salt and recovery ledgers, and time refinement. Its numerical
checks are separate from the native checks. The DARTS validation JSON includes a
5³ spatial sensitivity; it does not claim formal spatial convergence, because
the corner-cell completions move with the grid. Numerical checks do not validate
the assumed reaction rates.

Compile the thesis with an existing LaTeX installation, from `thesis/`:

```bash
tectonic --untrusted main.tex
# or latexmk -pdf main.tex
```

The `thesis/` folder is self-contained and can be uploaded as a LaTeX project.
After rerunning the model, synchronize its figures using the shared interpreter:

```bash
/Users/sanderbertdacosta/.local/share/python-envs/darts-py311/bin/python -B build_thesis.py --sync-only
```

Omit `--sync-only` to compile and replace the thesis PDF in the parent directory.
The entire model folder must remain available for rerunning the calculations. The cover illustration is
`thesis/figures/cover_concept.png`; it was produced with the built-in image tool.
Its exact prompt and conceptual scope are recorded in
`thesis/figures/cover_provenance.md`. All result graphs are computed from data or
simulation output, not image-generated.
