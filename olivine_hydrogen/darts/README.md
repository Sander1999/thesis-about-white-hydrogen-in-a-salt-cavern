# Assumed packed-bed DARTS alternative

This is a separate **flow-through packed-bed scenario**, not the thesis batch
reactor and not Darcy flow through an open salt cavern. Its permeability,
porosity, well layout, phase mobilities and conversion half-time are assumptions.
XRF constrains the iron inventory; it supplies no hydrogen production rate.

Run from the reservoir-simulation project root, with the existing shared DARTS
runtime selected by the launcher:

```bash
./run_darts.sh olivine_hydrogen/darts/run.py --suite
./run_darts.sh olivine_hydrogen/darts/test_model.py
```

`run_case(parameters=None, output=None, plots=True)` in `run.py` returns
`(summary, history)`. JSON overrides are supported by `--config path.json`.
The suite writes `output/reactive/` and `output/zero_source/`; each has a status
file, summary, accepted-timestep CSV, final spatial CSV, plots and native logs.
`spatial_final.csv` includes actual native Darcy-flow vectors and hydrogen
concentration; `packed_bed_3d.png` plots these fields.

## Model and assumptions

- Dry rock mass: 10,000 kg; grain density: 3,300 kg/m³; porosity: 0.35.
  A 2×2 m cross-section gives a 1.16550 m height, 4.66200 m³ bulk volume and
  1.63170 m³ pore volume. The initial pores contain approximately 1,512 kg water
  plus 398 kg NaCl. This initial inventory is **additional to** the 20,000 kg water
  feed over 72 h; salt accompanies that feed at 4.5 mol/kg water.
- A 3×3×3 finite-volume grid has an injector and producer in opposite corner
  cells. Permeability is 1,000 mD in every direction. Producer BHP is 120 bar
  absolute; reservoir pressure exceeds this under flow. Gravity is active;
  capillary pressure and molecular diffusion are omitted. Relative
  permeability is the square of each phase saturation. Spatial parameters
  are illustrative; the mesh is not a resolved description of real grain packing.
- Isothermal 50 °C; brine density 1,170 kg/m³ at 120 bar, compressibility
  4×10⁻⁵/bar, viscosity 1 cP are engineering assumptions. Hydrogen density uses
  the pressure-explicit NIST correlation; gas viscosity is 0.0094 cP.
  Water and salt are treated as nonvolatile.
- H₂/H₂O/NaCl are independent transported components. The empirical Chabab
  correlation partitions H₂ between brine and gas, using salt molality and
  **salt-free** aqueous mole fraction. Its stated domain is 323.15–373.15 K,
  10–230 bar and 0–5 molal. Accepted states are checked against these limits.
  No extra gas-fugacity multiplier is applied to the empirical fit.
- The default elemental Fe mass fraction 0.0591870406 is the fresh-Eifel scenario.
  All Fe is assumed ferrous and accessible, an explicit capacity ceiling.
  The balanced redox proxy is 3 FeO + H₂O → Fe₃O₄ + H₂. It yields 1/3 mol H₂ per mol Fe,
  giving 7.12173 kg H₂ capacity. It is not a complete olivine reaction network.
  The 365-day conversion half-time is exploratory, not measured at these
  conditions. pH, passivation, mineral equilibrium and reaction affinity
  are not modelled. Porosity/permeability remain fixed; solid oxygen uptake
  is recorded without predicting mineral-volume change.

## Conservation and interpretation

The source consumes water and finite accessible iron. The native source rate is
constant over an attempted timestep and capped by available water/iron. Only a
converged timestep commits mineral consumption. Its exact first-order integral
provides an independent check; rejected attempts consume nothing. A separate
solid oxygen ledger closes the mass transfer from water to reaction products.

Native mass balance uses interpolated accumulation and independently reconstructed
TPFA perforation fluxes, integrated at accepted backward-Euler timestep endpoints.
The legacy native component-rate display is not used for these balances.
Exact-property balances are also reported to expose interpolation error. A tiny
numerical H₂ composition floor is included explicitly in the initial and injected
inventories; the zero-source control records this background separately.

For the default 72 h case, approximately 40.46 g H₂ is generated, 38.59 g leaves in
brine and 1.87 g remains dissolved. No free gas is predicted. **Aqueous export is
not recovered gas**. Each accepted timestep independently flashes the total
exported stream at 10 bar and 50°C. A 90% gas-collection efficiency and 50%
fuel-cell efficiency at 33.33 kWh/kg H₂ apply only to the gas released. The
default stream is below saturation even at 10 bar: collected gas and gross
electricity are both zero. This flow-through separator calculation is separate
from the batch reactor, whose product inventory and fluid volume differ. Gross
electricity, when positive, excludes pumping, heating and gas-processing demands.

`test_model.py` checks published solubility data, analytical flash/salt closure,
rejected-step inventory handling, native component/Fe/total-mass balances,
finite inventory and gas appearance under a deliberately rapid reaction stress
case, and timestep/interpolation refinement. The stress case is a numerical test,
not a rate forecast. These checks do not establish spatial convergence or field
calibration. A 5×5×5 spatial sensitivity reports export, storage and pressure
changes. It also changes corner-cell well locations and completion lengths, so
it is not a controlled spatial-convergence demonstration.

The reaction equations are
`capacity_H2_kmol = accessible_Fe_kmol / 3` and
`delta_extent = remaining_capacity * (1 - exp(-ln(2) * dt / half_time))`.
The per-cell source is capped by available water and then divided by bulk cell
volume and timestep to obtain kmol/(m³ day). The native sink vector is
`[-1, +1, 0]` for H₂/H₂O/NaCl per kmol of reaction extent. Only an accepted
timestep subtracts `3 * delta_extent` from accessible Fe. Dissolved hydrogen
uses `xH2 = nH2 / (nH2 + nH2O)`; NaCl is excluded from this denominator but
remains a transported and conserved component.

References: Chabab et al. (2020),
[doi:10.1016/j.ijhydene.2020.08.192](https://doi.org/10.1016/j.ijhydene.2020.08.192);
Lemmon et al. (2008),
[NIST hydrogen-density correlation](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=832233).
