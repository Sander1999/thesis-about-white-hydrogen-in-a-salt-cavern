# Reduced expansion, cracking and flow comparison

This is an explicitly assumed **well-mixed compartment model**, separate from
the native spatial DARTS model in `../darts`. It asks how selected reaction,
surface-access and hydrogen-removal hypotheses change the comparison between
stationary brine, closed recirculation and once-through brine renewal. It does
not compute fracture stresses, crack propagation, mineral equilibrium or actual
salt-cavern deformation. XRF constrains the total-iron ceiling, not these rates.

Run from the repository root using the existing shared environment:

```sh
PYTHONDONTWRITEBYTECODE=1 /Users/sanderbertdacosta/.local/share/python-envs/darts-py311/bin/python -m olivine_hydrogen.feedback.run
PYTHONDONTWRITEBYTECODE=1 /Users/sanderbertdacosta/.local/share/python-envs/darts-py311/bin/python -m olivine_hydrogen.feedback.validate
```

Use `--days 3` or `--days 1095` for another duration; default is 365 days.
For controlled scenarios, import `Config`, use `dataclasses.replace`, and call
`simulate(mode, cfg)` or `run_suite(output, cfg)`. No additional installation is
needed. `simulate` returns arrays plus a summary; `run_suite` writes individual
CSV trajectories, `comparison_endpoints.csv`, `summary.json` and two PNG figures.
`validate` records its tests and timestep comparison in `output/validation.json`.

## Shared inventory and operating boundaries

All three modes start with 10,000 kg rock, 20,000 kg water and 4.5 mol NaCl per kg
water at imposed 323.15 K and 120 bar absolute. Initial brine is about 25.26 t.
The fresh-Eifel mean Fe mass fraction is 0.0591870406. Assuming all Fe is ferrous
and accessible gives a 7.121734 kg H2 ceiling through the redox proxy
`3 FeO + H2O -> Fe3O4 + H2`. That assumption is deliberately favourable.

- **Stationary:** no external flow; all hydrogen stays in the reactor until a
  separately evaluated terminal separator operation.
- **Recirculating:** the same closed inventory circulates at a reference rate
  equivalent to 20 t water per three days. It can change the prescribed film
  factor. It does not renew water, strip hydrogen, or dilute the inventory.
- **Once-through:** starts with the same 20 t water and receives **additional**
  20 t water every three days, plus feed salt. Equal water mass is withdrawn;
  the associated outlet salt follows the current reactor molality. Thus total
  brine mass flow is not assumed identical when salinities differ. Additional
  water reaches **2.433 million kg in one year**, with about 0.640 million kg
  additional salt. This is not a comparison using equal lifetime water demand.

The liquid pool includes the surrounding reactor/loop fluid. At rock density
3300 kg/m3 and initial bed porosity 0.35, the reactive bed occupies only 4.662 m3,
with 1.632 m3 initial pores; the roughly 21.59 m3 liquid inventory does not all
fit in those pores. Effective porosity governs the assumed **reactive-bed**
access and hydraulic resistance. Pressure remains externally imposed and the
model does not solve rigid-vessel gas pressure or expel liquid as pores close.

## Equations and scenario parameters

Let `Nmax` be the H2-producing reaction capacity in mol and `xi` its extent.
Only in this model we stipulate congruent whole-rock alteration
`f = xi/Nmax`. Neither that equivalence nor the fraction of actual olivine is
measured by XRF. Rock alteration, Fe oxidation and Mg hydration are not generally
equivalent processes.

```text
e = expansion * (1 - phi0) * f
opening = 1 - exp(-max(e - threshold, 0) / transition)
phi = max(phi_floor, phi0 - e*(1 - relief*opening))
area = 1 + (area_max - 1)*opening
passivation = exp(-passivation_coefficient*f)
access = ((phi - phi_floor)/(phi0 - phi_floor))**access_exponent
K/K0 = (phi/phi0)**3 * ((1-phi0)/(1-phi))**2
film = 1/(1 + film_Da/mixing)
inhibition = 1/(1 + inhibition_coefficient * dissolved_H2/saturation_capacity)
dxi/dt = ln(2)/half_time * (Nmax - xi) * area*passivation*access*film*inhibition
```

The intrinsic half-time is 365 days before these multipliers. The `rate_multiplier`
parameter is zero for the zero-source control. Defaults are expansion 0.4,
`phi0=0.35`, `phi_floor=0.01`, opening threshold 0.03 and transition width 0.05
(both fractions of initial bulk volume), relief 0.5, area cap 5, passivation
coefficient 2 and access exponent 1. Mixing is 1 when stationary and 4 when
flowing. These are scenario choices, not stress-derived laws or fitted values.
Cracking only relieves some expansion; it cannot increase porosity above its
initial value. The area cap is independent of the porosity relation.

The intrinsic-limited default sets `film_Da=0` and inhibition coefficient zero:
flow does not accelerate the source in that comparison. Other retained runs are:

| Scenario | Changes relative to the default |
|---|---|
| `no_feedback_null` | Expansion, passivation, access feedback disabled; area 1; analytic first-order source |
| `film_only_hypothesis` | Film Da 1 |
| `inhibition_only_hypothesis` | Product inhibition coefficient 1 |
| `contact_inhibition_hypothesis` | Both coefficients 1 |
| `closure_hypothesis` | Both coefficients 1; no relief, area 1, passivation 4 |
| `opening_hypothesis` | Both coefficients 1; relief 0.75, passivation 0.5 |
| `zero_source` | Source multiplier 0 |

Useful exploratory brackets are expansion 0–0.6, relief 0–1, area factor 1–5,
passivation 0–4, film Da 0–3, mixing 1–4, and inhibition coefficient 0–3.
They are neither confidence intervals nor validated parameter bounds. Pure
forsterite expansion is a geometric motivation; no specific product mineral
assemblage is calculated. The onset and width are especially unmeasured.

Water loss combines redox water and an extra hydration scenario:

```text
water_consumed = xi*0.01801528 + 0.15*rock_kg*f
solid_mass_gain = water_consumed - generated_H2_kg
```

The default extra 0.15 kg water/kg fully altered rock represents hypothetical
hydration in addition to the Fe-redox proxy, not a proven balanced whole-mineral
reaction. Set it to zero to omit additional hydration. Other `water_mol_per_h2`
choices are supported explicitly; the default redox stoichiometry is one.
NaCl is conserved except for measured model inlet/outlet fluxes. The model raises
an error if water loss takes brine outside the adopted solubility domain.

## Partition, collection and costs

Partition imports the primary Chabab et al. (2020) correlation from `../batch.py`
([DOI](https://doi.org/10.1016/j.ijhydene.2020.08.192)), valid only at
323.15–373.15 K, 10–230 bar and 0–5 molal. Gas is dry H2; water vapour and other
dissolved gases are omitted. The phenomenological inhibition uses dissolved
saturation fraction, not reaction Gibbs energy or a calculated pH/Eh.

The once-through outlet withdraws **liquid only**. Free gas stays in the reactor.
Each instantaneous effluent parcel equilibrates separately at 10 bar and 50 C;
90% of any liberated H2 is collected. Early dilute effluent is never pooled to
create an artificial rich stream. Exported dissolved hydrogen left in the
separator brine remains in the loss/recovery ledger. There is no membrane,
vacuum degassing or gas-entrainment mechanism.

At each reported endpoint, the remaining reactor contents can receive a single
10-bar terminal flash. `endpoint_collected_h2_kg` is that independent potential,
not ongoing delivery or a cumulative sum over earlier endpoints.
`total_deliverable_h2_kg` adds prior outlet collection to this terminal potential.
The 1 kg crossing times in summaries are interpolated from sampled trajectories.

Imposed flow requires `loop_dp=5 bar/(K/K0)`. Renewal additionally charges fresh
brine from 1 to 120 bar; recirculation pays only loop resistance. This is a
homogeneous-bed Kozeny-Carman cost estimate, not calibrated fracture hydraulics.
The model reports if the loop-head diagnostic exceeds 200 bar; it does not claim
a pump can sustain that head or automatically reduce the requested flow. Initial
charging is counted once. New feed heating from 20 to 50 C is charged only to
renewal; the closed loop returns fluid at reactor conditions. Pump efficiency is
70%, heater efficiency 95%, and default heat recovery is zero. Initial fluid and
rock heating follows the engineering assumptions of the original batch model.
Fuel-cell electricity uses 50% of the 33.33 kWh/kg LHV. Reported net balances
exclude heat loss, mining, crushing, solids handling, chemical treatment, pressure
vessel costs, reaction heat and other unquantified duties. No claim of net-positive
field operation follows from these scenario calculations.

## Verification and interpretation

Checks cover H2/water/salt balances, the solid mass gain, finite Fe capacity,
zero source, no iron, bounded area and porosity, the analytic no-feedback case,
the additional-water ledger, parcel-by-parcel recovery, property-domain guards
and timestep halving. Successful numerical checks do not validate feedback laws.

Flow can increase **generation** when film resistance or product inhibition is
assumed. The same renewal can decrease **retained concentration** and prevent
10-bar gas recovery by exporting dilute brine. Compare generation, retention,
collection and costs separately; the model contains no universal winner.
