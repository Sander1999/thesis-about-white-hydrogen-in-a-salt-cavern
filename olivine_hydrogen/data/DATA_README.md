# Olivine composition and conditional hydrogen capacity

These data describe **elemental composition**, not hydrogen production. The
56 readings from seven sample IDs can support composition summaries and
conditional stoichiometric capacity estimates. They cannot calibrate a
hydrogen-generation rate or activation energy.

## Inputs and reproduction

`raw/Olivine_data.csv` is an unchanged copy of the supplied thesis file
`bachelor/year 3 TU delft/scriptie/Olivine data.csv`. Its SHA-256 is
`81586a9f77f0a6c7797a7ba4fbd447f4828f7c157696b34e6a210bc535f6d16c`.
The original notebooks and exports were inspected read-only and not executed.

From the reservoir simulation project:

```bash
./run_darts.sh olivine_hydrogen/data/process_xrf.py
./run_darts.sh -m unittest olivine_hydrogen.tests.test_xrf -v
```

The launcher uses the existing shared interpreter
`/Users/sanderbertdacosta/.local/share/python-envs/darts-py311/bin/python`
and removes its temporary plotting cache. No installation was required.
Verified imports: Python 3.11.16, NumPy 2.4.6, SciPy 1.17.1,
Matplotlib 3.11.2 and pandas 2.3.3; `pip check` found no broken requirements.
The processor itself uses standard Python, NumPy and Matplotlib.

The default parser validates the 56 readings, known sample counts,
33 concentration channels and their separate uncertainty fields, unique
read IDs, and ppm units. `--allow-different-readings` permits a different
dataset; unfamiliar IDs stay unclassified. Raw concentration values are
never modified. `--output` chooses another result directory.

## What was measured

The export records readings 4458–4513 on 20 June 2025 from 03:47:06 to
04:45:54. `Reading Type` and `Mode` are `Mining`; `Units` is `ppm`.
The elemental concentrations are interpreted on a mass basis. The file
does not identify the instrument model, calibration standards or sample
origin. Its acquisition timestamps and roughly 20-unit `Duration` values
are not reaction exposure times; the duration unit is not stated.

There are 32 elemental channels plus `Bal`. Fe, Cr, Mn, Ni, Pb and Zn
are detected in every reading. Mg, Si, O and Fe oxidation state are absent.
Fe ranges from **4.0883 to 7.8708 wt%** across individual readings. These
are spot/repeat readings grouped by sample ID, not 56 independent reaction
experiments.

The original notebook assigns these groups by row order:

| Sample | Original read IDs | Label |
|---|---|---|
| sample1 | 4458–4462 | all |
| sample2 | 4463–4467 | all |
| sample3 | 4468–4473 | all |
| sample4 | 4474–4478 / 4479–4483 / 4484–4488 | big / small / mixed grains |
| sample5 | 4489–4493 / 4494–4498 | weathered / fresh surface |
| sample6 | 4499–4503 | all |
| sample7 | 4504–4508 / 4509–4513 | dark big / light green grains |

These labels come from zero-based cell 6 of
`Olivine_data_transformation.ipynb`, not the raw instrument metadata.
The processor maps the original labels to stable read IDs so file reordering
does not silently change them. Their geological interpretation still requires
the laboratory notes. Big/small are qualitative labels, not grain diameters
or measured reactive surface areas.

## Corrections and uncertainty

The original notebook's cell 1 melts both concentration and `2-Sigma`
columns, removes the suffix, and averages them as if both were concentrations.
Cell 5 reuses this table. For sample1 Fe this changes the correct mean
58,452.34844 ppm to 29,311.94285 ppm in the plotted quantity. Uncertainty
values can even create positive bars for elements reported entirely below
detection. The corrected processor keeps the two channel types separate.
The original ratio-plot loop also reuses its last sample variable and displays
only sample7 groups. All groups are exported here.

The Fe/Bal calculations in original cells 2/4/6 avoid that channel-mixing
error. Their eleven means reproduce to floating-point precision. The new
tables add read counts, IDs and uncertainty information:

- `sample_sd_ppm`: sample standard deviation across detected readings,
  using `ddof=1`. This includes heterogeneity and measurement variation.
- `sem_ppm`: SD/√n. It describes the mean under an independent-reading
  assumption; it is not independent geological replication.
- `mean_reported_2sigma_ppm`: the instrument export's named uncertainty
  channel, preserved separately.
- `nominal_instrument_1sigma_mean_ppm`: √Σ(u₂/2)²/n. This explicitly
  assumes that the column header means 2σ and that errors are independent.
  The export also gives `Sigma Value = 1.5`; that convention remains
  unresolved. Calibration bias, covariance and representativeness are not
  included. This nominal estimate is not an independently established
  confidence interval.

`<LOD` remains a censored value with no imputed concentration. This export
does not supply numerical detection limits. If a future input contains a
literal threshold such as `<3.5`, the processor retains that limit separately.
Trace-element means are explicitly **means of detected readings only** and
may therefore be biased. A channel with no detections has no mean, even if
its uncertainty column contains positive numbers. No zero imputation is used.

## Conditional capacities, not a fitted yield

Elemental Fe mass fraction is `Fe_ppm × 10⁻⁶`. With Fe molar mass
0.055845 kg/mol, the iron inventory per kg rock is:

\[
n_{Fe}=\frac{Fe_{ppm}\,10^{-6}}{0.055845}\quad\mathrm{mol/kg}.
\]

Two transparent end-member calculations are supplied:

- **Magnetite pathway:** `3FeO + H₂O → Fe₃O₄ + H₂`, giving
  `nH₂ = nFe/3` for the Fe participating in that reaction.
- **Complete ferric electron-balance ceiling:** all Fe(II) becomes Fe(III),
  giving `nH₂ = nFe/2`. This is an electron-budget ceiling, not a specified
  mineral reaction or evidence that the yield is attainable.

The default calculation assumes **all reported elemental Fe is Fe(II) and
accessible**. Neither assumption is measured here. The reusable
`iron_capacity` function accepts independent ferrous and accessibility
fractions in [0, 1] for stated sensitivities. These cannot be fitted from
the supplied composition file.

Multiplication by 0.00201588 kg/mol and 1000 kg rock/tonne converts mol/kg
to kg H₂/tonne rock. As an independent units check, 55,845 ppm Fe is exactly
one mole Fe/kg rock: the two capacities are **0.67196** and **1.00794 kg H₂
per tonne rock**. The figures propagate Fe reading scatter linearly, while
leaving unknown mineralogical and chemical uncertainty unquantified. Total
Fe is not automatically ferrous olivine; it can occur in accessory or already
oxidized phases. No output has a time denominator.

## Legacy Bal-based Mg number

`Bal` plus detected elemental channels totals one million ppm to within
0.054 ppm, consistent with an instrument balance channel. It is not a
measured Mg, Si, O or olivine fraction. `Fe/Bal` is therefore a descriptive
ratio to an unspecified remainder, not an independent mineralogical ratio.

For comparison with the original thesis, the output includes only the
clearly named `legacy_assumption_Mg_number_mean`. Its counterfactual
assumption is that **all Bal has pure-forsterite composition**:

\[
f_{Mg}=\frac{2(24.305)}{140.6931}=0.345504,\quad
f_{SiO_2}=\frac{60.0843}{140.6931}=0.427059,
\]

\[
Mg_{assumed}=Bal\,f_{Mg},\qquad
Mg\#_{assumed}=\frac{Mg_{assumed}/24.305}
{Mg_{assumed}/24.305+Fe/55.845}.
\]

The SiO₂ mass fraction is about **0.427, not 0.406**. These calculations do
not resolve accessory phases, iron valence or oxide oxygen. No measured Mg
column is generated, and this legacy estimate does not enter the H₂ capacity.

## Retained outputs

- `processed/readings.csv`: each concentration, reported uncertainty,
  detection status, optional limit, read ID and conditional Fe capacity.
- `processed/subgroup_summary.csv`: eleven groups, Fe/Bal summaries,
  capacity end members and separate variability/precision estimates.
- `processed/sample_summary.csv`: seven sample-ID summaries, pooling the
  within-sample subgroups explicitly as `all readings`.
- `processed/element_summary.csv`: all concentration channels by subgroup,
  including detection counts and means with no imputation.
- `figures/fe_and_conditional_capacity.png` and `.pdf`: composition and
  capacity with reading variability; the PDF retains vector graphics.
- `figures/fe_reading_repeatability.png` and `.pdf`: individual readings,
  reported uncertainty and observed scatter compared with nominal precision.
- `audit.json`: input fingerprint, counts, assumptions and limits.

Ten tests cover the true Fe mean, channel separation, censoring, original
group counts and stable labels, rejected units/duplicates, independent
stoichiometric conversions, capacity scaling and absence of invented
measured Mg. Test fixtures are created in temporary directories and removed.

A kinetic model additionally needs H₂ versus elapsed reaction time, sample
and fluid masses, reaction temperature/pressure history, gas-volume and
sampling corrections, mineralogy and Fe valence, surface area, fluid
chemistry, blanks and independent repeat experiments. The supplied Gantt
notebook is a project schedule and supplies none of those measurements.
