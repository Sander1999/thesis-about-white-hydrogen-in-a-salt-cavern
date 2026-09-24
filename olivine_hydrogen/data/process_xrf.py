"""Process thesis elemental readings without treating uncertainty as concentration.

The CSV contains compositional readings, not hydrogen-production experiments.
Reported elemental Fe can constrain conditional stoichiometric capacities;
it cannot identify ferrous fraction, reactive mineral abundance or kinetics.
No magnesium or mineralogical measurement is inferred from the Bal channel.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import re
import textwrap

import numpy as np

HERE = Path(__file__).resolve().parent
RAW_FILE = HERE / "raw" / "Olivine_data.csv"
REFERENCE_SHA256 = "81586a9f77f0a6c7797a7ba4fbd447f4828f7c157696b34e6a210bc535f6d16c"
EXPECTED_COUNTS = {"sample1": 5, "sample2": 5, "sample3": 6, "sample4": 15,
                   "sample5": 10, "sample6": 5, "sample7": 10}
REFERENCE_CHANNELS = ["Ag", "As", "Au", "Bal", "Bi", "Cd", "Co", "Cr", "Cu", "Fe", "Hf", "Hg", "Mn", "Mo", "Nb", "Ni", "Pb", "Pd", "Rb", "Re", "Sb", "Se", "Sn", "Sr", "Ta", "Th", "Ti", "U", "V", "W", "Y", "Zn", "Zr"]
METADATA = ["Reading No", "Reading Type", "Duration", "Main", "Low", "High", "Light",
            "Time", "User", "Sigma Value", "Units", "Batch", "Heat", "Lot", "Note", "Sample", "Mode"]
FE_KG_MOL = 0.055845
H2_KG_MOL = 0.00201588
MG_KG_MOL, SI_KG_MOL, O_KG_MOL = 0.024305, 0.0280855, 0.0159994
FORSTERITE_SIO2_MASS_FRACTION = (SI_KG_MOL + 2 * O_KG_MOL) / (2 * MG_KG_MOL + SI_KG_MOL + 4 * O_KG_MOL)


def parse_concentration(value: str) -> tuple[float | None, str, float | None]:
    """Return value, detection status and a stated numerical detection limit.

    '<LOD' has no numerical limit in this export. It remains censored with
    a missing value; it is never zero-filled or replaced by its uncertainty.
    """
    value = value.strip()
    if not value:
        return None, "missing", None
    if value.upper() == "<LOD":
        return None, "below_detection", None
    if value.startswith("<"):
        limit = float(value[1:].strip())
        if not np.isfinite(limit) or limit < 0:
            raise ValueError(f"Invalid detection limit: {value!r}")
        return None, "below_detection", limit
    number = float(value)
    if not np.isfinite(number) or number < 0:
        raise ValueError(f"Concentration must be finite and nonnegative: {value!r}")
    return number, "detected", None


def subgroup_for_reading(sample: str, reading_no: int) -> tuple[str, str]:
    """Reproduce original positional labels using stable read IDs.

    The notebook's row-order labels were not present in the instrument export.
    Mapping them to the original reading numbers makes reorderings safe but
    does not establish that the geological labels have been independently verified.
    """
    groups = {
        "sample1": [(4458, 4462, "all")], "sample2": [(4463, 4467, "all")],
        "sample3": [(4468, 4473, "all")],
        "sample4": [(4474, 4478, "big grains"), (4479, 4483, "small grains"), (4484, 4488, "mixed grains")],
        "sample5": [(4489, 4493, "weathered surface"), (4494, 4498, "fresh surface")],
        "sample6": [(4499, 4503, "all")],
        "sample7": [(4504, 4508, "dark big grains"), (4509, 4513, "light green grains")],
    }
    for first, last, label in groups.get(sample, []):
        if first <= reading_no <= last:
            return label, "Original notebook cell 6 row-order assumption, mapped to Reading No; not independently verified"
    return "unclassified", "No subgroup provenance for this reading"


def load_readings(path=RAW_FILE, *, validate_reference=True):
    """Read concentrations and their paired uncertainty fields separately.

    Reference mode checks the known 56-read sample counts, IDs and declared
    ppm units. Disable the count check for a deliberately different dataset;
    unknown read IDs remain unclassified rather than acquiring positional labels.
    """
    path = Path(path)
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        header = reader.fieldnames or []
        if len(set(header)) != len(header):
            raise ValueError("Duplicate CSV column names")
        missing = set(METADATA) - set(header)
        if missing:
            raise ValueError(f"Missing metadata columns: {sorted(missing)}")
        channels = [c for c in header if c and c not in METADATA and not c.endswith(" 2-Sigma")]
        if not {"Fe", "Bal"} <= set(channels):
            raise ValueError("Expected elemental Fe and Bal concentration channels")
        if any(not re.fullmatch(r"[A-Z][a-z]?|Bal", c) for c in channels):
            raise ValueError("Unrecognised concentration column; do not silently treat metadata as an element")
        if any(c + " 2-Sigma" not in header for c in channels):
            raise ValueError("Each concentration requires a separate reported 2-Sigma column")
        rows = []
        for line_number, source in enumerate(reader, 2):
            if None in source or any(v is None for v in source.values()):
                raise ValueError(f"Inconsistent CSV field count at line {line_number}")
            if not any(v.strip() for v in source.values()):
                continue
            if source["Units"].strip().lower() != "ppm":
                raise ValueError(f"Expected ppm units at line {line_number}")
            reading = int(source["Reading No"])
            sample = source["Sample"].strip()
            if not sample:
                raise ValueError(f"Missing sample at line {line_number}")
            subgroup, provenance = subgroup_for_reading(sample, reading)
            result = {"reading_no": reading, "csv_line": line_number, "sample": sample,
                      "subgroup": subgroup, "subgroup_provenance": provenance,
                      "timestamp": source["Time"].strip(), "units": "ppm",
                      "sigma_value_metadata": float(source["Sigma Value"]),
                      "duration_reported": float(source["Duration"]),
                      "duration_unit": "not stated by export; not reaction exposure",
                      "concentrations_ppm": {}, "reported_2sigma_ppm": {},
                      "detection_status": {}, "detection_limit_ppm": {}}
            for channel in channels:
                number, status, limit = parse_concentration(source[channel])
                uncertainty, uncertainty_status, _ = parse_concentration(source[channel + " 2-Sigma"])
                if uncertainty_status == "below_detection":
                    raise ValueError("Uncertainty fields cannot be used as censored concentrations")
                result["concentrations_ppm"][channel] = number
                result["reported_2sigma_ppm"][channel] = uncertainty
                result["detection_status"][channel] = status
                result["detection_limit_ppm"][channel] = limit
            rows.append(result)
    if not rows:
        raise ValueError("No readings found")
    ids = [r["reading_no"] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate Reading No values would double-count measurements")
    if validate_reference:
        if set(channels) != set(REFERENCE_CHANNELS):
            raise ValueError("Expected the 33 documented concentration channels with separate uncertainty pairs")
        if Counter(r["sample"] for r in rows) != EXPECTED_COUNTS:
            raise ValueError("Expected the reference 56 readings with the seven documented sample counts")
        if set(ids) != set(range(4458, 4514)) or any(r["subgroup"] == "unclassified" for r in rows):
            raise ValueError("Reading numbers or sample assignment differ from the subgroup provenance map")
    return rows, channels


def iron_capacity(fe_ppm, *, ferrous_fraction=1.0, accessible_fraction=1.0):
    """Conditional H2 capacities, never a measured yield or generation rate.

    Input is elemental Fe in mass ppm, not FeO. All Fe is assumed Fe(II)
    and accessible by default. Magnetite pathway: 3FeO+H2O -> Fe3O4+H2,
    so n(H2)=n(Fe)/3. Complete Fe(II)->Fe(III) electron-balance ceiling:
    n(H2)=n(Fe)/2. The latter is not a specified reaction or attainable yield.
    """
    fe = np.asarray(fe_ppm, dtype=float)
    if np.any(~np.isfinite(fe)) or np.any((fe < 0) | (fe > 1e6)):
        raise ValueError("Elemental Fe must lie between 0 and 1,000,000 mass ppm")
    for value in (ferrous_fraction, accessible_fraction):
        if not np.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("Ferrous/accessibility fractions must lie in [0, 1]")
    mol_fe_per_kg = fe * 1e-6 / FE_KG_MOL * ferrous_fraction * accessible_fraction
    return {
        "reactive_ferrous_Fe_mol_per_kg_rock": mol_fe_per_kg,
        "magnetite_H2_mol_per_kg_rock": mol_fe_per_kg / 3,
        "ferric_ceiling_H2_mol_per_kg_rock": mol_fe_per_kg / 2,
        "magnetite_H2_kg_per_tonne_rock": mol_fe_per_kg / 3 * H2_KG_MOL * 1000,
        "ferric_ceiling_H2_kg_per_tonne_rock": mol_fe_per_kg / 2 * H2_KG_MOL * 1000,
    }


def legacy_assumed_mg_number(fe_ppm, bal_ppm):
    """Counterfactual legacy Bal-derived Mg#, not a mineral measurement.

    Assumes ALL Bal is MgO+SiO2 with pure-forsterite mass proportions. In
    Mg2SiO4 the SiO2 fraction is 60.0843/140.6931=0.42706, not 0.406.
    No measured Mg is created: output is named as a legacy assumption only.
    Fe valence, accessories, oxide oxygen and the unknown Bal composition
    prevent interpreting this quantity as measured olivine Mg#.
    """
    fe, bal = np.broadcast_arrays(np.asarray(fe_ppm, float), np.asarray(bal_ppm, float))
    if np.any(~np.isfinite(fe + bal)) or np.any(fe < 0) or np.any(bal < 0):
        raise ValueError("Fe and Bal must be finite and nonnegative")
    assumed_mgo_mass = bal * 1e-6 * (1 - FORSTERITE_SIO2_MASS_FRACTION)
    assumed_mg_moles = assumed_mgo_mass / (MG_KG_MOL + O_KG_MOL)
    fe_moles = fe * 1e-6 / FE_KG_MOL
    denominator = assumed_mg_moles + fe_moles
    return np.divide(assumed_mg_moles, denominator,
                     out=np.full(fe.shape, np.nan), where=denominator > 0)


def element_statistics(rows, channel):
    """Keep read-to-read scatter and nominal instrument precision distinct.

    SEM=sample SD/sqrt(n_detected). The nominal instrument uncertainty on
    the mean assumes the header means 2σ and independent readings, dividing
    each reported value by two before propagation. Sigma Value=1.5 remains
    unresolved metadata. No covariance, calibration bias or independent
    geological replication is implied by this nominal estimate.
    """
    detected = [r for r in rows if r["detection_status"][channel] == "detected"]
    values = np.array([r["concentrations_ppm"][channel] for r in detected], dtype=float)
    n = len(values)
    sigmas = [r["reported_2sigma_ppm"][channel] for r in detected]
    instrument = np.array(sigmas, float) if n and all(v is not None for v in sigmas) else None
    sd = float(values.std(ddof=1)) if n > 1 else None
    return {
        "element": channel, "n_readings": len(rows), "n_detected": n,
        "n_below_detection": sum(r["detection_status"][channel] == "below_detection" for r in rows),
        "n_missing": sum(r["detection_status"][channel] == "missing" for r in rows),
        "mean_basis": "detected readings only; no imputation",
        "mean_ppm": float(values.mean()) if n else None,
        "sample_sd_ppm": sd, "sem_ppm": sd / np.sqrt(n) if sd is not None else None,
        "min_ppm": float(values.min()) if n else None, "max_ppm": float(values.max()) if n else None,
        "mean_reported_2sigma_ppm": float(instrument.mean()) if instrument is not None else None,
        "nominal_instrument_1sigma_mean_ppm": float(np.sqrt(np.sum((instrument / 2)**2)) / n) if instrument is not None else None,
    }


def summarize(rows, channels, *, by_subgroup=True):
    groups = defaultdict(list)
    for row in rows:
        groups[(row["sample"], row["subgroup"] if by_subgroup else "all readings")].append(row)
    summaries, long_table = [], []
    # Order by original read sequence, preserving physical subgroup order.
    for (sample, subgroup), group in sorted(groups.items(), key=lambda item: min(r["reading_no"] for r in item[1])):
        base = {"sample": sample, "subgroup": subgroup, "n_readings": len(group),
                "reading_numbers": " ".join(str(r["reading_no"]) for r in sorted(group, key=lambda r: r["reading_no"]))}
        all_stats = {c: element_statistics(group, c) for c in channels}
        for stats in all_stats.values():
            long_table.append(base | stats)
        for channel in ["Fe", "Bal"]:
            for key, value in all_stats[channel].items():
                if key not in ("element", "n_readings", "mean_basis"):
                    base[f"{channel}_{key}"] = value
        fe = all_stats["Fe"]
        base["Fe_mean_wt_percent"] = fe["mean_ppm"] / 1e4 if fe["mean_ppm"] is not None else None
        for key in ("sample_sd_ppm", "sem_ppm", "nominal_instrument_1sigma_mean_ppm"):
            base[f"Fe_{key.replace('_ppm', '_wt_percent')}"] = fe[key] / 1e4 if fe[key] is not None else None
        paired = [r for r in group if r["concentrations_ppm"]["Fe"] is not None and r["concentrations_ppm"]["Bal"] is not None and r["concentrations_ppm"]["Bal"] > 0]
        base["mean_readwise_Fe_to_Bal"] = float(np.mean([r["concentrations_ppm"]["Fe"] / r["concentrations_ppm"]["Bal"] for r in paired])) if paired else None
        base["legacy_assumption_Mg_number_mean"] = float(np.mean([legacy_assumed_mg_number(r["concentrations_ppm"]["Fe"], r["concentrations_ppm"]["Bal"]) for r in paired])) if paired else None
        for key in iron_capacity(0):
            base[key] = float(iron_capacity(fe["mean_ppm"])[key]) if fe["mean_ppm"] is not None else None
        for kind, stat in [("sample_sd", "sample_sd_ppm"), ("sem", "sem_ppm"),
                           ("nominal_instrument_1sigma_mean", "nominal_instrument_1sigma_mean_ppm")]:
            for path in ("magnetite", "ferric_ceiling"):
                base[f"{path}_H2_{kind}_kg_per_tonne_rock"] = (float(iron_capacity(fe[stat])[f"{path}_H2_kg_per_tonne_rock"])
                    if fe[stat] is not None else None)
        base["capacity_assumption"] = "All elemental Fe is ferrous and accessible; conditional stoichiometry, not observed H2"
        base["subgroup_provenance"] = "Original notebook positional labels mapped to Reading No; no independent laboratory confirmation"
        summaries.append(base)
    return summaries, long_table


def write_table(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_results(rows, summaries, output):
    """Save scientific figures as vector PDF and 180-dpi PNG."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with plt.rc_context({"font.size": 14, "axes.titlesize": 15, "axes.labelsize": 14,
                         "xtick.labelsize": 13, "ytick.labelsize": 13,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "pdf.fonttype": 42, "ps.fonttype": 42}):
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        labels = [s["sample"].replace("sample", "S") + ("\n" + textwrap.fill(s["subgroup"], 9) if s["subgroup"] != "all" else "") for s in summaries]
        position = np.arange(len(summaries))
        fig, axes = plt.subplots(2, 1, figsize=(11.7, 7.8), sharex=True, layout="constrained")
        fe = np.array([s["Fe_mean_wt_percent"] for s in summaries])
        axes[0].errorbar(position, fe, yerr=[s["Fe_sample_sd_wt_percent"] for s in summaries],
                         fmt="o", color="#176b87", elinewidth=1, capsize=4, label="Mean ± reading SD")
        axes[0].errorbar(position, fe, yerr=[s["Fe_nominal_instrument_1sigma_mean_wt_percent"] for s in summaries],
                         fmt="none", ecolor="#bc572b", elinewidth=2.5, capsize=2, label="Nominal instrument 1σ on mean")
        axes[0].set(ylabel="Reported elemental Fe [wt%]", title="Elemental composition and conditional hydrogen capacity")
        axes[0].legend(fontsize=13, loc="upper left", ncol=2)
        for offset, pathway, color, marker, label in [(-.13, "magnetite", "#176b87", "o", "Magnetite: n(H₂) = n(Fe)/3"),
                                                    (.13, "ferric_ceiling", "#aa5b31", "s", "Ferric ceiling: n(H₂) = n(Fe)/2")]:
            axes[1].errorbar(position+offset, [s[f"{pathway}_H2_kg_per_tonne_rock"] for s in summaries],
                             yerr=[s[f"{pathway}_H2_sample_sd_kg_per_tonne_rock"] for s in summaries],
                             fmt=marker, color=color, elinewidth=1, capsize=3, label=label)
        axes[1].set(ylabel="Conditional H₂ capacity\n[kg / tonne rock]", xticks=position, xticklabels=labels)
        axes[1].legend(fontsize=13, loc="upper left", ncol=2)
        capacity_top = max(s["ferric_ceiling_H2_kg_per_tonne_rock"] + s["ferric_ceiling_H2_sample_sd_kg_per_tonne_rock"] for s in summaries)
        axes[1].set_ylim(0, capacity_top * 1.25)
        for axis in axes:
            axis.grid(axis="y", alpha=.2)
            axis.margins(x=.06, y=.25)
        fig.supxlabel("All Fe assumed Fe(II) and accessible. Error bars show reading scatter;\nkinetic and mineralogical uncertainty are not included.", fontsize=13)
        for suffix in ("png", "pdf"):
            fig.savefig(output / f"fe_and_conditional_capacity.{suffix}", dpi=180)
        plt.close(fig)

        fig, axes = plt.subplots(2, 1, figsize=(11.7, 7.4), sharex=True, layout="constrained")
        for i, summary in enumerate(summaries):
            group = sorted([r for r in rows if r["sample"] == summary["sample"] and r["subgroup"] == summary["subgroup"]], key=lambda r:r["reading_no"])
            x = i + np.linspace(-.25, .25, len(group))
            fe = np.array([r["concentrations_ppm"]["Fe"] / 1e4 for r in group])
            sigma = np.array([r["reported_2sigma_ppm"]["Fe"] / 1e4 for r in group])
            axes[0].errorbar(x, fe, yerr=sigma, fmt="o", ms=4, color="#176b87", elinewidth=.9, capsize=2)
            axes[0].plot([i-.32, i+.32], [fe.mean()]*2, color="#aa5b31", lw=1.5)
            axes[1].bar(i-.15, summary["Fe_sample_sd_wt_percent"], width=.28, color="#176b87",
                        label="Reading SD" if i == 0 else None)
            axes[1].bar(i+.15, summary["Fe_mean_reported_2sigma_ppm"] / 2 / 1e4, width=.28, color="#aa5b31",
                        label="Mean nominal instrument 1σ per reading" if i == 0 else None)
        axes[0].set(ylabel="Elemental Fe [wt%]", title="Repeatability and variation between individual readings")
        axes[0].plot([], [], "o", color="#176b87", label="Reading ± reported 2-Sigma")
        axes[0].plot([], [], color="#aa5b31", label="Subgroup mean")
        axes[0].legend(fontsize=13, loc="upper left", ncol=2)
        axes[0].margins(y=.3)
        axes[1].set(ylabel="Fe spread / precision [wt%]", xticks=position, xticklabels=labels)
        axes[1].set_ylim(0, max(s["Fe_sample_sd_wt_percent"] for s in summaries) * 1.3)
        axes[1].legend(fontsize=13, loc="upper left", ncol=2)
        for axis in axes:
            axis.grid(axis="y", alpha=.2)
        fig.supxlabel("Subgroups follow the recorded reading order. These are repeated readings,\nnot independent reaction experiments.", fontsize=13)
        for suffix in ("png", "pdf"):
            fig.savefig(output / f"fe_reading_repeatability.{suffix}", dpi=180)
        plt.close(fig)


def process(path=RAW_FILE, output=HERE, *, make_plots=True, validate_reference=True):
    rows, channels = load_readings(path, validate_reference=validate_reference)
    output = Path(output)
    subgroup_summary, element_summary = summarize(rows, channels)
    sample_summary, _ = summarize(rows, channels, by_subgroup=False)
    flattened = []
    for row in rows:
        flat = {key: value for key, value in row.items() if not isinstance(value, dict)}
        for channel in channels:
            flat.update({f"{channel}_ppm": row["concentrations_ppm"][channel],
                         f"{channel}_reported_2sigma_ppm": row["reported_2sigma_ppm"][channel],
                         f"{channel}_status": row["detection_status"][channel],
                         f"{channel}_detection_limit_ppm": row["detection_limit_ppm"][channel]})
        fe = row["concentrations_ppm"]["Fe"]
        flat.update({key:float(iron_capacity(fe)[key]) if fe is not None else None for key in iron_capacity(0)})
        flattened.append(flat)
    for name, table in [("readings", flattened), ("subgroup_summary", subgroup_summary),
                         ("sample_summary", sample_summary), ("element_summary", element_summary)]:
        write_table(output / "processed" / f"{name}.csv", table)
    concentration_sum = [sum(value or 0 for value in row["concentrations_ppm"].values()) for row in rows]
    audit = {
        "input": str(Path(path).resolve()), "input_sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
        "reference_sha256": REFERENCE_SHA256, "reference_count_validation": validate_reference,
        "rows": len(rows), "sample_counts": dict(Counter(r["sample"] for r in rows)),
        "subgroup_counts": {s["sample"] + ":" + s["subgroup"]:s["n_readings"] for s in subgroup_summary},
        "concentration_channels": channels, "units": "ppm (mass basis assumed from elemental export)",
        "sum_reported_ppm_including_Bal": [min(concentration_sum), max(concentration_sum)],
        "Bal_status": "Balance channel of unspecified composition; never a measured Mg/mineral fraction",
        "reported_sigma_metadata_values": sorted(set(r["sigma_value_metadata"] for r in rows)),
        "uncertainty_convention": "2-Sigma headers divided by two ONLY for nominal instrument estimates; Sigma Value=1.5 unresolved. Independence/covariance/calibration not established.",
        "capacity_stoichiometry": {"magnetite_H2_per_Fe_mol": 1/3, "full_ferric_ceiling_H2_per_Fe_mol": 1/2,
                                  "Fe_molar_mass_kg_mol": FE_KG_MOL, "H2_molar_mass_kg_mol": H2_KG_MOL,
                                  "ferrous_fraction": 1.0, "accessible_fraction": 1.0},
        "legacy_Mg_number_status": "Assumes all Bal is MgO+SiO2 in pure-forsterite proportions; counterfactual diagnostic, not measured Mg#",
        "legacy_pure_forsterite_SiO2_mass_fraction": FORSTERITE_SIO2_MASS_FRACTION,
        "corrected_notebook_issues": ["Concentration and 2-Sigma channels were pooled in original cells1/5; now separate.",
                                      "Original cell1 ratio plot selected only the previous loop's final sample; now all groups exported.",
                                      "Original subgroup labels depended on current row order; now stable original read-ID mapping with explicit unverified provenance.",
                                      "<LOD remains censored with no invented concentration or numerical limit.",
                                      "SD, SEM and nominal instrument uncertainty are reported separately."],
        "missing_for_kinetic_calibration": ["H2 versus reaction time", "reaction T/P history", "rock/fluid mass and volumes", "Fe valence and mineral phases", "reactive surface area", "blanks/standards/independent experiments"],
        "inference_limit": "Theoretical capacity only. No measured H2 yield, fitted reaction rate or reservoir source flux is obtained from this dataset.",
    }
    (output / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    if make_plots:
        plot_results(rows, subgroup_summary, output / "figures")
    return audit, subgroup_summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=RAW_FILE)
    parser.add_argument("--output", type=Path, default=HERE)
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--allow-different-readings", action="store_true",
                        help="Disable reference sample/reading count checks; unknown IDs remain unclassified")
    args = parser.parse_args()
    audit, _ = process(args.input, args.output, make_plots=not args.no_plots,
                       validate_reference=not args.allow_different_readings)
    print(f"Processed {audit['rows']} readings into {len(audit['subgroup_counts'])} explicit subgroups: {args.output.resolve()}")
    print("Capacity assumes all elemental Fe is ferrous and accessible; no H2 rate or measured yield was fitted.")


if __name__ == "__main__":
    main()
