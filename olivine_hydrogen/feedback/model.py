"""Finite-inventory well-mixed reactor with explicitly assumed feedbacks.

The existing native DARTS bed resolves spatial flow. This separate compartment
model examines uncalibrated surface/closure/contact hypotheses, not stresses,
fracture propagation, equilibrium mineral assemblages or cavern mechanics.
Units: mol, kg, m3, bar, days, kWh. Pressure/temperature are imposed boundaries.
"""
from dataclasses import dataclass, asdict
import math
import numpy as np
from scipy.integrate import solve_ivp
from olivine_hydrogen.batch import (Setup, capacity_mol, solubility_mol_kg,
                                   M_H2, M_WATER, M_SALT, M_FE)

MODES = ("stationary", "recirculating", "once_through")
STATE_NAMES = ("extent_mol", "reactor_h2_mol", "water_kg", "salt_mol",
               "export_h2_mol", "export_water_kg", "export_salt_mol",
               "outlet_separator_gas_mol", "flow_pump_kwh", "feed_water_kg",
               "feed_salt_mol", "circulated_water_kg")


@dataclass(frozen=True)
class Config(Setup):
    fe_mass_fraction: float = 0.0591870406
    half_time_days: float = 365.
    rate_multiplier: float = 1.
    duration_days: float = 365.
    max_step_days: float = 1.
    relative_tolerance: float = 2e-9
    absolute_tolerance: float = 1e-9
    initial_porosity: float = .35
    porosity_floor: float = .01
    rock_density_kg_m3: float = 3300.
    solid_expansion_fraction: float = .4
    # Water for hypothetical Mg-silicate hydration IN ADDITION to Fe-redox water.
    hydration_water_kg_per_kg_rock: float = .15
    crack_threshold_bulk_fraction: float = .03
    crack_transition_bulk_fraction: float = .05
    crack_relief_fraction: float = .5
    maximum_area_factor: float = 5.
    passivation_coefficient: float = 2.
    porosity_access_exponent: float = 1.
    film_damkohler: float = 0.
    flow_mixing_factor: float = 4.
    product_inhibition_coefficient: float = 0.
    water_flow_kg_day: float = 20_000./3.
    reference_loop_dp_bar: float = 5.
    diagnostic_head_limit_bar: float = 200.

    def validate(self):
        super().validate()
        if not 0 < self.porosity_floor < self.initial_porosity < 1:
            raise ValueError("Require 0 < porosity floor < initial porosity < 1")
        if not 1 <= self.maximum_area_factor <= 5:
            raise ValueError("Scenario area factor must stay between one and five")
        if not 0 <= self.crack_relief_fraction <= 1:
            raise ValueError("Cracking relief must lie between zero and one")
        positive = ("half_time_days", "duration_days", "max_step_days", "rock_density_kg_m3",
                    "crack_transition_bulk_fraction", "relative_tolerance", "absolute_tolerance",
                    "diagnostic_head_limit_bar")
        nonnegative = ("rate_multiplier", "solid_expansion_fraction", "hydration_water_kg_per_kg_rock",
                       "crack_threshold_bulk_fraction", "passivation_coefficient", "film_damkohler",
                       "porosity_access_exponent", "product_inhibition_coefficient",
                       "water_flow_kg_day", "reference_loop_dp_bar")
        if any(getattr(self, key) <= 0 for key in positive) or any(getattr(self, key) < 0 for key in nonnegative):
            raise ValueError("Invalid feedback, integration or flow parameter")
        if self.flow_mixing_factor < 1:
            raise ValueError("Mixing factor must be at least one")
        return self


def feedback_factors(f, cfg):
    """f is a stipulated congruent whole-mineral alteration fraction.

    Cracking relieves part of expansion; it cannot create porosity above phi0.
    The area cap is an independent hypothesis, not inferred from porosity.
    """
    f = np.clip(np.asarray(f), 0., 1.)
    loading = cfg.solid_expansion_fraction*(1-cfg.initial_porosity)*f
    opening = -np.expm1(-np.maximum(loading-cfg.crack_threshold_bulk_fraction, 0.) /
                       cfg.crack_transition_bulk_fraction)
    porosity = np.maximum(cfg.porosity_floor, cfg.initial_porosity -
                          loading*(1-cfg.crack_relief_fraction*opening))
    area = 1+(cfg.maximum_area_factor-1)*opening
    passivation = np.exp(-cfg.passivation_coefficient*f)
    access = ((porosity-cfg.porosity_floor)/(cfg.initial_porosity-cfg.porosity_floor))**cfg.porosity_access_exponent
    # Kozeny-Carman ratio is only a homogeneous-bed engineering diagnostic.
    phi0 = cfg.initial_porosity
    permeability_ratio = (porosity/phi0)**3*((1-phi0)/(1-porosity))**2
    return {"altered_fraction_proxy": f, "expansion_loading": loading,
            "crack_opening_fraction": opening, "porosity": porosity,
            "area_factor": area, "passivation_factor": passivation,
            "pore_access_factor": access, "permeability_ratio": permeability_ratio}


def water_per_extent(cfg):
    nmax = capacity_mol(cfg)
    extra = cfg.hydration_water_kg_per_kg_rock*cfg.rock_kg/nmax if nmax > 0 else 0.
    return cfg.water_mol_per_h2*M_WATER + extra


def _partition(nh2, water, salt, cfg, pressure=None):
    if water <= 0:
        raise ValueError("Water exhausted: scenario no longer has a valid brine phase")
    molality = max(salt, 0.)/water
    concentration = float(solubility_mol_kg(cfg.temperature_k,
                          cfg.pressure_bar if pressure is None else pressure, molality))
    dissolved = min(max(nh2, 0.), water*concentration)
    return dissolved, max(nh2-dissolved, 0.), concentration, molality


def _rhs(mode, cfg):
    nmax = capacity_mol(cfg)
    water_stoich = water_per_extent(cfg)
    k = math.log(2.)/cfg.half_time_days*cfg.rate_multiplier
    flowing = mode != "stationary" and cfg.water_flow_kg_day > 0
    mixing = cfg.flow_mixing_factor if flowing else 1.
    film = 1/(1+cfg.film_damkohler/mixing)

    def rhs(time, y):
        extent, nh2, water, salt = y[:4]
        dissolved, _, csat, molality = _partition(nh2, water, salt, cfg)
        factors = feedback_factors(extent/nmax if nmax else 0., cfg)
        inhibition = 1/(1+cfg.product_inhibition_coefficient*dissolved/(water*csat))
        reaction = k*max(nmax-extent, 0.)*factors["area_factor"]*factors["passivation_factor"]*factors["pore_access_factor"]*film*inhibition
        q = cfg.water_flow_kg_day if flowing else 0.
        inlet_water = q if mode == "once_through" else 0.
        outlet_h2 = inlet_water*dissolved/water  # brine-only withdrawal; retain free gas
        outlet_salt = inlet_water*molality
        inlet_salt = inlet_water*cfg.salt_molality
        csat_separator = float(solubility_mol_kg(cfg.temperature_k, cfg.separator_bar, molality))
        separator_gas = inlet_water*max(dissolved/water-csat_separator, 0.)
        loop_dp = cfg.reference_loop_dp_bar/float(factors["permeability_ratio"]) if flowing else 0.
        feed_dp = cfg.pressure_bar-cfg.feed_pressure_bar if mode == "once_through" else 0.
        # Equal WATER inflow/outflow, with their respective salt concentration.
        volume_flow = q*(1+molality*M_SALT)/cfg.brine_density_kg_m3
        feed_volume_flow = inlet_water*(1+cfg.salt_molality*M_SALT)/cfg.brine_density_kg_m3
        pumping = (loop_dp*volume_flow+feed_dp*feed_volume_flow)*1e5/(3.6e6*cfg.pump_efficiency)
        return np.array([reaction, reaction-outlet_h2, -reaction*water_stoich,
                         inlet_salt-outlet_salt, outlet_h2, inlet_water, outlet_salt,
                         separator_gas, pumping, inlet_water, inlet_salt,
                         q if mode == "recirculating" else 0.])
    return rhs


def simulate(mode="stationary", cfg=Config(), times_days=None):
    """Return (history, summary). Endpoint recovery is potential, not summed in time.

    Once-through outlet parcels each receive their own instantaneous 10-bar
    flash. Dilute early effluent is never pooled into an artificial later flash.
    Closed recirculation has no gas-stripping unit and no water renewal.
    """
    cfg.validate()
    if mode not in MODES:
        raise ValueError(f"Unknown mode: {mode}")
    if times_days is None:
        times_days = np.unique(np.r_[np.linspace(0, cfg.duration_days, 366),
                                     [t for t in (0., 3., 30., 90., 365.) if t <= cfg.duration_days]])
    times = np.asarray(times_days, float)
    if times.ndim != 1 or len(times) < 2 or times[0] != 0 or np.any(np.diff(times) <= 0) or times[-1] > cfg.duration_days:
        raise ValueError("Output times must increase from zero within the duration")
    y0 = np.zeros(len(STATE_NAMES)); y0[2:4] = cfg.water_kg, cfg.water_kg*cfg.salt_molality
    solution = solve_ivp(_rhs(mode,cfg), (0,times[-1]), y0, t_eval=times,
                         rtol=cfg.relative_tolerance, atol=cfg.absolute_tolerance,
                         max_step=cfg.max_step_days, method="DOP853")
    if not solution.success:
        raise RuntimeError(solution.message)
    h = dict(zip(STATE_NAMES, solution.y)); h["time_days"] = times
    nmax = capacity_mol(cfg)
    f = h["extent_mol"]/nmax if nmax else np.zeros_like(times)
    h.update(feedback_factors(f,cfg))
    parts = [_partition(n,w,s,cfg) for n,w,s in zip(h["reactor_h2_mol"],h["water_kg"],h["salt_mol"])]
    dissolved, free, csat, molality = np.array(parts).T
    separator_capacity = h["water_kg"]*solubility_mol_kg(cfg.temperature_k,cfg.separator_bar,molality)
    endpoint_gas = np.maximum(h["reactor_h2_mol"]-separator_capacity,0.)
    h.update({"generated_h2_kg": h["extent_mol"]*M_H2,
              "retained_h2_kg": h["reactor_h2_mol"]*M_H2,
              "dissolved_h2_kg": dissolved*M_H2, "free_h2_kg": free*M_H2,
              "exported_h2_kg": h["export_h2_mol"]*M_H2,
              "outlet_collected_h2_kg": h["outlet_separator_gas_mol"]*M_H2*cfg.collection_efficiency,
              "outlet_uncollected_gas_h2_kg": h["outlet_separator_gas_mol"]*M_H2*(1-cfg.collection_efficiency),
              "outlet_residual_dissolved_h2_kg": (h["export_h2_mol"]-h["outlet_separator_gas_mol"])*M_H2,
              "endpoint_collected_h2_kg": endpoint_gas*M_H2*cfg.collection_efficiency,
              "endpoint_residual_h2_kg": (h["reactor_h2_mol"]-endpoint_gas)*M_H2,
              "endpoint_uncollected_gas_h2_kg": endpoint_gas*M_H2*(1-cfg.collection_efficiency),
              "salt_molality": molality,
              "water_consumed_kg": h["extent_mol"]*water_per_extent(cfg),
              "remaining_reactive_fe_kg": np.maximum(nmax-h["extent_mol"],0.)*3*M_FE,
              "flow_head_bar": np.where(mode!="stationary" and cfg.water_flow_kg_day>0,cfg.reference_loop_dp_bar/h["permeability_ratio"],0.)})
    h["solid_mass_gain_kg"] = h["water_consumed_kg"]-h["generated_h2_kg"]
    h["total_deliverable_h2_kg"] = h["outlet_collected_h2_kg"]+h["endpoint_collected_h2_kg"]
    h["gross_electricity_kwh"] = h["total_deliverable_h2_kg"]*cfg.lhv_kwh_kg*cfg.fuel_cell_efficiency_lhv
    h["outlet_delivered_electricity_kwh"] = h["outlet_collected_h2_kg"]*cfg.lhv_kwh_kg*cfg.fuel_cell_efficiency_lhv
    brine_initial = cfg.water_kg*(1+cfg.salt_molality*M_SALT)
    initial_pump = (cfg.pressure_bar-cfg.feed_pressure_bar)*1e5*brine_initial/cfg.brine_density_kg_m3/(3.6e6*cfg.pump_efficiency)
    initial_heat = max(cfg.temperature_k-cfg.feed_temperature_k,0)*(brine_initial*cfg.brine_cp_kj_kg_k+cfg.rock_kg*cfg.rock_cp_kj_kg_k)/3600/cfg.heater_efficiency*(1-cfg.heat_recovery_fraction)
    h["new_feed_heater_kwh"] = (h["feed_water_kg"]+h["feed_salt_mol"]*M_SALT)*cfg.brine_cp_kj_kg_k*max(cfg.temperature_k-cfg.feed_temperature_k,0)/3600/cfg.heater_efficiency*(1-cfg.heat_recovery_fraction)
    h["standalone_engineering_balance_kwh"] = h["gross_electricity_kwh"]-h["flow_pump_kwh"]-h["new_feed_heater_kwh"]-initial_pump-initial_heat
    h["h2_balance_error_kg"] = (h["extent_mol"]-h["reactor_h2_mol"]-h["export_h2_mol"])*M_H2
    h["water_balance_error_kg"] = h["water_kg"]+h["water_consumed_kg"]+h["export_water_kg"]-cfg.water_kg-h["feed_water_kg"]
    h["salt_balance_error_kg"] = (h["salt_mol"]+h["export_salt_mol"]-cfg.water_kg*cfg.salt_molality-h["feed_salt_mol"])*M_SALT
    h["hydrogen_recovery_balance_error_kg"] = h["generated_h2_kg"]-(h["total_deliverable_h2_kg"]+h["outlet_uncollected_gas_h2_kg"]+h["outlet_residual_dissolved_h2_kg"]+h["endpoint_residual_h2_kg"]+h["endpoint_uncollected_gas_h2_kg"])
    summary = {"mode":mode,"config":asdict(cfg),"hydrogen_capacity_kg":nmax*M_H2,
               "scope":"Uncalibrated reduced feedback scenario, not fracture mechanics or native DARTS",
               "final":{key:float(value[-1]) for key,value in h.items()},
               "max_abs_balance_errors":{key:float(np.max(np.abs(value))) for key,value in h.items() if "balance_error" in key},
               "initial_pump_kwh":initial_pump,"initial_heater_kwh":initial_heat,
               "head_diagnostic_limit_exceeded":bool(np.max(h["flow_head_bar"])>cfg.diagnostic_head_limit_bar),
               "function_evaluations":int(solution.nfev)}
    summary["sampled_first_time_to_1kg_days"] = {}
    for key in ("generated_h2_kg","retained_h2_kg","total_deliverable_h2_kg"):
        indices=np.flatnonzero(h[key]>=1.)
        if not len(indices):
            value=None
        else:
            i=int(indices[0]); value=float(np.interp(1.,h[key][i-1:i+1],times[i-1:i+1]))
        summary["sampled_first_time_to_1kg_days"][key]=value
    return h, summary
