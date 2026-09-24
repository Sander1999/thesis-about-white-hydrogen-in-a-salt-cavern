"""Finite-inventory, prescribed-pressure olivine/brine reaction scenarios.

This is a reduced redox model, not a geochemical equilibrium calculation.
Measured total Fe bounds capacity. Neither Fe(II), accessibility nor kinetics
is measured by XRF. All reaction scenarios are assumptions, including zero.
Units here are kg, mol, m3, bar, kelvin and days (DARTS uses kmol separately).
"""
from dataclasses import dataclass, asdict, replace
import math
import numpy as np

M_FE = 0.055845
M_H2 = 0.00201588
M_WATER = 0.01801528
M_SALT = 0.0584428
R = 8.314462618


@dataclass(frozen=True)
class Setup:
    rock_kg: float = 10_000.
    fe_mass_fraction: float = 0.059187
    ferrous_fraction: float = 1.
    accessible_fraction: float = 1.
    water_kg: float = 20_000.
    salt_molality: float = 4.5
    temperature_k: float = 323.15
    pressure_bar: float = 120.
    separator_bar: float = 10.
    brine_density_kg_m3: float = 1170.
    collection_efficiency: float = 0.90
    fuel_cell_efficiency_lhv: float = 0.50
    lhv_kwh_kg: float = 33.33
    feed_temperature_k: float = 293.15
    brine_cp_kj_kg_k: float = 3.3
    rock_cp_kj_kg_k: float = 0.85
    pump_efficiency: float = 0.70
    feed_pressure_bar: float = 1.
    warm_system_extra_dp_bar: float = 5.
    heater_efficiency: float = 0.95
    heat_recovery_fraction: float = 0.
    # Extra hydration water beyond the redox reaction is not measured.
    # The baseline uses only the balanced redox water: one mole per mole H2.
    water_mol_per_h2: float = 1.

    def validate(self):
        if any(not math.isfinite(v) for v in asdict(self).values()):
            raise ValueError('Setup values must be finite')
        for name in ('ferrous_fraction', 'accessible_fraction', 'collection_efficiency',
                     'fuel_cell_efficiency_lhv', 'heat_recovery_fraction'):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f'{name} must lie between zero and one')
        if self.rock_kg < 0 or not 0 <= self.fe_mass_fraction <= 1 or self.water_kg <= 0:
            raise ValueError('Invalid rock, Fe or water inventory')
        if self.water_mol_per_h2 < 1:
            raise ValueError('At least redox water must be consumed')
        if not 0 < self.pump_efficiency <= 1 or not 0 < self.heater_efficiency <= 1:
            raise ValueError('Invalid pump or heater efficiency')
        if min(self.brine_density_kg_m3, self.lhv_kwh_kg, self.brine_cp_kj_kg_k,
               self.rock_cp_kj_kg_k) <= 0:
            raise ValueError('Density, LHV and heat capacities must be positive')
        if not self.feed_pressure_bar <= self.separator_bar <= self.pressure_bar:
            raise ValueError('Pressure order must be feed <= separator <= reactor')
        if self.feed_pressure_bar <= 0 or self.feed_temperature_k <= 0 or self.warm_system_extra_dp_bar < 0:
            raise ValueError('Absolute feed pressure/temperature must be positive; pumping differential nonnegative')
        solubility_mol_kg(self.temperature_k, self.pressure_bar, self.salt_molality)
        solubility_mol_kg(self.temperature_k, self.separator_bar, self.salt_molality)
        return self


def solubility_mol_kg(temperature_k, pressure_bar, molality):
    """Chabab et al. (2020), empirical correlation, salt-free mole fraction.

    DOI:10.1016/j.ijhydene.2020.08.192. Validated domain only: 323.15–373.15K,
    10–230bar, 0–5mol NaCl/kg water. No fugacity multiplier is added.
    """
    t, p, m = np.broadcast_arrays(temperature_k, pressure_bar, molality)
    if np.any((t < 323.15-1e-8) | (t > 373.15) | (p < 10) | (p > 230) |
              (m < 0) | (m > 5)):
        raise ValueError('Outside Chabab fit: 323.15–373.15K, 10–230bar, 0–5molal')
    x0 = p*(3.338844e-7*t + 0.0363161/t - 0.00020734) - 2.1301815e-9*p*p
    x = x0*np.exp(0.018519*m*m - 0.30185103*m)
    return x/(M_WATER*(1-x))


def capacity_mol(cfg):
    """Magnetite pathway: 3FeO + H2O -> Fe3O4 + H2 (FeO is a proxy)."""
    return (cfg.rock_kg*cfg.fe_mass_fraction/M_FE * cfg.ferrous_fraction *
            cfg.accessible_fraction/3.)


def reaction_extent(times_days, half_time_days, cfg):
    """Analytic finite-inventory first order conversion; infinite half-time = zero."""
    times = np.asarray(times_days, float)
    if np.any(~np.isfinite(times)) or np.any(times < 0) or half_time_days <= 0 or math.isnan(half_time_days):
        raise ValueError('Nonnegative times and positive half-time required')
    iron_extent = capacity_mol(cfg)*(-np.expm1(-math.log(2.)*times/half_time_days))
    return np.minimum(iron_extent, cfg.water_kg/(M_WATER*cfg.water_mol_per_h2))


def evaluate(times_days, half_time_days, cfg=Setup()):
    """Batch endpoints, each independently recovered once, never cumulatively summed.

    Total pressure is imposed by an external pressure boundary; variable gas
    volume is diagnosed. No claim of rigid-vessel pressure prediction is made.
    Separator performs one equilibrium flash at 10bar. A dry-gas collection
    factor is applied afterwards; unrecovered gas and dissolved H2 stay in ledger.
    """
    cfg.validate()
    time = np.atleast_1d(np.asarray(times_days, float))
    produced = reaction_extent(time, half_time_days, cfg)
    water = cfg.water_kg - produced*M_WATER*cfg.water_mol_per_h2
    if np.any(water <= 0):
        raise ValueError('Water exhausted: aqueous partition is no longer defined')
    salt_mol = cfg.water_kg*cfg.salt_molality
    molality = salt_mol/water
    if np.any(molality > 5):
        raise ValueError('Reaction concentrated brine beyond solubility fit; reduce initial salt')
    dissolved = np.minimum(produced, water*solubility_mol_kg(
        cfg.temperature_k, cfg.pressure_bar, molality))
    free = produced-dissolved
    retained = np.minimum(produced, water*solubility_mol_kg(
        cfg.temperature_k, cfg.separator_bar, molality))
    flashed = produced-retained
    collected = cfg.collection_efficiency*flashed
    uncollected = (1-cfg.collection_efficiency)*flashed
    salt_kg = salt_mol*M_SALT
    brine_kg = cfg.water_kg+salt_kg
    brine_volume = brine_kg/cfg.brine_density_kg_m3
    # Engineering boundary costs, not a coupled thermal simulation.
    initial_pump = (cfg.pressure_bar-cfg.feed_pressure_bar)*1e5*brine_volume/(3.6e6*cfg.pump_efficiency)
    warm_pump = cfg.warm_system_extra_dp_bar*1e5*brine_volume/(3.6e6*cfg.pump_efficiency)
    heat = max(0., cfg.temperature_k-cfg.feed_temperature_k)*(
        brine_kg*cfg.brine_cp_kj_kg_k+cfg.rock_kg*cfg.rock_cp_kj_kg_k)/3600
    heat_electric = heat*(1-cfg.heat_recovery_fraction)/cfg.heater_efficiency
    gross = collected*M_H2*cfg.lhv_kwh_kg*cfg.fuel_cell_efficiency_lhv
    ones = np.ones(time.shape)
    return {
        'time_days': time, 'generated_h2_kg': produced*M_H2,
        'dissolved_reactor_h2_kg': dissolved*M_H2, 'free_reactor_h2_kg': free*M_H2,
        'retained_separator_h2_kg': retained*M_H2,
        'uncollected_separator_h2_kg': uncollected*M_H2,
        'collected_h2_kg': collected*M_H2,
        'remaining_reactive_fe_kg': (capacity_mol(cfg)-produced)*3*M_FE,
        'water_remaining_kg': water, 'salt_kg': ones*salt_kg, 'salt_molality': molality,
        'free_h2_volume_m3_ideal': free*R*cfg.temperature_k/(cfg.pressure_bar*1e5),
        'gross_electricity_kwh': gross,
        'average_gross_power_kw': np.divide(gross,time*24,out=np.zeros_like(gross),where=time>0),
        'initial_pump_electricity_kwh': ones*initial_pump,
        'sensible_heat_kwh_thermal': ones*heat,
        'heater_electricity_kwh': ones*heat_electric,
        'standalone_balance_kwh': gross-initial_pump-heat_electric,
        'warm_pressurized_balance_kwh': gross-warm_pump,
        'h2_ledger_error_kg': (produced-collected-uncollected-retained)*M_H2,
    }


SCENARIOS = {
    'zero': math.inf,
    'ten_year_assumption': 3650.,
    'exploratory_assumption': 365.,
    'accelerated_target': 3.,
}


def metadata(cfg=Setup()):
    cfg.validate()
    return {'setup': asdict(cfg), 'magnetite_h2_ceiling_kg': capacity_mol(cfg)*M_H2,
            'absolute_ferric_h2_ceiling_kg': 1.5*capacity_mol(cfg)*M_H2,
            'salinity_g_l': cfg.salt_molality*M_SALT/(1+cfg.salt_molality*M_SALT)*cfg.brine_density_kg_m3,
            'scenarios_half_time_days': {k: (None if math.isinf(v) else v) for k,v in SCENARIOS.items()},
            'kinetics_status': 'Uncalibrated scenarios, not measurements or probability bounds',
            'model_scope': 'Finite redox proxy, equilibrium partition, external engineering energy budget'}
