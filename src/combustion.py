"""
燃烧计算模块  –  Combustion calculations
==========================================
Handles fuel properties, stoichiometry, flue-gas composition and
heat-release for the fired heater simulation.

Reference:
    API STD 560, 4th ed. (Fired Heaters for General Refinery Service)
    Ganapathy V., "Industrial Boilers and Heat Recovery Steam Generators"
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict

# ---------------------------------------------------------------------------
# Fuel database
# ---------------------------------------------------------------------------
# Composition keys for gaseous fuels: component symbol → mole fraction
# Composition keys for liquid fuels:  element symbol   → mass fraction

FUEL_DB: Dict[str, dict] = {
    "natural_gas": {
        "name": "天然气 / Natural Gas",
        "phase": "gas",
        "LHV": 47_141.0,       # kJ / kg  (lower heating value)
        "HHV": 52_225.0,       # kJ / kg
        "stoic_air_mass": 17.2,  # kg air / kg fuel  (stoichiometric)
        "MW": 17.4,             # kg / kmol
        "composition_mol": {    # typical pipeline natural gas
            "CH4": 0.900,
            "C2H6": 0.050,
            "C3H8": 0.020,
            "N2": 0.020,
            "CO2": 0.010,
        },
    },
    "fuel_oil": {
        "name": "重油 / Fuel Oil",
        "phase": "liquid",
        "LHV": 41_800.0,
        "HHV": 44_600.0,
        "stoic_air_mass": 13.8,
        "MW": 250.0,
        "composition_mass": {   # typical No. 6 fuel oil (mass fractions)
            "C": 0.86,
            "H": 0.12,
            "S": 0.015,
            "N": 0.005,
        },
    },
    "refinery_gas": {
        "name": "炼厂气 / Refinery Gas",
        "phase": "gas",
        "LHV": 44_000.0,
        "HHV": 48_000.0,
        "stoic_air_mass": 15.5,
        "MW": 30.0,
        "composition_mol": {
            "CH4": 0.60,
            "C2H6": 0.15,
            "C3H8": 0.10,
            "H2": 0.10,
            "N2": 0.05,
        },
    },
    "lpg": {
        "name": "液化石油气 / LPG",
        "phase": "gas",
        "LHV": 46_100.0,
        "HHV": 50_000.0,
        "stoic_air_mass": 15.7,
        "MW": 44.1,
        "composition_mol": {
            "C3H8": 0.60,
            "C4H10": 0.40,
        },
    },
}

# Dry air composition (mole fractions)
AIR_COMPOSITION_MOL = {"N2": 0.7809, "O2": 0.2095, "Ar": 0.0096}
MW_AIR = 28.97   # kg / kmol
CP_AIR = 1.005   # kJ / (kg·K)  – treated as constant


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

@dataclass
class CombustionResult:
    """Results from a combustion calculation."""
    fuel_type: str
    LHV: float           # kJ / kg fuel
    stoic_air: float     # kg air / kg fuel  (theoretical)
    actual_air: float    # kg air / kg fuel  (with excess)
    excess_air_frac: float   # 0.20 = 20 %
    flue_gas_per_fuel: float  # kg flue-gas / kg fuel
    flue_comp_mass: Dict[str, float]  # mass fractions of flue gas
    T_ref_K: float = 288.15  # reference temperature (15 °C)


def get_fuel(fuel_type: str) -> dict:
    """Return fuel data dict; raise KeyError for unknown types."""
    if fuel_type not in FUEL_DB:
        raise KeyError(
            f"Unknown fuel type '{fuel_type}'. "
            f"Available: {list(FUEL_DB.keys())}"
        )
    return FUEL_DB[fuel_type]


def combustion_calc(
    fuel_type: str = "natural_gas",
    excess_air_frac: float = 0.15,
) -> CombustionResult:
    """
    Perform stoichiometric combustion calculation.

    Parameters
    ----------
    fuel_type : str
        Key from FUEL_DB.
    excess_air_frac : float
        Fractional excess air, e.g. 0.15 for 15 %.

    Returns
    -------
    CombustionResult
    """
    fuel = get_fuel(fuel_type)
    stoic = fuel["stoic_air_mass"]
    actual = stoic * (1.0 + excess_air_frac)
    flue_per_fuel = 1.0 + actual  # kg flue / kg fuel

    # ------------------------------------------------------------------
    # Approximate flue-gas mass composition
    # We split the flue gas into CO2, H2O, N2, O2 (and SO2 for oil)
    # ------------------------------------------------------------------
    if fuel["phase"] == "gas":
        comp = _flue_comp_from_gas(fuel, excess_air_frac, actual)
    else:
        comp = _flue_comp_from_liquid(fuel, excess_air_frac, actual)

    return CombustionResult(
        fuel_type=fuel_type,
        LHV=fuel["LHV"],
        stoic_air=stoic,
        actual_air=actual,
        excess_air_frac=excess_air_frac,
        flue_gas_per_fuel=flue_per_fuel,
        flue_comp_mass=comp,
    )


# ---------------------------------------------------------------------------
# Flue-gas enthalpy and heat functions
# ---------------------------------------------------------------------------

def flue_gas_cp(T_K: float, excess_air_frac: float = 0.15) -> float:
    """
    Mean specific heat of combustion flue gas at temperature *T_K* [K].

    Polynomial fit valid 300–2000 K for typical natural-gas/oil combustion
    products with 10–30 % excess air.

    Returns
    -------
    cp : float  [kJ / (kg·K)]
    """
    T_C = T_K - 273.15
    # Ganapathy (2003) polynomial – dry + wet flue gas blend
    cp = 1.0047 + 1.668e-4 * T_C + 1.5e-8 * T_C ** 2
    # Slight correction for excess air level
    cp += 0.002 * (excess_air_frac - 0.15)
    return cp


def flue_gas_enthalpy(
    T_K: float,
    T_ref_K: float = 288.15,
    excess_air_frac: float = 0.15,
) -> float:
    """
    Specific enthalpy of flue gas [kJ/kg] relative to *T_ref_K*.

    Uses numerical integration of :func:`flue_gas_cp`.
    """
    from scipy.integrate import quad
    h, _ = quad(
        lambda T: flue_gas_cp(T, excess_air_frac),
        T_ref_K,
        T_K,
        limit=50,
    )
    return h


def flue_gas_enthalpy_approx(
    T_K: float,
    T_ref_K: float = 288.15,
    excess_air_frac: float = 0.15,
) -> float:
    """
    Fast closed-form enthalpy approximation (no numerical integration).

    H(T) = integral of cp dT from T_ref to T
         = a*(T-T_ref) + b/2*(T²-T_ref²) + c/3*(T³-T_ref³)
    (T in °C, result in kJ/kg)
    """
    T_C = T_K - 273.15
    T_C_ref = T_ref_K - 273.15
    a = 1.0047 + 0.002 * (excess_air_frac - 0.15)
    b = 1.668e-4
    c = 1.5e-8
    h = (
        a * (T_C - T_C_ref)
        + b / 2.0 * (T_C ** 2 - T_C_ref ** 2)
        + c / 3.0 * (T_C ** 3 - T_C_ref ** 3)
    )
    return h


def adiabatic_flame_temp(
    fuel_type: str = "natural_gas",
    excess_air_frac: float = 0.15,
    T_air_K: float = 288.15,
    T_ref_K: float = 288.15,
) -> float:
    """
    Estimate adiabatic flame temperature [K] via energy balance.

    Energy available = LHV + air sensible heat above reference
    That energy heats the flue gas from T_ref to T_adiab.
    """
    from scipy.optimize import brentq

    fuel = get_fuel(fuel_type)
    LHV = fuel["LHV"]
    stoic = fuel["stoic_air_mass"]
    actual_air = stoic * (1.0 + excess_air_frac)
    flue_per_fuel = 1.0 + actual_air

    q_air = actual_air * CP_AIR * (T_air_K - T_ref_K)
    q_total = LHV + q_air  # kJ / kg fuel

    def residual(T_K: float) -> float:
        h = flue_gas_enthalpy_approx(T_K, T_ref_K, excess_air_frac)
        return h * flue_per_fuel - q_total

    T_aft = brentq(residual, 800.0, 2800.0, xtol=0.5)
    return T_aft


def available_heat_fraction(
    fuel_type: str = "natural_gas",
    excess_air_frac: float = 0.15,
    T_stack_K: float = 450.0,
    T_air_K: float = 288.15,
    T_ref_K: float = 288.15,
) -> float:
    """
    Available heat fraction η_avail = (LHV – stack loss + air preheat) / LHV.

    Stack loss = specific enthalpy of flue gas at stack temperature.
    """
    fuel = get_fuel(fuel_type)
    LHV = fuel["LHV"]
    stoic = fuel["stoic_air_mass"]
    actual_air = stoic * (1.0 + excess_air_frac)
    flue_per_fuel = 1.0 + actual_air

    h_stack = flue_gas_enthalpy_approx(T_stack_K, T_ref_K, excess_air_frac)
    q_stack = h_stack * flue_per_fuel  # kJ / kg fuel

    q_air = actual_air * CP_AIR * (T_air_K - T_ref_K)

    eta = (LHV + q_air - q_stack) / LHV
    return max(0.0, min(1.0, eta))


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _flue_comp_from_gas(fuel: dict, excess_air_frac: float, actual_air: float) -> dict:
    """
    Approximate mass-fraction composition of flue gas for a gaseous fuel.
    Assumes complete combustion.
    """
    # Molecular weights
    MW = {"CH4": 16, "C2H6": 30, "C3H8": 44, "H2": 2,
          "N2": 28, "CO2": 44, "H2O": 18, "O2": 32}

    comp_mol = fuel.get("composition_mol", {})
    MW_fuel = sum(MW.get(k, 28) * v for k, v in comp_mol.items())
    MW_fuel = MW_fuel if MW_fuel > 0 else fuel.get("MW", 28)

    # Per kg of fuel (in mol):
    mol_fuel = 1000.0 / MW_fuel  # kmol / kg fuel   (× 1000 for kg→g conversion)

    # Products per kmol of each component
    products = {"CO2": 0.0, "H2O": 0.0, "N2_fuel": 0.0}
    for comp, mol_frac in comp_mol.items():
        mol_comp = mol_fuel * mol_frac  # kmol / kg fuel
        if comp == "CH4":
            products["CO2"] += mol_comp * 1
            products["H2O"] += mol_comp * 2
        elif comp == "C2H6":
            products["CO2"] += mol_comp * 2
            products["H2O"] += mol_comp * 3
        elif comp == "C3H8":
            products["CO2"] += mol_comp * 3
            products["H2O"] += mol_comp * 4
        elif comp == "H2":
            products["H2O"] += mol_comp * 1
        elif comp == "N2":
            products["N2_fuel"] += mol_comp
        elif comp == "CO2":
            products["CO2"] += mol_comp

    # Convert products to kg / kg fuel
    kg_CO2 = products["CO2"] * 44
    kg_H2O = products["H2O"] * 18
    kg_N2_fuel = products["N2_fuel"] * 28

    # Air components (kg / kg fuel)
    stoic = fuel["stoic_air_mass"]
    # air = 23.2 % O2 + 76.8 % N2 by mass
    kg_N2_air = actual_air * 0.768
    kg_O2_air = actual_air * 0.232

    # Excess O2 in flue gas
    kg_O2_excess = kg_O2_air - (actual_air - actual_air * excess_air_frac / (1 + excess_air_frac)) * 0.0
    # Simpler: excess O2 = (excess_air / (1+excess_air)) * total O2
    kg_O2_total = actual_air * 0.232
    kg_O2_consumed = stoic * 0.232   # stoichiometric O2
    kg_O2_excess = kg_O2_total - kg_O2_consumed

    # Total mass of flue gas per kg fuel
    total = kg_CO2 + kg_H2O + (kg_N2_fuel + kg_N2_air) + max(0.0, kg_O2_excess)

    if total <= 0:
        return {"CO2": 0.13, "H2O": 0.12, "N2": 0.72, "O2": 0.03}

    return {
        "CO2": kg_CO2 / total,
        "H2O": kg_H2O / total,
        "N2": (kg_N2_fuel + kg_N2_air) / total,
        "O2": max(0.0, kg_O2_excess) / total,
    }


def _flue_comp_from_liquid(fuel: dict, excess_air_frac: float, actual_air: float) -> dict:
    """
    Approximate mass-fraction composition of flue gas for a liquid fuel.
    """
    comp = fuel.get("composition_mass", {"C": 0.86, "H": 0.12, "S": 0.015})
    C = comp.get("C", 0.86)
    H = comp.get("H", 0.12)
    S = comp.get("S", 0.015)

    # Products per kg fuel (assuming complete combustion)
    kg_CO2 = C * 44.0 / 12.0        # C → CO2
    kg_H2O = H * 18.0 / 2.0         # H2 → H2O  (H is mass fraction of H atoms)
    kg_SO2 = S * 64.0 / 32.0        # S → SO2

    stoic = fuel["stoic_air_mass"]
    kg_O2_consumed = (
        C * 32.0 / 12.0             # C + O2 → CO2
        + H * 16.0 / 2.0            # H2 + ½O2 → H2O  (H here = mass fraction H)
        + S * 32.0 / 32.0           # S + O2 → SO2
    )
    kg_N2_air = actual_air * 0.768
    kg_O2_total = actual_air * 0.232
    kg_O2_excess = kg_O2_total - kg_O2_consumed

    total = kg_CO2 + kg_H2O + kg_SO2 + kg_N2_air + max(0.0, kg_O2_excess)

    if total <= 0:
        return {"CO2": 0.13, "H2O": 0.11, "N2": 0.72, "SO2": 0.01, "O2": 0.03}

    return {
        "CO2": kg_CO2 / total,
        "H2O": kg_H2O / total,
        "SO2": kg_SO2 / total,
        "N2": kg_N2_air / total,
        "O2": max(0.0, kg_O2_excess) / total,
    }
