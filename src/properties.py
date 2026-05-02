"""
热物性模块  –  Thermophysical properties
==========================================
Flue-gas transport properties, tube-metal thermal conductivity,
and simple process-fluid properties.

All temperatures in Kelvin unless noted.

References:
    Perry's Chemical Engineers' Handbook, 8th ed.
    Incropera & DeWitt, "Fundamentals of Heat and Mass Transfer"
"""
from __future__ import annotations

import math


# ---------------------------------------------------------------------------
# Flue-gas transport properties
# ---------------------------------------------------------------------------

def flue_gas_viscosity(T_K: float) -> float:
    """
    Dynamic viscosity of flue gas [Pa·s] at temperature T_K.

    Sutherland-type polynomial fit for typical combustion products.
    Valid 300–2000 K.
    """
    T_C = T_K - 273.15
    mu = (1.46e-5 + 4.1e-8 * T_C) * (1.0 + 0.0005 * (T_C / 1000.0) ** 2)
    return max(1e-6, mu)


def flue_gas_thermal_conductivity(T_K: float) -> float:
    """
    Thermal conductivity of flue gas [W / (m·K)] at temperature T_K.

    Polynomial fit valid 300–2000 K.
    """
    T_C = T_K - 273.15
    lam = 0.02442 + 7.2e-5 * T_C - 2.0e-8 * T_C ** 2
    return max(0.01, lam)


def flue_gas_density(T_K: float, P_kPa: float = 101.325, MW: float = 28.5) -> float:
    """
    Density of flue gas [kg / m³].

    Uses ideal-gas law.  MW is the mean molecular weight of the flue-gas
    mixture (kg / kmol), default 28.5 for typical natural-gas combustion
    products with ~15 % excess air.
    """
    R = 8.314  # kJ / (kmol·K)
    rho = (P_kPa * MW) / (R * T_K)
    return rho


def flue_gas_prandtl(T_K: float, excess_air_frac: float = 0.15) -> float:
    """
    Prandtl number of flue gas [-] at temperature T_K.

    Pr = mu * cp / lambda  (using SI units throughout)
    """
    from .combustion import flue_gas_cp
    mu = flue_gas_viscosity(T_K)             # Pa·s
    lam = flue_gas_thermal_conductivity(T_K)  # W/(m·K)
    cp = flue_gas_cp(T_K, excess_air_frac) * 1000.0  # J/(kg·K)
    return mu * cp / lam


# ---------------------------------------------------------------------------
# Tube-metal thermal conductivity
# ---------------------------------------------------------------------------

TUBE_MATERIALS = {
    "carbon_steel": {
        "name": "碳钢 / Carbon Steel",
        "k_300K": 50.0,   # W/(m·K) at ~300 K
        "dk_dT": -0.030,  # W/(m·K) per K (linear fit)
    },
    "Cr5Mo": {
        "name": "5Cr½Mo (P5)",
        "k_300K": 38.0,
        "dk_dT": -0.015,
    },
    "Cr9Mo": {
        "name": "9Cr1Mo (P9/T9)",
        "k_300K": 28.0,
        "dk_dT": -0.010,
    },
    "SS304": {
        "name": "304不锈钢 / 304 Stainless Steel",
        "k_300K": 14.9,
        "dk_dT": 0.012,
    },
    "SS316": {
        "name": "316不锈钢 / 316 Stainless Steel",
        "k_300K": 13.4,
        "dk_dT": 0.012,
    },
    "Incoloy800": {
        "name": "Incoloy 800H",
        "k_300K": 11.5,
        "dk_dT": 0.018,
    },
}


def tube_metal_k(material: str, T_K: float = 700.0) -> float:
    """
    Thermal conductivity of tube material [W / (m·K)].

    Parameters
    ----------
    material : str
        Key from TUBE_MATERIALS.
    T_K : float
        Mean tube-wall temperature in K.
    """
    mat = TUBE_MATERIALS.get(material, TUBE_MATERIALS["carbon_steel"])
    k = mat["k_300K"] + mat["dk_dT"] * (T_K - 300.0)
    return max(5.0, k)


# ---------------------------------------------------------------------------
# Air properties (used in combustion pre-heat)
# ---------------------------------------------------------------------------

def air_density(T_K: float, P_kPa: float = 101.325) -> float:
    """Density of dry air [kg/m³]."""
    return (P_kPa * 28.97) / (8.314 * T_K)


def air_viscosity(T_K: float) -> float:
    """Dynamic viscosity of dry air [Pa·s]."""
    return 1.716e-5 * (T_K / 273.15) ** 1.5 * (273.15 + 110.4) / (T_K + 110.4)


def air_thermal_conductivity(T_K: float) -> float:
    """Thermal conductivity of dry air [W/(m·K)]."""
    return 0.0241 * (T_K / 273.15) ** 0.82


# ---------------------------------------------------------------------------
# Simple process-fluid properties (liquid hydrocarbon approximations)
# ---------------------------------------------------------------------------

PROCESS_FLUID_PRESETS = {
    "crude_oil": {
        "name": "原油 / Crude Oil",
        "rho_20C": 850.0,   # kg/m³ at 20 °C
        "cp_avg": 2.2,       # kJ/(kg·K)  average
        "mu_50C": 5e-3,      # Pa·s  at 50 °C
        "k_liq": 0.135,      # W/(m·K)
    },
    "naphtha": {
        "name": "石脑油 / Naphtha",
        "rho_20C": 720.0,
        "cp_avg": 2.4,
        "mu_50C": 0.5e-3,
        "k_liq": 0.120,
    },
    "gas_oil": {
        "name": "蜡油 / Gas Oil",
        "rho_20C": 870.0,
        "cp_avg": 2.1,
        "mu_50C": 3e-3,
        "k_liq": 0.130,
    },
    "vacuum_residue": {
        "name": "减压渣油 / Vacuum Residue",
        "rho_20C": 990.0,
        "cp_avg": 2.0,
        "mu_50C": 500e-3,
        "k_liq": 0.145,
    },
    "water_steam": {
        "name": "水/蒸汽 / Water-Steam",
        "rho_20C": 998.0,
        "cp_avg": 4.18,
        "mu_50C": 0.55e-3,
        "k_liq": 0.640,
    },
}


def process_fluid_cp(fluid: str) -> float:
    """Mean specific heat [kJ/(kg·K)] of process fluid."""
    return PROCESS_FLUID_PRESETS.get(fluid, PROCESS_FLUID_PRESETS["crude_oil"])["cp_avg"]


def process_fluid_k(fluid: str) -> float:
    """Thermal conductivity [W/(m·K)] of process fluid."""
    return PROCESS_FLUID_PRESETS.get(fluid, PROCESS_FLUID_PRESETS["crude_oil"])["k_liq"]


def process_fluid_mu(fluid: str) -> float:
    """Dynamic viscosity [Pa·s] of process fluid at ~50 °C."""
    return PROCESS_FLUID_PRESETS.get(fluid, PROCESS_FLUID_PRESETS["crude_oil"])["mu_50C"]


def process_fluid_prandtl(fluid: str) -> float:
    """Approximate Prandtl number of process fluid."""
    cp = process_fluid_cp(fluid) * 1000.0  # J/(kg·K)
    mu = process_fluid_mu(fluid)            # Pa·s
    lam = process_fluid_k(fluid)            # W/(m·K)
    return mu * cp / lam


def process_fluid_mu_at_T(fluid: str, T_K: float) -> float:
    """
    Temperature-corrected dynamic viscosity [Pa·s] of process fluid.

    Uses the Andrade equation (Arrhenius-type):
        mu(T) = mu_ref * exp(B * (1/T - 1/T_ref))

    where T_ref = 323.15 K (50 °C) and B is derived from the fact that
    petroleum viscosity roughly halves every 20–25 °C.
    """
    mu_ref = process_fluid_mu(fluid)       # Pa·s at 50 °C (323.15 K)
    T_ref = 323.15
    # B calibrated so viscosity drops ~10× over 150 °C (typical petroleum)
    B = 2300.0   # K  (Andrade constant, approximate for petroleum fractions)
    if fluid == "water_steam":
        # Water viscosity uses a different correlation
        B = 1700.0
    mu_T = mu_ref * math.exp(B * (1.0 / T_K - 1.0 / T_ref))
    # Clamp: keep within a physically reasonable range
    return max(1.0e-5, min(mu_ref * 5.0, mu_T))


def inner_film_htc(
    fluid: str,
    mass_flux: float,   # kg/(m²·s)  through tube inner area
    d_i: float,         # m  tube inner diameter
    T_avg_K: float = 500.0,
) -> float:
    """
    Process-side (inner) heat transfer coefficient [W / (m²·K)].

    Uses the Dittus-Boelter correlation:
      Nu = 0.023 * Re^0.8 * Pr^0.4  (heating)

    Parameters
    ----------
    fluid : str
        Key from PROCESS_FLUID_PRESETS.
    mass_flux : float
        Mass velocity G = W / A_c  [kg/(m²·s)].
    d_i : float
        Tube inner diameter [m].
    T_avg_K : float
        Average process fluid temperature [K].
    """
    mu = process_fluid_mu_at_T(fluid, T_avg_K)   # temperature-corrected
    cp = process_fluid_cp(fluid) * 1000.0          # J/(kg·K)
    lam = process_fluid_k(fluid)                   # W/(m·K)

    Re = mass_flux * d_i / mu
    Pr = mu * cp / lam

    if Re < 2300:
        # Laminar – use Sieder-Tate simplified
        Nu = 3.66 + (0.0668 * (d_i / 10.0) * Re * Pr) / (
            1.0 + 0.04 * ((d_i / 10.0) * Re * Pr) ** (2.0 / 3.0)
        )
    else:
        # Turbulent – Dittus-Boelter
        Nu = 0.023 * Re ** 0.8 * Pr ** 0.4

    h = Nu * lam / d_i
    return max(100.0, h)
