"""
辐射室计算模块  –  Radiant section heat transfer
==================================================
Implements the Lobo-Evans single-zone model for the radiant (firebox) section
of a process fired heater.

The two governing equations are:

  (1) Heat balance:
        Q_radiant = Q_released - W_fg * h_fg(T_BWT) - Q_losses

  (2) Radiation exchange:
        Q_radiant = F * sigma * A_cp * (T_BWT^4 - T_tube^4)

These are solved simultaneously for the bridgewall temperature T_BWT and the
radiant heat absorption Q_radiant.

Tube skin temperature is estimated from the absorbed heat flux and the
combined film + wall resistance.

References:
    Lobo W.E. & Evans J.E. (1939) Trans. AIChE 35, 743.
    API STD 560, 4th ed.
    API RP 530 (tube design temperature / allowable stress)
    Wimpress R.N. (1963) Hydrocarbon Processing 42(10).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

from scipy.optimize import brentq

from .combustion import flue_gas_enthalpy_approx, combustion_calc, FUEL_DB
from .geometry import FireboxGeometry
from .properties import tube_metal_k, inner_film_htc

# Stefan-Boltzmann constant  [W / (m²·K⁴)]
SIGMA = 5.67e-8


@dataclass
class RadiantResult:
    """Results of the radiant-section calculation."""

    # ---- temperatures [°C] ------------------------------------------------
    T_BWT_C: float          # Bridgewall temperature (flue gas leaving radiant)
    T_AFT_C: float          # Adiabatic flame temperature
    T_tube_avg_C: float     # Average tube outer-wall temperature
    T_tube_max_C: float     # Estimated peak tube skin temperature
    T_process_out_C: float  # Process fluid outlet temperature from radiant

    # ---- heat duties [kW] -------------------------------------------------
    Q_released_kW: float    # Total heat released by combustion
    Q_radiant_kW: float     # Net heat absorbed by radiant tubes
    Q_losses_kW: float      # Heat losses (refractory conduction, openings)

    # ---- heat flux [kW / m²] ----------------------------------------------
    q_avg_kW_m2: float      # Average heat flux on tube surface
    q_max_kW_m2: float      # Peak heat flux (circumferential peak factor)

    # ---- exchange factor and geometry ------------------------------------
    F_exchange: float
    A_cp_m2: float

    # ---- fuel and flue gas -----------------------------------------------
    W_fuel_kg_s: float      # Fuel mass flow [kg/s]
    W_flue_kg_s: float      # Flue-gas mass flow [kg/s]

    def to_dict(self) -> dict:
        return {k: round(v, 4) if isinstance(v, float) else v
                for k, v in self.__dict__.items()}


def calc_radiant(
    # ---- fuel / combustion ------------------------------------------------
    fuel_type: str = "natural_gas",
    excess_air_frac: float = 0.15,
    W_fuel_kg_s: float = 1.0,
    T_air_K: float = 288.15,
    # ---- firebox geometry -------------------------------------------------
    geom: Optional[FireboxGeometry] = None,
    # ---- process fluid (radiant tubes) ------------------------------------
    process_fluid: str = "crude_oil",
    W_process_kg_s: float = 50.0,
    T_process_in_C: float = 200.0,    # inlet to radiant (from convection)
    T_process_out_C: Optional[float] = None,  # if None, calculated
    # ---- tube metallurgy --------------------------------------------------
    tube_material: str = "Cr5Mo",
    # ---- losses & factors -------------------------------------------------
    loss_fraction: float = 0.02,       # fraction of Q_released lost (2 %)
    peak_flux_factor: float = 1.8,     # circumferential peak / average
    T_ref_K: float = 288.15,
) -> RadiantResult:
    """
    Calculate the radiant-section heat transfer.

    Parameters
    ----------
    fuel_type : str         Key from FUEL_DB.
    excess_air_frac : float Fractional excess air (0.15 = 15 %).
    W_fuel_kg_s : float     Fuel mass-flow rate [kg/s].
    T_air_K : float         Combustion air temperature [K].
    geom : FireboxGeometry  Firebox / tube-bundle geometry object.
    process_fluid : str     Key from PROCESS_FLUID_PRESETS.
    W_process_kg_s : float  Process fluid mass flow [kg/s].
    T_process_in_C : float  Process fluid inlet temperature to radiant [°C].
    T_process_out_C : float Process fluid outlet temperature [°C]; if None,
                            calculated from Q_radiant and Cp.
    tube_material : str     Key from TUBE_MATERIALS.
    loss_fraction : float   Fraction of Q_released lost through refractory.
    peak_flux_factor : float  Ratio of peak to average heat flux.
    T_ref_K : float         Reference temperature for enthalpy [K].

    Returns
    -------
    RadiantResult
    """
    if geom is None:
        geom = FireboxGeometry()

    # ------------------------------------------------------------------
    # Combustion data
    # ------------------------------------------------------------------
    comb = combustion_calc(fuel_type, excess_air_frac)
    W_fg = W_fuel_kg_s * comb.flue_gas_per_fuel      # kg/s flue gas
    Q_released = W_fuel_kg_s * comb.LHV              # kW

    # Adiabatic flame temperature
    from .combustion import adiabatic_flame_temp
    T_AFT_K = adiabatic_flame_temp(fuel_type, excess_air_frac, T_air_K, T_ref_K)

    # Heat losses
    Q_losses = loss_fraction * Q_released

    # ------------------------------------------------------------------
    # Mean tube outer-wall temperature (first estimate: assume process
    # outlet + film ΔT; refined after Q_radiant is known)
    # ------------------------------------------------------------------
    from .properties import process_fluid_cp as pf_cp
    cp_proc = pf_cp(process_fluid)  # kJ/(kg·K)

    # Initial guess for T_process_out (if not given)
    # Will be updated after Q_radiant is found
    T_proc_in_K = T_process_in_C + 273.15

    # Approximate tube wall temperature for first iteration
    # Use an estimate: T_tube ≈ T_process_out + 50 K (film + wall ΔT)
    T_tube_est_K = T_proc_in_K + 100.0 + 50.0  # rough first guess

    # ------------------------------------------------------------------
    # Solve heat balance + radiation equation simultaneously for T_BWT
    # ------------------------------------------------------------------
    F = geom.F_exchange
    A_cp = geom.A_cp
    T_tube_K = T_tube_est_K  # will iterate if needed

    def heat_balance_residual(T_BWT_K: float) -> float:
        """
        residual = (Q_released - Q_losses - Q_fg_at_BWT) - Q_radiation
        At convergence this should be zero.
        """
        # Heat leaving with flue gas at T_BWT
        h_fg = flue_gas_enthalpy_approx(T_BWT_K, T_ref_K, excess_air_frac)
        Q_fg = W_fg * h_fg  # kW

        # Heat absorbed from heat-balance perspective
        Q_rad_hb = Q_released - Q_losses - Q_fg
        if Q_rad_hb < 0:
            return -1e6  # signal that T_BWT is too high

        # Radiation equation: Q_rad_rad
        Q_rad_rad = F * SIGMA * A_cp * (T_BWT_K ** 4 - T_tube_K ** 4) * 1e-3  # kW

        return Q_rad_hb - Q_rad_rad

    # Bracket T_BWT between T_tube_K + 10 and T_AFT_K - 10
    T_lo = T_tube_K + 100.0
    T_hi = min(T_AFT_K - 10.0, 1800.0)

    # Clamp search range
    T_lo = max(T_lo, T_ref_K + 50.0)
    T_hi = max(T_hi, T_lo + 200.0)

    # Two-pass iteration: solve for T_BWT, update T_tube, re-solve
    for _iteration in range(6):
        try:
            T_BWT_K = brentq(
                heat_balance_residual, T_lo, T_hi, xtol=0.1, maxiter=200
            )
        except ValueError:
            # If brentq fails (e.g. monotone residual), fall back to a midpoint
            T_BWT_K = (T_lo + T_hi) / 2.0

        # Compute Q_radiant from heat balance
        h_fg_bwt = flue_gas_enthalpy_approx(T_BWT_K, T_ref_K, excess_air_frac)
        Q_radiant = Q_released - Q_losses - W_fg * h_fg_bwt
        Q_radiant = max(0.0, Q_radiant)

        # Update process outlet temperature
        if T_process_out_C is None:
            dT_proc = Q_radiant / max(W_process_kg_s * cp_proc, 0.001)
            T_proc_out_K = T_proc_in_K + dT_proc
        else:
            T_proc_out_K = T_process_out_C + 273.15

        # Average tube outer-wall temperature
        # T_tube ≈ T_process_avg + ΔT_film + ΔT_wall
        T_proc_avg_K = (T_proc_in_K + T_proc_out_K) / 2.0

        # Average heat flux on tube surface
        A_tube_outer = geom.A_tube_outer
        if A_tube_outer > 0 and Q_radiant > 0:
            q_avg = Q_radiant * 1000.0 / A_tube_outer  # W/m²
        else:
            q_avg = 0.0

        # Inner film ΔT
        A_c_inner = geom.n_tubes_radiant * math.pi * geom.tube_id ** 2 / 4.0
        G = W_process_kg_s / max(A_c_inner, 1e-6)  # kg/(m²·s)
        h_i = inner_film_htc(process_fluid, G, geom.tube_id, T_proc_avg_K)
        # Flux referred to outer surface: q_o ≈ q_avg * (d_o / d_i)
        q_outer = q_avg * (geom.tube_od / geom.tube_id)
        dT_film = q_outer / max(h_i, 1.0)

        # Tube-wall ΔT
        k_m = tube_metal_k(tube_material, T_proc_avg_K + 50.0)
        r_o = geom.tube_od / 2.0
        r_i = geom.tube_id / 2.0
        dT_wall = q_outer * r_o * math.log(r_o / r_i) / k_m

        T_tube_avg_K = T_proc_avg_K + dT_film + dT_wall
        T_tube_K = T_tube_avg_K  # update for next iteration

    # ------------------------------------------------------------------
    # Final quantities
    # ------------------------------------------------------------------
    T_BWT_C = T_BWT_K - 273.15
    T_AFT_C = T_AFT_K - 273.15
    T_tube_avg_C = T_tube_avg_K - 273.15

    q_avg_kW_m2 = q_avg * 1e-3  # kW/m²
    q_max_kW_m2 = q_avg_kW_m2 * peak_flux_factor

    # Peak tube skin temperature (at maximum flux location)
    T_tube_max_K = (
        T_proc_out_K  # worst case: outlet process temp
        + (q_max_kW_m2 * 1000.0 / max(h_i, 1.0))   # film ΔT at peak flux
        + q_max_kW_m2 * 1000.0 * r_o * math.log(r_o / r_i) / k_m  # wall ΔT
    )
    T_tube_max_C = T_tube_max_K - 273.15

    T_proc_out_C_final = T_proc_out_K - 273.15

    return RadiantResult(
        T_BWT_C=round(T_BWT_C, 1),
        T_AFT_C=round(T_AFT_C, 1),
        T_tube_avg_C=round(T_tube_avg_C, 1),
        T_tube_max_C=round(T_tube_max_C, 1),
        T_process_out_C=round(T_proc_out_C_final, 1),
        Q_released_kW=round(Q_released, 1),
        Q_radiant_kW=round(Q_radiant, 1),
        Q_losses_kW=round(Q_losses, 1),
        q_avg_kW_m2=round(q_avg_kW_m2, 2),
        q_max_kW_m2=round(q_max_kW_m2, 2),
        F_exchange=round(F, 4),
        A_cp_m2=round(A_cp, 2),
        W_fuel_kg_s=round(W_fuel_kg_s, 4),
        W_flue_kg_s=round(W_fg, 4),
    )
