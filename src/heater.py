"""
加热炉主模型  –  Fired heater integration model
=================================================
Integrates the radiant and convection sections into a complete fired-heater
simulation, analogous to the FRNC-5PC calculation sequence.

Calculation sequence
--------------------
1.  Combustion data (flue-gas flow, adiabatic flame temperature).
2.  Radiant section: bridgewall temperature, Q_radiant, tube temperatures.
3.  Convection section (optional): bank-by-bank heat transfer.
4.  Overall balance: stack temperature, thermal efficiency.
5.  Summary report.

References:
    API STD 560, 4th ed.
    Lobo & Evans (1939) Trans. AIChE 35, 743.
    Wimpress (1963) Hydrocarbon Processing 42(10).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any

from .combustion import combustion_calc, adiabatic_flame_temp, available_heat_fraction, FUEL_DB
from .geometry import FireboxGeometry
from .radiant import calc_radiant, RadiantResult
from .convection import TubeBank, ConvectionBankResult, calc_convection_section
from .properties import process_fluid_cp, PROCESS_FLUID_PRESETS, TUBE_MATERIALS


# ---------------------------------------------------------------------------
# Input data class
# ---------------------------------------------------------------------------

@dataclass
class HeaterInput:
    """
    Complete specification for a fired-heater simulation.
    """
    # ---- identification --------------------------------------------------
    case_name: str = "案例1 / Case 1"
    heater_service: str = "原油加热炉 / Crude Preheater"

    # ---- fuel & combustion -----------------------------------------------
    fuel_type: str = "natural_gas"
    excess_air_pct: float = 15.0       # %  (will be divided by 100)
    W_fuel_kg_s: float = 1.0           # kg/s  fuel mass flow
    T_air_C: float = 15.0              # °C  combustion air temperature

    # ---- firebox geometry ------------------------------------------------
    firebox_type: str = "box"          # "box" | "cabin" | "cylindrical"
    firebox_length_m: float = 10.0
    firebox_width_m: float = 4.0
    firebox_height_m: float = 6.0

    # ---- radiant tubes ---------------------------------------------------
    n_tubes_radiant: int = 24
    rad_tube_od_mm: float = 168.3      # mm → converted to m
    rad_tube_wt_mm: float = 8.0        # mm  wall thickness
    rad_tube_pitch_mm: float = 336.6   # mm  (2 × OD default)
    n_radiant_rows: int = 1
    rad_tube_length_m: float = 10.0
    rad_tube_material: str = "Cr5Mo"

    # ---- process fluid ---------------------------------------------------
    process_fluid: str = "crude_oil"
    W_process_kg_s: float = 50.0
    T_proc_in_C: float = 100.0         # overall cold inlet
    T_proc_out_C: Optional[float] = None  # design outlet; None = calculate

    # ---- convection section banks (list of bank dicts) -------------------
    # Each dict has keys matching TubeBank fields (with _mm suffix for dims)
    conv_banks: List[Dict[str, Any]] = field(default_factory=list)

    # ---- losses & factors ------------------------------------------------
    loss_fraction: float = 0.02        # 2 % heat loss
    peak_flux_factor: float = 1.8


# ---------------------------------------------------------------------------
# Overall result data class
# ---------------------------------------------------------------------------

@dataclass
class HeaterResult:
    """
    Complete fired-heater simulation results.
    """
    # -- identification ----------------------------------------------------
    case_name: str
    heater_service: str

    # -- combustion summary ------------------------------------------------
    fuel_type: str
    fuel_name: str
    LHV_kJ_kg: float
    excess_air_pct: float
    W_fuel_kg_s: float
    W_flue_kg_s: float
    T_AFT_C: float                     # adiabatic flame temperature

    # -- heat duties [kW] --------------------------------------------------
    Q_released_kW: float               # total heat released by fuel
    Q_radiant_kW: float                # absorbed by radiant tubes
    Q_convection_kW: float             # absorbed by convection banks
    Q_total_process_kW: float          # total to process fluid
    Q_losses_kW: float                 # firebox shell losses
    Q_stack_kW: float                  # heat leaving with stack gas

    # -- temperatures [°C] -------------------------------------------------
    T_BWT_C: float                     # bridgewall temperature
    T_stack_C: float                   # flue-gas temperature at stack
    T_proc_out_rad_C: float            # process outlet from radiant
    T_proc_out_conv_C: float           # process outlet from convection (= overall)

    # -- performance -------------------------------------------------------
    thermal_efficiency_pct: float      # (Q_total_process / Q_released) × 100
    available_heat_fraction: float

    # -- radiant section details ------------------------------------------
    radiant: RadiantResult

    # -- convection section details ----------------------------------------
    convection_banks: List[ConvectionBankResult]

    def to_dict(self) -> dict:
        d = {k: v for k, v in self.__dict__.items()
             if k not in ("radiant", "convection_banks")}
        d["radiant"] = self.radiant.to_dict()
        d["convection_banks"] = [b.to_dict() for b in self.convection_banks]
        return d


# ---------------------------------------------------------------------------
# Main simulation function
# ---------------------------------------------------------------------------

def simulate_heater(inp: HeaterInput) -> HeaterResult:
    """
    Run a complete fired-heater simulation.

    Parameters
    ----------
    inp : HeaterInput

    Returns
    -------
    HeaterResult
    """
    excess_air_frac = inp.excess_air_pct / 100.0
    T_air_K = inp.T_air_C + 273.15
    T_ref_K = 288.15  # 15 °C reference

    # ------------------------------------------------------------------
    # Build FireboxGeometry
    # ------------------------------------------------------------------
    geom = FireboxGeometry(
        length=inp.firebox_length_m,
        width=inp.firebox_width_m,
        height=inp.firebox_height_m,
        n_tubes_radiant=inp.n_tubes_radiant,
        tube_od=inp.rad_tube_od_mm / 1000.0,
        tube_wall_t=inp.rad_tube_wt_mm / 1000.0,
        tube_pitch=inp.rad_tube_pitch_mm / 1000.0,
        n_rows=inp.n_radiant_rows,
        tube_length=inp.rad_tube_length_m,
        firebox_type=inp.firebox_type,
    )

    # ------------------------------------------------------------------
    # Process fluid temperature entering radiant section
    # For a heater WITH convection, process fluid goes through
    # convection first (cold side), then radiant (hot side).
    # For a heater WITHOUT convection, process goes directly to radiant.
    # ------------------------------------------------------------------
    has_conv = len(inp.conv_banks) > 0
    if has_conv:
        # Rough estimate of convection duty to find T_process at radiant inlet
        # Will be refined after convection section is calculated.
        T_proc_in_radiant_C = inp.T_proc_in_C + 50.0  # rough guess
    else:
        T_proc_in_radiant_C = inp.T_proc_in_C

    # ------------------------------------------------------------------
    # Radiant section
    # ------------------------------------------------------------------
    rad_result = calc_radiant(
        fuel_type=inp.fuel_type,
        excess_air_frac=excess_air_frac,
        W_fuel_kg_s=inp.W_fuel_kg_s,
        T_air_K=T_air_K,
        geom=geom,
        process_fluid=inp.process_fluid,
        W_process_kg_s=inp.W_process_kg_s,
        T_process_in_C=T_proc_in_radiant_C,
        T_process_out_C=inp.T_proc_out_C,
        tube_material=inp.rad_tube_material,
        loss_fraction=inp.loss_fraction,
        peak_flux_factor=inp.peak_flux_factor,
        T_ref_K=T_ref_K,
    )

    # ------------------------------------------------------------------
    # Convection section
    # ------------------------------------------------------------------
    conv_results: List[ConvectionBankResult] = []
    Q_convection_kW = 0.0
    T_proc_out_conv_C = inp.T_proc_in_C  # default (no convection)

    if has_conv:
        # Build TubeBank objects from input dicts
        banks = []
        for bd in inp.conv_banks:
            tb = TubeBank(
                tube_od=bd.get("tube_od_mm", 114.3) / 1000.0,
                tube_wall_t=bd.get("tube_wt_mm", 8.0) / 1000.0,
                tube_pitch_trans=bd.get("pitch_trans_mm", 220.0) / 1000.0,
                tube_pitch_long=bd.get("pitch_long_mm", 220.0) / 1000.0,
                tube_length=bd.get("tube_length_m", 10.0),
                n_tubes_per_row=bd.get("n_tubes_per_row", 12),
                n_rows=bd.get("n_rows", 4),
                arrangement=bd.get("arrangement", "staggered"),
                fin_height=bd.get("fin_height_mm", None) and bd["fin_height_mm"] / 1000.0,
                fin_thickness=bd.get("fin_thickness_mm", None) and bd["fin_thickness_mm"] / 1000.0,
                fin_pitch=bd.get("fin_pitch_per_m", None),
                fin_material=bd.get("fin_material", "carbon_steel"),
                tube_material=bd.get("tube_material", "carbon_steel"),
                process_fluid=inp.process_fluid,
                W_process_kg_s=inp.W_process_kg_s,
            )
            banks.append(tb)

        conv_results = calc_convection_section(
            banks=banks,
            T_gas_in_K=rad_result.T_BWT_C + 273.15,
            W_flue_kg_s=rad_result.W_flue_kg_s,
            T_proc_in_C=inp.T_proc_in_C,
            excess_air_frac=excess_air_frac,
        )

        Q_convection_kW = sum(r.Q_duty_kW for r in conv_results)

        if conv_results:
            # Process outlet is from the hottest bank (last bank in hot→cold gas order
            # = first bank that process fluid exits toward radiant).
            # In counter-flow: process exits from the bank closest to firebox.
            T_proc_out_conv_C = conv_results[0].T_proc_out_C

    # ------------------------------------------------------------------
    # Overall heat balance
    # ------------------------------------------------------------------
    Q_released = rad_result.Q_released_kW
    Q_radiant = rad_result.Q_radiant_kW
    Q_total_process = Q_radiant + Q_convection_kW
    Q_losses = rad_result.Q_losses_kW
    Q_stack = max(0.0, Q_released - Q_total_process - Q_losses)

    # Stack temperature (from flue-gas enthalpy at stack)
    comb = combustion_calc(inp.fuel_type, excess_air_frac)
    W_flue = rad_result.W_flue_kg_s
    from .combustion import flue_gas_enthalpy_approx
    from scipy.optimize import brentq

    if W_flue > 0 and Q_stack > 0:
        h_stack_per_kg = Q_stack / W_flue  # kJ/kg
        try:
            T_stack_K = brentq(
                lambda T: flue_gas_enthalpy_approx(T, T_ref_K, excess_air_frac) - h_stack_per_kg,
                T_ref_K + 10, T_ref_K + 1500,
                xtol=0.5,
            )
        except ValueError:
            T_stack_K = T_ref_K + 250.0
    else:
        T_stack_K = T_ref_K + 200.0

    T_stack_C = T_stack_K - 273.15

    # If convection section exists, override stack temperature with the
    # exit temperature of the last (coldest) convection bank
    if conv_results:
        T_stack_C = conv_results[-1].T_gas_out_C

    # Thermal efficiency
    eta = Q_total_process / max(Q_released, 1.0)

    # Available heat fraction (stack-loss basis)
    eta_av = available_heat_fraction(
        inp.fuel_type, excess_air_frac, T_stack_K, T_air_K, T_ref_K
    )

    fuel_info = FUEL_DB.get(inp.fuel_type, {})

    return HeaterResult(
        case_name=inp.case_name,
        heater_service=inp.heater_service,
        fuel_type=inp.fuel_type,
        fuel_name=fuel_info.get("name", inp.fuel_type),
        LHV_kJ_kg=fuel_info.get("LHV", 0.0),
        excess_air_pct=inp.excess_air_pct,
        W_fuel_kg_s=round(rad_result.W_fuel_kg_s, 4),
        W_flue_kg_s=round(W_flue, 4),
        T_AFT_C=round(rad_result.T_AFT_C, 1),
        Q_released_kW=round(Q_released, 1),
        Q_radiant_kW=round(Q_radiant, 1),
        Q_convection_kW=round(Q_convection_kW, 1),
        Q_total_process_kW=round(Q_total_process, 1),
        Q_losses_kW=round(Q_losses, 1),
        Q_stack_kW=round(Q_stack, 1),
        T_BWT_C=rad_result.T_BWT_C,
        T_stack_C=round(T_stack_C, 1),
        T_proc_out_rad_C=rad_result.T_process_out_C,
        T_proc_out_conv_C=round(T_proc_out_conv_C, 1),
        thermal_efficiency_pct=round(eta * 100.0, 2),
        available_heat_fraction=round(eta_av, 4),
        radiant=rad_result,
        convection_banks=conv_results,
    )
