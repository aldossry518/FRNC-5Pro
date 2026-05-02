"""
对流室计算模块  –  Convection section heat transfer
=====================================================
Calculates heat transfer for one or more tube banks in the convection section
of a process fired heater.

Tube banks may be:
    - bare tubes (cross-flow correlation)
    - internally finned (not addressed here)
    - externally finned – segmented or solid fins

Outer (gas-side) heat transfer coefficient:
    Grimison / Zukauskas correlation for tube banks in cross-flow.
    Churchill-Bernstein alternative for a single tube row.

Overall heat transfer coefficient (referred to outer area):
    1/U_o = 1/h_o + (A_o/A_i) * 1/h_i + (A_o * r_w) / A_lm
where r_w = tube-wall resistance = ln(r_o/r_i) / (2π k_m L)  per unit length.

Fin efficiency (for finned tubes):
    η_f = tanh(m * l_f) / (m * l_f)
    m   = sqrt(2 h_o / (k_fin * t_fin))
    Effective outer area: A_o_eff = A_bare + η_f * A_fins

References:
    Grimison E.D. (1937) Trans. ASME 59, 583.
    Zukauskas A. (1972) Adv. Heat Transfer 8, 93.
    API STD 560, 4th ed.
    Kern D.Q. (1950) "Process Heat Transfer", McGraw-Hill.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, List

from .properties import (
    flue_gas_viscosity,
    flue_gas_thermal_conductivity,
    flue_gas_density,
    flue_gas_prandtl,
    inner_film_htc,
    tube_metal_k,
    TUBE_MATERIALS,
)
from .combustion import flue_gas_cp, flue_gas_enthalpy_approx


# ---------------------------------------------------------------------------
# Tube-bank data class
# ---------------------------------------------------------------------------

@dataclass
class TubeBank:
    """
    Single tube bank (row group) in the convection section.

    Dimensions in SI units (m, m², W, kW …).
    """
    # ---- tube geometry -------------------------------------------------
    tube_od: float = 0.114          # m  (4-inch NPS ≈ 0.1143 m)
    tube_wall_t: float = 0.008      # m
    tube_pitch_trans: float = 0.200  # m  transverse pitch  S_T
    tube_pitch_long: float = 0.200   # m  longitudinal pitch S_L
    tube_length: float = 10.0       # m  effective heated length
    n_tubes_per_row: int = 12       # tubes per row (across duct width)
    n_rows: int = 4                 # number of rows deep

    # ---- tube arrangement: "staggered" or "inline" ---------------------
    arrangement: str = "staggered"

    # ---- fin geometry (set to None for bare tubes) ----------------------
    fin_height: Optional[float] = None  # m  (radial fin height)
    fin_thickness: Optional[float] = None  # m
    fin_pitch: Optional[float] = None    # fins / m  (linear density)
    fin_material: str = "carbon_steel"

    # ---- metallurgy ----------------------------------------------------
    tube_material: str = "carbon_steel"

    # ---- process-side (shell-and-tube inside tubes) --------------------
    process_fluid: str = "crude_oil"
    W_process_kg_s: float = 50.0

    # ---- computed fields -----------------------------------------------
    tube_id: float = field(init=False)
    n_tubes_total: int = field(init=False)
    A_tube_outer: float = field(init=False)  # m²  bare outer surface
    A_tube_inner: float = field(init=False)  # m²  inner surface
    A_outer_eff: float = field(init=False)   # m²  effective outer (with fins)
    fin_efficiency: float = field(init=False)

    def __post_init__(self) -> None:
        self.tube_id = self.tube_od - 2.0 * self.tube_wall_t
        self.n_tubes_total = self.n_tubes_per_row * self.n_rows
        # Bare outer area
        self.A_tube_outer = (
            self.n_tubes_total * math.pi * self.tube_od * self.tube_length
        )
        # Inner area
        self.A_tube_inner = (
            self.n_tubes_total * math.pi * self.tube_id * self.tube_length
        )
        # Fin area computed later (needs h_o)
        self.fin_efficiency = 1.0
        self.A_outer_eff = self.A_tube_outer


@dataclass
class ConvectionBankResult:
    """Heat transfer result for a single tube bank."""
    bank_index: int
    T_gas_in_C: float
    T_gas_out_C: float
    T_proc_in_C: float
    T_proc_out_C: float
    Q_duty_kW: float
    h_o_W_m2K: float          # gas-side (outer) HTC  [W/(m²·K)]
    h_i_W_m2K: float          # process-side (inner) HTC  [W/(m²·K)]
    U_o_W_m2K: float          # overall HTC referred to outer area
    LMTD_K: float
    A_outer_eff_m2: float
    T_tube_max_C: float       # peak tube skin temperature

    def to_dict(self) -> dict:
        return {k: round(v, 2) if isinstance(v, float) else v
                for k, v in self.__dict__.items()}


# ---------------------------------------------------------------------------
# Main convection calculation function
# ---------------------------------------------------------------------------

def calc_convection_bank(
    bank: TubeBank,
    T_gas_in_K: float,
    W_flue_kg_s: float,
    T_proc_in_C: float,
    excess_air_frac: float = 0.15,
    bank_index: int = 0,
) -> ConvectionBankResult:
    """
    Calculate heat transfer for one tube bank.

    Uses the effectiveness-NTU method to solve for both outlet temperatures.

    Parameters
    ----------
    bank : TubeBank
        Tube bank geometry and conditions.
    T_gas_in_K : float
        Flue-gas temperature entering this bank [K].
    W_flue_kg_s : float
        Flue-gas mass flow [kg/s].
    T_proc_in_C : float
        Process fluid temperature entering this bank [°C].
    excess_air_frac : float
        Excess air fraction for flue-gas properties.
    bank_index : int
        For reporting.

    Returns
    -------
    ConvectionBankResult
    """
    T_gas_in_C = T_gas_in_K - 273.15
    T_proc_in_K = T_proc_in_C + 273.15

    # ---- Flue-gas properties at mean temperature (estimate) ----------
    T_gas_mean_K = T_gas_in_K - 50.0  # rough estimate; refined below

    rho_fg = flue_gas_density(T_gas_mean_K)            # kg/m³
    mu_fg = flue_gas_viscosity(T_gas_mean_K)           # Pa·s
    lam_fg = flue_gas_thermal_conductivity(T_gas_mean_K)  # W/(m·K)
    cp_fg = flue_gas_cp(T_gas_mean_K, excess_air_frac) * 1000.0  # J/(kg·K)
    Pr_fg = mu_fg * cp_fg / lam_fg

    # ---- Gas-side velocity through minimum free section ---------------
    # Duct width ≈ n_tubes_per_row * S_T (transverse pitch)
    duct_width = bank.n_tubes_per_row * bank.tube_pitch_trans
    duct_height = bank.tube_length   # along tube axis = gas flow height
    # Minimum free area (between tubes in transverse direction)
    A_min_free = (
        bank.n_tubes_per_row
        * (bank.tube_pitch_trans - bank.tube_od)
        * bank.tube_length
    )
    if A_min_free <= 0:
        A_min_free = duct_width * duct_height * 0.5

    # Mass velocity through minimum free section
    G_fg = W_flue_kg_s / A_min_free   # kg/(m²·s)
    Re_fg = G_fg * bank.tube_od / mu_fg

    # ---- Outer (gas-side) heat transfer coefficient -------------------
    h_o = _grimison_htc(
        Re=Re_fg,
        Pr=Pr_fg,
        lam=lam_fg,
        d_o=bank.tube_od,
        arrangement=bank.arrangement,
        S_T=bank.tube_pitch_trans,
        S_L=bank.tube_pitch_long,
        n_rows=bank.n_rows,
    )

    # ---- Fin area and efficiency (if finned) --------------------------
    A_outer_eff = bank.A_tube_outer
    eta_fin = 1.0
    if bank.fin_height is not None and bank.fin_thickness is not None and bank.fin_pitch is not None:
        A_outer_eff, eta_fin = _fin_area(bank, h_o)
    bank.A_outer_eff = A_outer_eff
    bank.fin_efficiency = eta_fin

    # ---- Inner (process-side) HTC ------------------------------------
    A_c_inner = (
        bank.n_tubes_total * math.pi * bank.tube_id ** 2 / 4.0
    )
    G_proc = bank.W_process_kg_s / max(A_c_inner, 1e-9)
    h_i = inner_film_htc(
        bank.process_fluid, G_proc, bank.tube_id, T_proc_in_K + 50.0
    )

    # ---- Tube-wall resistance ----------------------------------------
    r_o = bank.tube_od / 2.0
    r_i = bank.tube_id / 2.0
    k_m = tube_metal_k(bank.tube_material, T_proc_in_K + 50.0)
    # Wall resistance referred to outer area [m²·K/W]:
    R_wall = (r_o * math.log(r_o / r_i)) / k_m
    # Area ratio
    A_ratio = bank.A_tube_outer / bank.A_tube_inner  # ≈ d_o / d_i

    # ---- Overall HTC referred to effective outer area ----------------
    U_o = 1.0 / (1.0 / h_o + R_wall + A_ratio / h_i)

    # ---- Heat capacities ---------------------------------------------
    from .properties import process_fluid_cp as pf_cp
    cp_proc = pf_cp(bank.process_fluid) * 1000.0  # J/(kg·K)

    C_hot = W_flue_kg_s * cp_fg                   # W/K  flue gas
    C_cold = bank.W_process_kg_s * cp_proc          # W/K  process fluid
    C_min = min(C_hot, C_cold)
    C_max = max(C_hot, C_cold)
    C_ratio = C_min / C_max if C_max > 0 else 0.0

    # ---- Effectiveness-NTU (counter-flow) ---------------------------
    UA = U_o * A_outer_eff  # W/K
    NTU = UA / max(C_min, 1.0)
    eps = _effectiveness_counterflow(NTU, C_ratio)

    Q_max = C_min * (T_gas_in_K - T_proc_in_K)  # W
    Q_duty = eps * Q_max                          # W

    # Outlet temperatures
    T_gas_out_K = T_gas_in_K - Q_duty / max(C_hot, 1.0)
    T_proc_out_K = T_proc_in_K + Q_duty / max(C_cold, 1.0)

    # ---- LMTD (counter-flow) -----------------------------------------
    dT1 = T_gas_in_K - T_proc_out_K
    dT2 = T_gas_out_K - T_proc_in_K
    if dT1 > 0 and dT2 > 0 and abs(dT1 - dT2) > 0.1:
        LMTD = (dT1 - dT2) / math.log(dT1 / dT2)
    else:
        LMTD = (dT1 + dT2) / 2.0

    # ---- Peak tube skin temperature ----------------------------------
    q_outer = Q_duty / max(A_outer_eff, 1e-3)  # W/m²
    T_tube_max_K = T_proc_out_K + q_outer / max(h_i, 1.0) + q_outer * R_wall

    return ConvectionBankResult(
        bank_index=bank_index,
        T_gas_in_C=round(T_gas_in_C, 1),
        T_gas_out_C=round(T_gas_out_K - 273.15, 1),
        T_proc_in_C=round(T_proc_in_C, 1),
        T_proc_out_C=round(T_proc_out_K - 273.15, 1),
        Q_duty_kW=round(Q_duty * 1e-3, 1),
        h_o_W_m2K=round(h_o, 1),
        h_i_W_m2K=round(h_i, 1),
        U_o_W_m2K=round(U_o, 1),
        LMTD_K=round(LMTD, 1),
        A_outer_eff_m2=round(A_outer_eff, 2),
        T_tube_max_C=round(T_tube_max_K - 273.15, 1),
    )


def calc_convection_section(
    banks: List[TubeBank],
    T_gas_in_K: float,           # bridgewall temperature from radiant
    W_flue_kg_s: float,
    T_proc_in_C: float,           # process fluid inlet to convection (cold end)
    excess_air_frac: float = 0.15,
) -> List[ConvectionBankResult]:
    """
    Calculate the complete convection section (multiple banks).

    Flue gas flows top-to-bottom (hot to cold).
    Process fluid flows bottom-to-top (cold to hot) – counter-current.

    Parameters
    ----------
    banks : list of TubeBank
        Ordered hot-end → cold-end (flue gas direction).
    T_gas_in_K : float
        Flue-gas inlet temperature (= bridgewall temperature) [K].
    W_flue_kg_s : float
        Flue-gas mass flow [kg/s].
    T_proc_in_C : float
        Process fluid cold-end inlet temperature [°C].
    excess_air_frac : float

    Returns
    -------
    list of ConvectionBankResult
        One entry per bank, in flue-gas direction (hot → cold).
    """
    results: List[ConvectionBankResult] = []
    T_gas_K = T_gas_in_K

    # Process fluid flows counter-current: start at cold end.
    # We need to find process inlet to each bank.  Use simple sequential
    # calculation, marching in the flue-gas direction (hot→cold) and
    # accumulating process heat from the cold end.

    # --- collect total process W and cp ----------------------------------
    # Assume same process fluid and flow through all banks for simplicity
    from .properties import process_fluid_cp as pf_cp
    first_bank = banks[0] if banks else None
    if first_bank is None:
        return results

    cp_proc = pf_cp(first_bank.process_fluid) * 1000.0  # J/(kg·K)
    W_proc = first_bank.W_process_kg_s
    C_cold_total = W_proc * cp_proc  # W/K

    # First pass: calculate duty for each bank assuming a rough T_proc profile
    # (counter-flow multi-pass approximation: iterate once)
    T_proc_K = T_proc_in_C + 273.15  # cold end

    # For counter-flow, process flow is opposite to gas flow.
    # Process enters at cold end of last bank, exits hot end of first bank.
    # We march in flue-gas direction, but process temperature increases
    # in the opposite direction.

    # Simple approach: single-pass, assuming process inlet to each bank
    # is computed from cumulative heat added (cold→hot direction).
    # This is exact for equal flow rates; approximate otherwise.

    # Compute Q per bank from hot end to cold end
    bank_results_pass1 = []
    T_gas_K_cur = T_gas_K
    T_proc_K_cur = T_proc_in_C + 273.15  # this will track the cold-end side

    for idx, bank in enumerate(reversed(banks)):
        # Cold-end bank first
        res = calc_convection_bank(
            bank,
            T_gas_K_cur + (idx + 1) * 0.0,  # placeholder for sequential T_gas
            W_flue_kg_s,
            T_proc_K_cur - 273.15,
            excess_air_frac,
            bank_index=len(banks) - 1 - idx,
        )
        bank_results_pass1.append(res)
        T_proc_K_cur += res.Q_duty_kW * 1000.0 / max(C_cold_total, 1.0)

    # Now do the proper hot-end-first pass with correct process temperatures
    # Process outlet from convection = T_proc_K_cur (from cold end accumulation)
    T_proc_hot_exit_K = T_proc_K_cur

    T_proc_bank_inlet: List[float] = []
    T_proc_cur = T_proc_hot_exit_K
    for res in bank_results_pass1:
        T_proc_cur -= res.Q_duty_kW * 1000.0 / max(C_cold_total, 1.0)
        T_proc_bank_inlet.append(T_proc_cur)
    T_proc_bank_inlet.reverse()  # now in cold→hot order (matching reversed bank list)
    # Re-order to hot→cold
    T_proc_bank_inlet_hot_first = list(reversed(T_proc_bank_inlet))

    T_gas_K_cur = T_gas_K
    for idx, bank in enumerate(banks):
        T_proc_in_this = T_proc_bank_inlet_hot_first[idx] - 273.15
        res = calc_convection_bank(
            bank,
            T_gas_K_cur,
            W_flue_kg_s,
            T_proc_in_this,
            excess_air_frac,
            bank_index=idx,
        )
        results.append(res)
        T_gas_K_cur = res.T_gas_out_C + 273.15

    return results


# ---------------------------------------------------------------------------
# Heat-transfer correlations
# ---------------------------------------------------------------------------

def _grimison_htc(
    Re: float,
    Pr: float,
    lam: float,
    d_o: float,
    arrangement: str = "staggered",
    S_T: float = 0.20,
    S_L: float = 0.20,
    n_rows: int = 4,
) -> float:
    """
    Gas-side heat transfer coefficient using the Grimison / Zukauskas correlation.

    Valid for 1000 ≤ Re ≤ 2×10^5 for tube banks in cross-flow.

    Returns h_o in [W / (m²·K)].
    """
    if Re < 100:
        # Very low Re – use Churchill-Bernstein for a single cylinder
        Nu = 0.3 + 0.62 * Re ** 0.5 * Pr ** (1.0 / 3.0) / (
            1.0 + (0.4 / Pr) ** (2.0 / 3.0)
        ) ** 0.25
    else:
        x_T = S_T / d_o
        x_L = S_L / d_o

        if arrangement == "staggered":
            # Zukauskas staggered
            C1 = 0.35 * (x_T / x_L) ** 0.2 if x_T / x_L < 2.0 else 0.40
            if Re < 1000:
                m = 0.5
            elif Re < 2e5:
                m = 0.6
            else:
                m = 0.63
        else:
            # Inline
            C1 = 0.27
            if Re < 1000:
                m = 0.4
            elif Re < 2e5:
                m = 0.63
            else:
                m = 0.7

        Nu = C1 * Re ** m * Pr ** 0.36

        # Row correction for n_rows < 10
        if n_rows < 10:
            F_row = _row_correction(arrangement, n_rows)
            Nu *= F_row

    h_o = Nu * lam / d_o
    return max(20.0, h_o)


def _row_correction(arrangement: str, n_rows: int) -> float:
    """
    Correction factor for tube banks with fewer than 10 rows.

    Values from Incropera & DeWitt, Table 7.2.
    """
    # Tabulated correction factors (Zukauskas, 1972)
    staggered = {1: 0.68, 2: 0.75, 3: 0.83, 4: 0.89, 5: 0.92,
                 6: 0.95, 7: 0.97, 8: 0.98, 9: 0.99, 10: 1.00}
    inline = {1: 0.64, 2: 0.80, 3: 0.87, 4: 0.90, 5: 0.92,
              6: 0.94, 7: 0.96, 8: 0.98, 9: 0.99, 10: 1.00}
    table = staggered if arrangement == "staggered" else inline
    return table.get(min(n_rows, 10), 1.00)


def _fin_area(bank: TubeBank, h_o: float) -> tuple[float, float]:
    """
    Compute effective outer area and fin efficiency for a finned tube bank.

    Returns
    -------
    (A_outer_eff, eta_fin)
    """
    from .properties import tube_metal_k as k_metal
    l_f = bank.fin_height                        # m
    t_f = bank.fin_thickness                     # m
    fins_per_m = bank.fin_pitch                  # fins / m
    k_fin = k_metal(bank.fin_material, 700.0)    # W/(m·K)

    r_i = bank.tube_od / 2.0
    r_o = r_i + l_f

    # Fin efficiency (annular fin approximation, Kern)
    m = math.sqrt(2.0 * h_o / (k_fin * t_f))
    ml = m * l_f
    if ml < 0.01:
        eta_fin = 1.0
    else:
        eta_fin = math.tanh(ml) / ml

    # Areas per unit tube length
    L_tube = bank.tube_length
    n_t = bank.n_tubes_total

    # Bare area between fins (per tube per metre)
    fin_pitch_m = 1.0 / fins_per_m   # m per fin
    fin_spacing = fin_pitch_m - t_f   # m between fins
    A_bare_per_m = math.pi * bank.tube_od * fin_spacing * fins_per_m
    # Fin surface area per metre (two faces + tip)
    A_fin_per_m = fins_per_m * (
        2.0 * math.pi * (r_o ** 2 - r_i ** 2)
        + 2.0 * math.pi * r_o * t_f
    )

    A_bare_total = n_t * A_bare_per_m * L_tube
    A_fin_total = n_t * A_fin_per_m * L_tube

    A_outer_eff = A_bare_total + eta_fin * A_fin_total
    return A_outer_eff, eta_fin


def _effectiveness_counterflow(NTU: float, C_ratio: float) -> float:
    """
    Heat exchanger effectiveness for counter-flow arrangement.

    ε = [1 - exp(-NTU(1-C*))] / [1 - C* exp(-NTU(1-C*))]
    Special case C* = 1: ε = NTU / (1 + NTU)
    """
    if C_ratio >= 1.0 or abs(C_ratio - 1.0) < 1e-6:
        return NTU / (1.0 + NTU)
    exp_term = math.exp(-NTU * (1.0 - C_ratio))
    return (1.0 - exp_term) / (1.0 - C_ratio * exp_term)
