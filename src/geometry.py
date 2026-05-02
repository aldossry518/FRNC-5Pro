"""
几何计算模块  –  Furnace geometry helpers
==========================================
Firebox dimensions → areas, cold-plane areas, and exchange factors for the
Lobo-Evans single-zone radiant-section model.

Firebox types supported:
    "box"       – rectangular box (vertical or horizontal)
    "cabin"     – cabin (A-frame / double-fired)
    "cylindrical" – vertical cylindrical

References:
    Lobo W.E. & Evans J.E. (1939) Trans. AIChE 35, 743.
    API STD 560, 4th ed., Annex A.
    Wimpress R.N. (1963) Hydrocarbon Processing & Petroleum Refiner 42(10).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class FireboxGeometry:
    """
    Firebox geometry descriptor.

    All linear dimensions in metres.
    """
    # ---- dimensions --------------------------------------------------------
    length: float = 10.0    # m  (horizontal firebox: tube length + end-wall)
    width: float = 4.0      # m
    height: float = 6.0     # m

    # ---- tube arrangement (radiant section) --------------------------------
    n_tubes_radiant: int = 24       # total number of radiant tubes
    tube_od: float = 0.168          # m  (e.g. 6-inch NPS ≈ 0.1683 m)
    tube_wall_t: float = 0.008      # m  tube wall thickness
    tube_pitch: float = 0.336       # m  centre-to-centre spacing (≥ tube_od)
    n_rows: int = 1                 # number of tube rows (radiant)

    # ---- tube arrangement (convection section) -----------------------------
    # (handled separately in ConvectionBank, but stored here for completeness)
    tube_length: float = 10.0       # m  effective heated tube length

    # ---- surface emissivities ----------------------------------------------
    alpha_tube: float = 0.90        # tube-surface absorptivity (grey body)
    alpha_refractory: float = 0.50  # mean refractory-wall emissivity

    # ---- firebox type ------------------------------------------------------
    firebox_type: str = "box"       # "box" | "cabin" | "cylindrical"

    # ---- convection section offset from firebox ----------------------------
    has_convection: bool = True

    # ---- computed (filled by post_init) ------------------------------------
    tube_id: float = field(init=False)
    A_floor: float = field(init=False)
    A_roof: float = field(init=False)
    A_sidewall: float = field(init=False)
    A_endwall: float = field(init=False)
    A_refractory: float = field(init=False)
    A_cp: float = field(init=False)   # cold-plane area of radiant tubes [m²]
    F_exchange: float = field(init=False)  # Lobo-Evans exchange factor [-]

    def __post_init__(self) -> None:
        self.tube_id = self.tube_od - 2.0 * self.tube_wall_t
        self._compute_areas()

    def _compute_areas(self) -> None:
        L = self.length
        W = self.width
        H = self.height

        if self.firebox_type == "box":
            self.A_floor = L * W
            self.A_roof = L * W
            self.A_sidewall = 2.0 * L * H
            self.A_endwall = 2.0 * W * H
        elif self.firebox_type == "cabin":
            # Cabin (A-frame): two angled roof panels + side walls + ends
            roof_angle = math.radians(30.0)
            half_span = W / 2.0
            rafter = half_span / math.cos(roof_angle)
            self.A_floor = L * W
            self.A_roof = 2.0 * L * rafter   # two roof panels
            self.A_sidewall = 2.0 * L * H
            self.A_endwall = 2.0 * (W * H + half_span * math.tan(roof_angle) * W / 2.0)
        elif self.firebox_type == "cylindrical":
            # Vertical cylindrical: diameter = width, height = height
            D = W
            self.A_floor = math.pi * D ** 2 / 4.0
            self.A_roof = self.A_floor
            self.A_sidewall = math.pi * D * H   # lateral surface
            self.A_endwall = 0.0
        else:
            raise ValueError(f"Unknown firebox_type: {self.firebox_type!r}")

        A_total_shell = (
            self.A_floor + self.A_roof + self.A_sidewall + self.A_endwall
        )

        # ---- cold-plane area of the radiant tube bundle ---------------------
        # Cold-plane area = projected area of the tube row(s) facing the gas
        # For n rows against a refractory wall:
        #   A_cp = N_tubes × tube_length × tube_od   (single row)
        A_cp_1row = self.n_tubes_radiant * self.tube_length * self.tube_od
        self.A_cp = A_cp_1row  # consistent with single/multiple row correction below

        # ---- absorption factor of tube bundle (rows) -----------------------
        # Single-row cold-plane factor (fraction of pitch covered by tube):
        alpha_s = _single_row_absorptivity(
            self.tube_od, self.tube_pitch, self.alpha_tube
        )
        # For n_rows stacked:
        alpha_n = 1.0 - (1.0 - alpha_s) ** self.n_rows

        # ---- refractory area not shielded by tubes -------------------------
        # Tubes are typically arranged along one or two walls.
        # Simplified: tubes shield one side wall; rest is refractory.
        A_tube_shield = self.n_tubes_radiant * self.tube_length * self.tube_pitch
        # The remaining refractory surface:
        self.A_refractory = max(0.0, A_total_shell - A_tube_shield)

        # ---- Lobo-Evans exchange factor ------------------------------------
        self.F_exchange = _lobo_evans_F(
            alpha_n, self.A_cp, self.A_refractory, self.alpha_refractory
        )

    # ------------------------------------------------------------------ helpers

    @property
    def A_tube_outer(self) -> float:
        """Total outer surface area of radiant tubes [m²]."""
        return self.n_tubes_radiant * math.pi * self.tube_od * self.tube_length

    @property
    def A_tube_inner(self) -> float:
        """Total inner surface area of radiant tubes [m²]."""
        return self.n_tubes_radiant * math.pi * self.tube_id * self.tube_length

    def summary(self) -> dict:
        return {
            "firebox_type": self.firebox_type,
            "length_m": self.length,
            "width_m": self.width,
            "height_m": self.height,
            "A_floor_m2": round(self.A_floor, 2),
            "A_roof_m2": round(self.A_roof, 2),
            "A_sidewall_m2": round(self.A_sidewall, 2),
            "A_endwall_m2": round(self.A_endwall, 2),
            "A_refractory_m2": round(self.A_refractory, 2),
            "A_cp_m2": round(self.A_cp, 2),
            "A_tube_outer_m2": round(self.A_tube_outer, 2),
            "F_exchange": round(self.F_exchange, 4),
            "n_tubes_radiant": self.n_tubes_radiant,
            "tube_od_mm": round(self.tube_od * 1000, 1),
            "tube_id_mm": round(self.tube_id * 1000, 1),
            "tube_pitch_mm": round(self.tube_pitch * 1000, 1),
            "n_rows": self.n_rows,
        }


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _single_row_absorptivity(d: float, s: float, alpha_tube: float) -> float:
    """
    Effective absorptivity of a single row of tubes in front of a refractory wall.

    Uses the Truelove / Lobo-Evans formula for view-factor-weighted absorption:
      alpha_s = (d/s) * alpha_tube   (simplified plane-surface approximation)
    with correction for refractory re-radiation.

    Parameters
    ----------
    d : float  tube outer diameter [m]
    s : float  tube centre-to-centre pitch [m]
    alpha_tube : float  tube-surface absorptivity
    """
    if s <= 0 or d <= 0:
        return 0.0
    # Direct absorption fraction (fraction of flat plane subtended by tubes)
    f_direct = min(1.0, d / s)
    # With refractory backing (alpha_ref assumed ~ 0.5), some radiation is
    # re-emitted and absorbed by the tube row.
    # Effective single-row absorptivity (Wimpress, 1963):
    alpha_s = f_direct * alpha_tube + (1.0 - f_direct) * alpha_tube * 0.50
    return min(1.0, alpha_s)


def _lobo_evans_F(
    alpha_n: float,
    A_cp: float,
    A_ref: float,
    alpha_ref: float = 0.50,
) -> float:
    """
    Lobo-Evans / Hottel exchange factor F for the radiant section.

    F = 1 / [1/alpha_n + A_ref/A_cp * (1/alpha_ref - 1) * alpha_n/1]

    Simplified common form (API 560 Annex A, Eq. A.1):
        F = alpha_n / (alpha_n + A_ref/A_cp * (1 - alpha_n))

    Here alpha_ref is effectively lumped into the refractory area correction.

    Parameters
    ----------
    alpha_n : float  effective tube-row absorptivity
    A_cp    : float  cold-plane area [m²]
    A_ref   : float  refractory area [m²]
    alpha_ref : float  refractory emissivity
    """
    if A_cp <= 0:
        return 0.0
    # Lobo-Evans formula (with refractory contribution):
    ratio = A_ref / A_cp
    # Exchange factor accounting for refractory walls
    F = alpha_n / (alpha_n + ratio * (1.0 - alpha_n) * (1.0 - alpha_ref))
    return min(1.0, max(0.01, F))
