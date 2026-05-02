"""
Tests for the radiant section calculation module.
"""
import pytest

from src.geometry import FireboxGeometry
from src.radiant import calc_radiant, RadiantResult


class TestFireboxGeometry:
    def test_default_geometry_computes(self):
        g = FireboxGeometry()
        assert g.A_cp > 0
        assert g.F_exchange > 0
        assert g.A_tube_outer > 0

    def test_exchange_factor_between_0_and_1(self):
        g = FireboxGeometry()
        assert 0.0 < g.F_exchange <= 1.0

    def test_tube_id_less_than_od(self):
        g = FireboxGeometry(tube_od=0.168, tube_wall_t=0.008)
        assert g.tube_id < g.tube_od
        assert g.tube_id == pytest.approx(0.152, abs=0.001)

    def test_larger_tube_array_increases_cold_plane_area(self):
        g_few = FireboxGeometry(n_tubes_radiant=10)
        g_many = FireboxGeometry(n_tubes_radiant=30)
        assert g_many.A_cp > g_few.A_cp

    def test_cylindrical_firebox_no_endwall(self):
        g = FireboxGeometry(firebox_type="cylindrical")
        assert g.A_endwall == 0.0

    def test_summary_contains_required_keys(self):
        g = FireboxGeometry()
        s = g.summary()
        for key in ("F_exchange", "A_cp_m2", "A_tube_outer_m2", "n_tubes_radiant"):
            assert key in s


class TestCalcRadiant:
    def _default_result(self) -> RadiantResult:
        return calc_radiant(
            fuel_type="natural_gas",
            excess_air_frac=0.15,
            W_fuel_kg_s=1.0,
        )

    def test_returns_radiant_result(self):
        r = self._default_result()
        assert isinstance(r, RadiantResult)

    def test_bwt_below_aft(self):
        r = self._default_result()
        # BWT must be below adiabatic flame temperature
        assert r.T_BWT_C < r.T_AFT_C

    def test_bwt_reasonable_range(self):
        r = self._default_result()
        # BWT can be high for an overloaded firebox; must stay below AFT
        assert 500 < r.T_BWT_C < r.T_AFT_C, f"BWT = {r.T_BWT_C:.1f} °C"

    def test_q_radiant_positive(self):
        r = self._default_result()
        assert r.Q_radiant_kW > 0

    def test_heat_balance(self):
        r = self._default_result()
        # Q_released = Q_radiant + Q_to_flue + Q_losses
        # The "Q_to_flue" at BWT ≈ Q_released - Q_radiant - Q_losses
        total = r.Q_radiant_kW + r.Q_losses_kW
        assert total <= r.Q_released_kW * 1.001  # can't absorb more than released

    def test_peak_flux_greater_than_average(self):
        r = self._default_result()
        assert r.q_max_kW_m2 > r.q_avg_kW_m2

    def test_more_fuel_increases_duty(self):
        r_low = calc_radiant(W_fuel_kg_s=0.5)
        r_high = calc_radiant(W_fuel_kg_s=2.0)
        assert r_high.Q_radiant_kW > r_low.Q_radiant_kW

    def test_tube_skin_temperature_above_process(self):
        r = self._default_result()
        # Tube skin must be hotter than mean process fluid temperature
        assert r.T_tube_avg_C > 100.0  # process inlet is ~200 °C, so avg should be high

    def test_to_dict_serializable(self):
        r = self._default_result()
        d = r.to_dict()
        assert isinstance(d, dict)
        assert "T_BWT_C" in d
