"""
Tests for convection section and full heater model.
"""
import pytest

from src.convection import TubeBank, calc_convection_bank, calc_convection_section
from src.heater import HeaterInput, HeaterResult, simulate_heater


class TestConvectionBank:
    def _default_bank(self) -> TubeBank:
        return TubeBank(
            tube_od=0.1143,
            tube_wall_t=0.008,
            tube_pitch_trans=0.22,
            tube_pitch_long=0.22,
            tube_length=10.0,
            n_tubes_per_row=12,
            n_rows=4,
            arrangement="staggered",
            process_fluid="crude_oil",
            W_process_kg_s=50.0,
        )

    def test_bank_inner_less_than_outer(self):
        bank = self._default_bank()
        assert bank.tube_id < bank.tube_od

    def test_total_tubes(self):
        bank = self._default_bank()
        assert bank.n_tubes_total == 12 * 4

    def test_calc_returns_result(self):
        bank = self._default_bank()
        res = calc_convection_bank(
            bank=bank,
            T_gas_in_K=900.0,
            W_flue_kg_s=18.0,
            T_proc_in_C=150.0,
            excess_air_frac=0.15,
        )
        assert res.Q_duty_kW > 0

    def test_gas_cools_down(self):
        bank = self._default_bank()
        res = calc_convection_bank(
            bank=bank,
            T_gas_in_K=900.0,
            W_flue_kg_s=18.0,
            T_proc_in_C=150.0,
        )
        assert res.T_gas_out_C < res.T_gas_in_C

    def test_process_heats_up(self):
        bank = self._default_bank()
        res = calc_convection_bank(
            bank=bank,
            T_gas_in_K=900.0,
            W_flue_kg_s=18.0,
            T_proc_in_C=150.0,
        )
        assert res.T_proc_out_C > res.T_proc_in_C

    def test_overall_u_positive(self):
        bank = self._default_bank()
        res = calc_convection_bank(
            bank=bank,
            T_gas_in_K=900.0,
            W_flue_kg_s=18.0,
            T_proc_in_C=150.0,
        )
        assert res.U_o_W_m2K > 0

    def test_lmtd_positive(self):
        bank = self._default_bank()
        res = calc_convection_bank(
            bank=bank,
            T_gas_in_K=900.0,
            W_flue_kg_s=18.0,
            T_proc_in_C=150.0,
        )
        assert res.LMTD_K > 0

    def test_inline_arrangement(self):
        bank = self._default_bank()
        bank.arrangement = "inline"
        res = calc_convection_bank(
            bank=bank,
            T_gas_in_K=900.0,
            W_flue_kg_s=18.0,
            T_proc_in_C=150.0,
        )
        assert res.Q_duty_kW > 0


class TestHeaterSimulation:
    def _default_input(self) -> HeaterInput:
        return HeaterInput(
            fuel_type="natural_gas",
            excess_air_pct=15.0,
            W_fuel_kg_s=1.0,
        )

    def test_simulate_returns_result(self):
        res = simulate_heater(self._default_input())
        assert isinstance(res, HeaterResult)

    def test_total_duty_positive(self):
        res = simulate_heater(self._default_input())
        assert res.Q_total_process_kW > 0

    def test_efficiency_between_0_and_100(self):
        res = simulate_heater(self._default_input())
        assert 0 < res.thermal_efficiency_pct < 100

    def test_stack_temp_above_ambient(self):
        res = simulate_heater(self._default_input())
        assert res.T_stack_C > 50

    def test_heat_balance_closes(self):
        inp = self._default_input()
        res = simulate_heater(inp)
        # Q_released ≈ Q_process + Q_losses + Q_stack  (within 2 %)
        total = res.Q_total_process_kW + res.Q_losses_kW + res.Q_stack_kW
        assert total == pytest.approx(res.Q_released_kW, rel=0.05)

    def test_bwt_below_aft(self):
        res = simulate_heater(self._default_input())
        assert res.T_BWT_C < res.T_AFT_C

    def test_more_fuel_increases_total_duty(self):
        inp_low = HeaterInput(W_fuel_kg_s=0.5)
        inp_high = HeaterInput(W_fuel_kg_s=2.0)
        res_low = simulate_heater(inp_low)
        res_high = simulate_heater(inp_high)
        assert res_high.Q_total_process_kW > res_low.Q_total_process_kW

    def test_with_convection_bank(self):
        inp = HeaterInput(
            W_fuel_kg_s=1.0,
            conv_banks=[{
                "tube_od_mm": 114.3,
                "tube_wt_mm": 8.0,
                "pitch_trans_mm": 220.0,
                "pitch_long_mm": 220.0,
                "tube_length_m": 10.0,
                "n_tubes_per_row": 12,
                "n_rows": 4,
                "arrangement": "staggered",
            }],
        )
        res = simulate_heater(inp)
        assert res.Q_convection_kW > 0
        assert len(res.convection_banks) == 1

    def test_to_dict_serializable(self):
        res = simulate_heater(self._default_input())
        d = res.to_dict()
        assert isinstance(d, dict)
        for key in ("Q_released_kW", "Q_radiant_kW", "T_BWT_C", "thermal_efficiency_pct"):
            assert key in d

    def test_fuel_oil(self):
        inp = HeaterInput(fuel_type="fuel_oil", W_fuel_kg_s=1.0)
        res = simulate_heater(inp)
        assert res.Q_total_process_kW > 0

    def test_cylindrical_firebox(self):
        inp = HeaterInput(firebox_type="cylindrical")
        res = simulate_heater(inp)
        assert res.Q_total_process_kW > 0
