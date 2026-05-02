"""
Tests for the combustion calculation module.
"""
import math
import pytest

from src.combustion import (
    combustion_calc,
    flue_gas_cp,
    flue_gas_enthalpy_approx,
    adiabatic_flame_temp,
    available_heat_fraction,
    FUEL_DB,
)


class TestCombustionCalc:
    def test_known_fuel_returns_result(self):
        res = combustion_calc("natural_gas", 0.15)
        assert res.fuel_type == "natural_gas"
        assert res.LHV == pytest.approx(47_141.0, rel=0.01)

    def test_actual_air_greater_than_stoic(self):
        res = combustion_calc("natural_gas", 0.20)
        assert res.actual_air > res.stoic_air

    def test_excess_air_zero_equals_stoic(self):
        res = combustion_calc("natural_gas", 0.0)
        assert res.actual_air == pytest.approx(res.stoic_air, rel=1e-6)

    def test_flue_gas_per_fuel_mass_balance(self):
        res = combustion_calc("natural_gas", 0.15)
        # flue = fuel + air
        assert res.flue_gas_per_fuel == pytest.approx(1.0 + res.actual_air, rel=1e-6)

    def test_flue_composition_sums_to_one(self):
        for fuel_type in FUEL_DB:
            res = combustion_calc(fuel_type, 0.15)
            total = sum(res.flue_comp_mass.values())
            assert total == pytest.approx(1.0, abs=0.02), (
                f"Fuel {fuel_type}: flue comp sum = {total:.4f}"
            )

    def test_unknown_fuel_raises(self):
        with pytest.raises(KeyError):
            combustion_calc("dragon_gas")


class TestFlueGasProperties:
    def test_cp_increases_with_temperature(self):
        cp_low = flue_gas_cp(400.0)
        cp_high = flue_gas_cp(1200.0)
        assert cp_high > cp_low

    def test_cp_reasonable_range(self):
        cp = flue_gas_cp(800.0)
        assert 1.0 < cp < 1.4, f"cp = {cp:.4f} kJ/(kg·K) out of typical range"

    def test_enthalpy_zero_at_reference(self):
        h = flue_gas_enthalpy_approx(288.15, T_ref_K=288.15)
        assert h == pytest.approx(0.0, abs=0.1)

    def test_enthalpy_positive_above_reference(self):
        h = flue_gas_enthalpy_approx(800.0, T_ref_K=288.15)
        assert h > 0

    def test_enthalpy_increases_monotonically(self):
        temps = [300, 500, 700, 1000, 1500]
        hs = [flue_gas_enthalpy_approx(T) for T in temps]
        assert all(hs[i] < hs[i + 1] for i in range(len(hs) - 1))


class TestAdiabaticFlameTemp:
    def test_natural_gas_aft_reasonable(self):
        T_aft = adiabatic_flame_temp("natural_gas", 0.15, 288.15)
        # Typical AFT for natural gas at 15% excess air ≈ 1800–2000 °C (in K: ~2100–2300 K)
        assert 1700 < T_aft < 2400, f"AFT = {T_aft:.1f} K"

    def test_more_excess_air_lowers_aft(self):
        T_aft_low = adiabatic_flame_temp("natural_gas", 0.10)
        T_aft_high = adiabatic_flame_temp("natural_gas", 0.40)
        assert T_aft_low > T_aft_high

    def test_air_preheat_raises_aft(self):
        T_cold = adiabatic_flame_temp("natural_gas", 0.15, T_air_K=288.15)
        T_hot = adiabatic_flame_temp("natural_gas", 0.15, T_air_K=500.0)
        assert T_hot > T_cold


class TestAvailableHeat:
    def test_fraction_between_zero_and_one(self):
        eta = available_heat_fraction("natural_gas", 0.15, T_stack_K=450.0)
        assert 0.0 < eta < 1.0

    def test_lower_stack_temp_gives_higher_efficiency(self):
        eta_low = available_heat_fraction("natural_gas", 0.15, T_stack_K=400.0)
        eta_high = available_heat_fraction("natural_gas", 0.15, T_stack_K=600.0)
        assert eta_low > eta_high

    def test_more_excess_air_lowers_efficiency(self):
        eta_low = available_heat_fraction("natural_gas", 0.10)
        eta_high = available_heat_fraction("natural_gas", 0.40)
        assert eta_low > eta_high
