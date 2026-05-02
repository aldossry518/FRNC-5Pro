"""
FRNC-5Pro  –  管式加热炉传热计算软件  (Web Application)
==========================================================
Flask web application providing the front-end for the fired-heater
heat transfer calculation engine.
"""
from __future__ import annotations

import json
import math
from flask import Flask, render_template, request, jsonify

from src.heater import HeaterInput, simulate_heater
from src.combustion import FUEL_DB
from src.properties import PROCESS_FLUID_PRESETS, TUBE_MATERIALS
from src.geometry import FireboxGeometry

app = Flask(__name__)
app.config["JSON_SORT_KEYS"] = False


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

@app.route("/")
def index():
    """Landing page – input form."""
    return render_template(
        "index.html",
        fuels=FUEL_DB,
        fluids=PROCESS_FLUID_PRESETS,
        materials=TUBE_MATERIALS,
    )


@app.route("/calculate", methods=["POST"])
def calculate():
    """Run calculation from form data or JSON payload."""
    # Accept both form POST and JSON API calls
    if request.is_json:
        data = request.get_json(force=True)
    else:
        data = request.form.to_dict()

    try:
        inp = _parse_input(data)
        result = simulate_heater(inp)
        result_dict = result.to_dict()

        if request.is_json:
            return jsonify(result_dict)

        return render_template(
            "results.html",
            result=result,
            result_dict=result_dict,
            inp=inp,
        )
    except Exception as exc:
        error_msg = str(exc)
        if request.is_json:
            return jsonify({"error": error_msg}), 400
        return render_template("index.html",
                               error=error_msg,
                               fuels=FUEL_DB,
                               fluids=PROCESS_FLUID_PRESETS,
                               materials=TUBE_MATERIALS), 400


@app.route("/api/fuels")
def api_fuels():
    """Return fuel database."""
    return jsonify({k: {"name": v["name"], "LHV": v["LHV"]} for k, v in FUEL_DB.items()})


@app.route("/api/fluids")
def api_fluids():
    """Return process fluid database."""
    return jsonify({k: {"name": v["name"], "cp_avg": v["cp_avg"]}
                    for k, v in PROCESS_FLUID_PRESETS.items()})


@app.route("/api/defaults")
def api_defaults():
    """Return default input parameters."""
    return jsonify(_default_input_dict())


# ---------------------------------------------------------------------------
# Input parsing helper
# ---------------------------------------------------------------------------

def _f(data: dict, key: str, default):
    """Safely parse a float from form data."""
    val = data.get(key, default)
    if val is None or val == "":
        return float(default)
    try:
        return float(val)
    except (ValueError, TypeError):
        return float(default)


def _i(data: dict, key: str, default: int) -> int:
    """Safely parse an int from form data."""
    try:
        return int(_f(data, key, default))
    except Exception:
        return default


def _parse_input(data: dict) -> HeaterInput:
    """Convert form/JSON data dict into a HeaterInput object."""
    # Convection banks from JSON array or form fields
    conv_banks = []

    # JSON input may pass conv_banks directly
    if "conv_banks" in data and isinstance(data["conv_banks"], list):
        conv_banks = data["conv_banks"]
    else:
        # Build from numbered form fields conv_bank_0_*, conv_bank_1_*, …
        idx = 0
        while True:
            prefix = f"conv_bank_{idx}_"
            if not any(k.startswith(prefix) for k in data):
                break
            bank = {
                "tube_od_mm":      _f(data, prefix + "tube_od_mm", 114.3),
                "tube_wt_mm":      _f(data, prefix + "tube_wt_mm", 8.0),
                "pitch_trans_mm":  _f(data, prefix + "pitch_trans_mm", 220.0),
                "pitch_long_mm":   _f(data, prefix + "pitch_long_mm", 220.0),
                "tube_length_m":   _f(data, prefix + "tube_length_m", 10.0),
                "n_tubes_per_row": _i(data, prefix + "n_tubes_per_row", 12),
                "n_rows":          _i(data, prefix + "n_rows", 4),
                "arrangement":     data.get(prefix + "arrangement", "staggered"),
                "tube_material":   data.get(prefix + "tube_material", "carbon_steel"),
            }
            fin_h = _f(data, prefix + "fin_height_mm", 0)
            if fin_h > 0:
                bank["fin_height_mm"] = fin_h
                bank["fin_thickness_mm"] = _f(data, prefix + "fin_thickness_mm", 2.5)
                bank["fin_pitch_per_m"] = _f(data, prefix + "fin_pitch_per_m", 197.0)
                bank["fin_material"] = data.get(prefix + "fin_material", "carbon_steel")
            conv_banks.append(bank)
            idx += 1

    T_proc_out = None
    raw_out = data.get("T_proc_out_C", "")
    if raw_out and str(raw_out).strip():
        try:
            T_proc_out = float(raw_out)
        except ValueError:
            pass

    return HeaterInput(
        case_name=data.get("case_name", "Case 1"),
        heater_service=data.get("heater_service", "Process Heater"),
        fuel_type=data.get("fuel_type", "natural_gas"),
        excess_air_pct=_f(data, "excess_air_pct", 15.0),
        W_fuel_kg_s=_f(data, "W_fuel_kg_s", 1.0),
        T_air_C=_f(data, "T_air_C", 15.0),
        firebox_type=data.get("firebox_type", "box"),
        firebox_length_m=_f(data, "firebox_length_m", 10.0),
        firebox_width_m=_f(data, "firebox_width_m", 4.0),
        firebox_height_m=_f(data, "firebox_height_m", 6.0),
        n_tubes_radiant=_i(data, "n_tubes_radiant", 24),
        rad_tube_od_mm=_f(data, "rad_tube_od_mm", 168.3),
        rad_tube_wt_mm=_f(data, "rad_tube_wt_mm", 8.0),
        rad_tube_pitch_mm=_f(data, "rad_tube_pitch_mm", 336.6),
        n_radiant_rows=_i(data, "n_radiant_rows", 1),
        rad_tube_length_m=_f(data, "rad_tube_length_m", 10.0),
        rad_tube_material=data.get("rad_tube_material", "Cr5Mo"),
        process_fluid=data.get("process_fluid", "crude_oil"),
        W_process_kg_s=_f(data, "W_process_kg_s", 50.0),
        T_proc_in_C=_f(data, "T_proc_in_C", 100.0),
        T_proc_out_C=T_proc_out,
        conv_banks=conv_banks,
        loss_fraction=_f(data, "loss_fraction_pct", 2.0) / 100.0,
        peak_flux_factor=_f(data, "peak_flux_factor", 1.8),
    )


def _default_input_dict() -> dict:
    return {
        "case_name": "Case 1",
        "heater_service": "原油加热炉 / Crude Preheater",
        "fuel_type": "natural_gas",
        "excess_air_pct": 15.0,
        "W_fuel_kg_s": 1.0,
        "T_air_C": 15.0,
        "firebox_type": "box",
        "firebox_length_m": 10.0,
        "firebox_width_m": 4.0,
        "firebox_height_m": 6.0,
        "n_tubes_radiant": 24,
        "rad_tube_od_mm": 168.3,
        "rad_tube_wt_mm": 8.0,
        "rad_tube_pitch_mm": 336.6,
        "n_radiant_rows": 1,
        "rad_tube_length_m": 10.0,
        "rad_tube_material": "Cr5Mo",
        "process_fluid": "crude_oil",
        "W_process_kg_s": 50.0,
        "T_proc_in_C": 100.0,
        "loss_fraction_pct": 2.0,
        "peak_flux_factor": 1.8,
    }


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
