"""
Tests for the Flask web application routes.
"""
import json
import pytest

from app import app as flask_app


@pytest.fixture
def client():
    flask_app.config["TESTING"] = True
    with flask_app.test_client() as c:
        yield c


class TestIndexPage:
    def test_index_returns_200(self, client):
        rv = client.get("/")
        assert rv.status_code == 200

    def test_index_contains_form(self, client):
        rv = client.get("/")
        html = rv.data.decode()
        assert "heater-form" in html or 'method="POST"' in html

    def test_index_shows_fuel_options(self, client):
        rv = client.get("/")
        html = rv.data.decode()
        assert "natural_gas" in html


class TestCalculateEndpoint:
    def _default_form(self):
        return {
            "case_name": "Test Case",
            "heater_service": "Test Heater",
            "fuel_type": "natural_gas",
            "excess_air_pct": "15",
            "W_fuel_kg_s": "1.0",
            "T_air_C": "15",
            "firebox_type": "box",
            "firebox_length_m": "10",
            "firebox_width_m": "4",
            "firebox_height_m": "6",
            "n_tubes_radiant": "24",
            "rad_tube_od_mm": "168.3",
            "rad_tube_wt_mm": "8.0",
            "rad_tube_pitch_mm": "336.6",
            "n_radiant_rows": "1",
            "rad_tube_length_m": "10",
            "rad_tube_material": "Cr5Mo",
            "process_fluid": "crude_oil",
            "W_process_kg_s": "50",
            "T_proc_in_C": "100",
            "loss_fraction_pct": "2",
            "peak_flux_factor": "1.8",
        }

    def test_form_post_returns_200(self, client):
        rv = client.post("/calculate", data=self._default_form())
        assert rv.status_code == 200

    def test_form_post_shows_results(self, client):
        rv = client.post("/calculate", data=self._default_form())
        html = rv.data.decode()
        assert "BWT" in html or "T_BWT" in html or "过桥" in html

    def test_json_api_returns_result(self, client):
        payload = dict(self._default_form())
        rv = client.post(
            "/calculate",
            data=json.dumps(payload),
            content_type="application/json",
        )
        assert rv.status_code == 200
        data = rv.get_json()
        assert "Q_released_kW" in data
        assert "T_BWT_C" in data

    def test_json_api_with_convection(self, client):
        payload = dict(self._default_form())
        payload["conv_banks"] = [{
            "tube_od_mm": 114.3,
            "tube_wt_mm": 8.0,
            "pitch_trans_mm": 220.0,
            "pitch_long_mm": 220.0,
            "tube_length_m": 10.0,
            "n_tubes_per_row": 12,
            "n_rows": 4,
            "arrangement": "staggered",
        }]
        rv = client.post(
            "/calculate",
            data=json.dumps(payload),
            content_type="application/json",
        )
        assert rv.status_code == 200
        data = rv.get_json()
        assert data.get("Q_convection_kW", 0) > 0


class TestAPIEndpoints:
    def test_fuels_endpoint(self, client):
        rv = client.get("/api/fuels")
        assert rv.status_code == 200
        data = rv.get_json()
        assert "natural_gas" in data

    def test_fluids_endpoint(self, client):
        rv = client.get("/api/fluids")
        assert rv.status_code == 200
        data = rv.get_json()
        assert "crude_oil" in data

    def test_defaults_endpoint(self, client):
        rv = client.get("/api/defaults")
        assert rv.status_code == 200
        data = rv.get_json()
        assert "fuel_type" in data
        assert "W_fuel_kg_s" in data
