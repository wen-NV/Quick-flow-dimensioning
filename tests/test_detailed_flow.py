"""Regression checks for pressure-balance acceptance and solver failures."""
import math
import unittest
from types import SimpleNamespace

import general_function as functions
from detailed_flow import solve_detailed_flow


class FlowSolverTests(unittest.TestCase):
    def solve(self, path, pressure_drop=1000.0, fluid=None):
        return solve_detailed_flow(
            pressure_inlet_pa=100000.0 + pressure_drop,
            pressure_outlet_pa=100000.0,
            temperature_k=298.15,
            fluid=fluid if fluid is not None else SimpleNamespace(rho=1000.0, mu=0.001),
            path_components=path,
            functions=functions,
        )

    def tube(self, length=1.0, diameter=0.004):
        return {"type": "Straight tube", "id_m": diameter, "length_m": length}

    def test_laminar_tube_matches_poiseuille(self):
        result = self.solve([self.tube()], pressure_drop=10.0)
        expected_flow = math.pi * 10.0 * 0.004**4 / (128 * 0.001)
        self.assertAlmostEqual(result["flow_m3_s"] / expected_flow, 1.0, places=9)
        self.assertAlmostEqual(result["total_loss_pa"], 10.0, places=7)

    def test_kv_valve_matches_liquid_flow_relation(self):
        result = self.solve([{"type": "Valve Kv", "valve_coefficient": 1.5}])
        expected_flow = 1.5 / 3600 * math.sqrt(1000 / 100000 / (1000 / 997))
        self.assertAlmostEqual(result["flow_m3_s"] / expected_flow, 1.0, places=9)

    def test_default_path_balances_pressure(self):
        path = [
            self.tube(length=0.025),
            {"type": "Barb", "inlet_id_m": 0.0024, "outlet_id_m": 0.0024,
             "restrictor_length_m": 0.034},
            {"type": "Bend tube", "id_m": 0.003175, "length_m": 0.04,
             "bend_radius_m": 0.05, "bend_angle_deg": 90.0},
            self.tube(length=0.025),
        ]
        result = self.solve(path, pressure_drop=80000)
        self.assertGreater(result["flow_m3_s"], 0)
        self.assertAlmostEqual(result["total_loss_pa"], 80000, places=5)
        self.assertAlmostEqual(sum(row["ΔP (Pa)"] for row in result["breakdown"]),
                               result["total_loss_pa"], places=5)

    def test_zero_resistance_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unable to bracket"):
            self.solve([self.tube(length=0)])

    def test_solution_above_search_limit_rejected(self):
        with self.assertRaisesRegex(ValueError, "Unable to bracket"):
            self.solve([{"type": "Valve Kv", "valve_coefficient": 1e10}])

    def test_discontinuous_friction_does_not_produce_fake_solution(self):
        with self.assertRaisesRegex(ValueError, "Unable to converge to pressure balance"):
            self.solve([self.tube()], pressure_drop=1200)

    def test_nonfinite_loss_rejected(self):
        with self.assertRaisesRegex(ValueError, "pressure loss must be finite"):
            self.solve([self.tube(diameter=float("nan"))])

    def test_nonfinite_pressure_rejected(self):
        with self.assertRaisesRegex(ValueError, "pressures must be finite"):
            self.solve([self.tube()], pressure_drop=float("inf"))

    def test_invalid_fluid_properties_rejected(self):
        for rho, mu in [(float("nan"), 0.001), (1000, float("inf")), (-1000, 0.001)]:
            with self.subTest(rho=rho, mu=mu):
                with self.assertRaises(ValueError):
                    self.solve([self.tube()], fluid=SimpleNamespace(rho=rho, mu=mu))


if __name__ == "__main__":
    unittest.main()
