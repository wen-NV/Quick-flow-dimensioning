"""Detailed fluid-path solver using Process_Parameters hydraulic helpers."""

from __future__ import annotations

import math
from typing import Any

from thermo import Chemical


def solve_detailed_flow(
    *, pressure_inlet_pa: float, pressure_outlet_pa: float, temperature_k: float,
    fluid: str | Chemical, tube_sections: list[dict[str, Any]], functions: Any,
) -> dict[str, Any]:
    """Solve one path made of tube rows and optional per-row restrictions.

    Each section needs ``id_m`` and ``length_m``.  Its optional ``restriction``
    is one of: None, Bend, Barb, Orifice, Valve Cv, or Valve Kv.
    """
    pressure_drop_pa = pressure_inlet_pa - pressure_outlet_pa
    if pressure_drop_pa <= 0:
        raise ValueError("Inlet pressure must be greater than outlet pressure.")
    if not tube_sections:
        raise ValueError("Add at least one tube section.")
    if any(s["id_m"] <= 0 or s["length_m"] < 0 for s in tube_sections):
        raise ValueError("Tube IDs must be positive and lengths cannot be negative.")

    chemical = Chemical(fluid, T=temperature_k, P=(pressure_inlet_pa + pressure_outlet_pa) / 2) if isinstance(fluid, str) else fluid
    rho, mu = chemical.rho, chemical.mu
    if not rho or not mu:
        raise ValueError(f"Thermo could not determine density or viscosity for {fluid!r}.")

    def local_loss(section: dict[str, Any], flow: float, velocity: float, reynolds: float) -> float:
        kind = section.get("restriction", "None")
        count = max(int(section.get("restriction_count", 1)), 1)
        diameter = section["id_m"]
        dynamic_pressure = rho * velocity**2 / 2

        if kind == "None":
            return 0.0
        if kind == "Bend":
            ratio = float(section.get("bend_radius_ratio", 11.71))
            return count * functions.get_zeta_bend(max(reynolds, 10), ratio) * dynamic_pressure
        if kind in {"Barb", "Orifice"}:
            restriction_id = float(section.get("restriction_id_m", 0))
            if restriction_id <= 0 or restriction_id > diameter:
                raise ValueError(f"{kind} ID must be greater than zero and no larger than its tube ID.")
            restriction_velocity = flow / (math.pi * restriction_id**2 / 4)
            restriction_re = functions.reynolds_number(flow, restriction_id, rho, mu)
            zeta = functions.get_zeta_in(restriction_id, diameter, Re=max(restriction_re, 1))
            zeta += functions.OutletDrag_coefficient(restriction_id, diameter)
            if kind == "Barb":
                restriction_length = float(section.get("restriction_length_m", 0))
                zeta += functions.friction_factor(max(restriction_re, 1e-12)) * restriction_length / restriction_id
            return count * zeta * rho * restriction_velocity**2 / 2
        if kind in {"Valve Cv", "Valve Kv"}:
            coefficient = float(section.get("valve_coefficient", 0))
            if coefficient <= 0:
                raise ValueError("Valve Cv/Kv must be greater than zero.")
            specific_gravity = rho / 997.0
            if kind == "Valve Cv":
                flow_gpm = flow * 15_850.323
                return count * specific_gravity * (flow_gpm / coefficient) ** 2 * 6_894.757
            flow_m3_h = flow * 3600
            return count * specific_gravity * (flow_m3_h / coefficient) ** 2 * 100_000
        raise ValueError(f"Unknown restriction type: {kind}")

    def losses(flow: float) -> tuple[float, list[dict[str, float]]]:
        rows: list[dict[str, float]] = []
        total = 0.0
        for index, section in enumerate(tube_sections, 1):
            diameter = section["id_m"]
            velocity = flow / (math.pi * diameter**2 / 4)
            reynolds = functions.reynolds_number(flow, diameter, rho, mu)
            friction = functions.friction_factor(max(reynolds, 1e-12))
            tube_loss = friction * section["length_m"] / diameter * rho * velocity**2 / 2
            restriction_loss = local_loss(section, flow, velocity, reynolds)
            total += tube_loss + restriction_loss
            name = section.get("name") or f"Tube {index}"
            kind = section.get("restriction", "None")
            rows.append({
                "Component": name, "Restriction": kind, "Tube ΔP (Pa)": tube_loss,
                "Restriction ΔP (Pa)": restriction_loss, "Total ΔP (Pa)": tube_loss + restriction_loss,
                "Re": reynolds,
            })
        return total, rows

    lower, upper = 0.0, 0.01
    while losses(upper)[0] < pressure_drop_pa and upper < 10:
        upper *= 2
    for _ in range(80):
        midpoint = (lower + upper) / 2
        if losses(midpoint)[0] < pressure_drop_pa:
            lower = midpoint
        else:
            upper = midpoint
    flow = (lower + upper) / 2
    total, breakdown = losses(flow)
    return {
        "flow_m3_s": flow,
        "flow_l_min": flow * 60_000,
        "total_loss_pa": total,
        "density_kg_m3": rho,
        "viscosity_pa_s": mu,
        "breakdown": breakdown,
    }
