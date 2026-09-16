"""Ordered fluid-path solver using Process_Parameters hydraulic helpers."""

from __future__ import annotations

import math
from typing import Any

from thermo import Chemical


def solve_detailed_flow(
    *, pressure_inlet_pa: float, pressure_outlet_pa: float, temperature_k: float,
    fluid: str | Chemical, path_components: list[dict[str, Any]], functions: Any,
) -> dict[str, Any]:
    """Solve an ordered path of straight tubes and separate restrictions.

    A restriction (barb or orifice) sits between tube components. Its own ID
    is checked against both the preceding and following tube IDs.
    """
    pressure_drop_pa = pressure_inlet_pa - pressure_outlet_pa
    if pressure_drop_pa <= 0:
        raise ValueError("Inlet pressure must be greater than outlet pressure.")
    if not path_components:
        raise ValueError("Add at least one path component.")

    chemical = Chemical(fluid, T=temperature_k, P=(pressure_inlet_pa + pressure_outlet_pa) / 2) if isinstance(fluid, str) else fluid
    rho, mu = chemical.rho, chemical.mu
    if not rho or not mu:
        raise ValueError(f"Thermo could not determine density or viscosity for {fluid!r}.")

    def bend_zeta(reynolds: float, radius_ratio: float, angle_deg: float) -> float:
        """Interpolate available 90° bend data, then scale for bend angle."""
        if radius_ratio <= 0 or angle_deg <= 0:
            raise ValueError("Bend radius and angle must be positive.")
        ratios = [2.26, 3.04, 6.53, 11.71]
        zetas = [functions.get_zeta_bend(max(reynolds, 10), ratio) for ratio in ratios]
        if radius_ratio <= ratios[0]:
            zeta_90 = zetas[0]
        elif radius_ratio >= ratios[-1]:
            zeta_90 = zetas[-1]
        else:
            for lower, upper, zeta_lower, zeta_upper in zip(ratios, ratios[1:], zetas, zetas[1:]):
                if lower <= radius_ratio <= upper:
                    fraction = (radius_ratio - lower) / (upper - lower)
                    zeta_90 = zeta_lower + fraction * (zeta_upper - zeta_lower)
                    break
        return zeta_90 * angle_deg / 90

    def losses(flow: float) -> tuple[float, list[dict[str, float]]]:
        total = 0.0
        rows: list[dict[str, float]] = []
        previous_tube_id: float | None = None
        next_tube_ids: list[float | None] = [None] * len(path_components)
        next_id: float | None = None
        for index in range(len(path_components) - 1, -1, -1):
            next_tube_ids[index] = next_id
            if path_components[index]["type"] in {"Straight tube", "Bend tube"}:
                next_id = float(path_components[index]["id_m"])

        for index, component in enumerate(path_components, 1):
            kind = component["type"]
            name = component.get("name") or f"{kind} {index}"
            count = max(int(component.get("count", 1)), 1)
            loss, reynolds = 0.0, 0.0

            if kind == "Straight tube":
                diameter = float(component["id_m"])
                length = float(component["length_m"])
                if diameter <= 0 or length < 0:
                    raise ValueError("Straight tube ID must be positive and length cannot be negative.")
                
                loss = functions.tube_loss (flow,rho, mu,  diameter, length)
                previous_tube_id = diameter

            elif kind == "Bend tube":
                diameter = float(component["id_m"])
                if diameter <= 0:
                    raise ValueError("Bend ID must be positive.")
                velocity = flow / (math.pi * diameter**2 / 4)
                reynolds = functions.reynolds_number(flow, diameter, rho, mu)
                radius = float(component.get("bend_radius_m", 0))
                angle = float(component.get("bend_angle_deg", 90))
                ratio = radius / diameter
                arc_length = radius * math.radians(angle)
                friction = functions.friction_factor(max(reynolds, 1e-12))
                loss = count * (
                    friction * arc_length / diameter + bend_zeta(reynolds, ratio, angle) # need to be check for how many bends, the ratio will change, according to VDI Heat Atlas
                ) * rho * velocity**2 / 2
                previous_tube_id = diameter

            elif kind in {"Barb", "Orifice"}:
                if previous_tube_id is None:
                    raise ValueError(f"Add a straight or bend tube before the {kind.lower()} row.")
                next_tube_id = next_tube_ids[index - 1]
                if next_tube_id is None:
                    raise ValueError(f"Add a straight or bend tube after the {kind.lower()} row.")
                inlet_id = float(component["inlet_id_m"])
                outlet_id = float(component["outlet_id_m"])
                length = float(component.get("restrictor_length_m", 0))
                if inlet_id <= 0 or outlet_id <= 0 or length < 0:
                    raise ValueError(f"{kind} inlet ID, outlet ID, and length must be valid.")
                if inlet_id > previous_tube_id or outlet_id > next_tube_id:
                    raise ValueError(f"{kind} inlet ID must fit its preceding tube and outlet ID must fit its following tube.")
                #mean_id = (inlet_id + outlet_id) / 2
                #velocity = flow / (math.pi * mean_id**2 / 4)
                #reynolds = functions.reynolds_number(flow, mean_id, rho, mu)
                #zeta = functions.get_zeta_in(inlet_id, previous_tube_id, Re=max(reynolds, 1))
                #zeta += functions.OutletDrag_coefficient(outlet_id, next_tube_id)
                #zeta += functions.friction_factor(max(reynolds, 1e-12)) * length / mean_id
                #loss = count * zeta * rho * velocity**2 / 2
                loss = functions.barb_dp (flow,mu, rho, inlet_id,outlet_id, previous_tube_id, next_tube_id, length)


            elif kind in {"LBarb"}:
                            if previous_tube_id is None:
                                raise ValueError(f"Add a straight or bend tube before the {kind.lower()} row.")
                            next_tube_id = next_tube_ids[index - 1]
                            if next_tube_id is None:
                                raise ValueError(f"Add a straight or bend tube after the {kind.lower()} row.")
                            inlet_id = float(component["inlet_id_m"])
                            outlet_id = float(component["outlet_id_m"])
                            length = float(component.get("restrictor_length_m", 0))
                            if inlet_id <= 0 or outlet_id <= 0 or length < 0:
                                raise ValueError(f"{kind} inlet ID, outlet ID, and length must be valid.")
                            if inlet_id > previous_tube_id or outlet_id > next_tube_id:
                                raise ValueError(f"{kind} inlet ID must fit its preceding tube and outlet ID must fit its following tube.")
                            #mean_id = (inlet_id + outlet_id) / 2
                            #velocity = flow / (math.pi * mean_id**2 / 4)
                            #reynolds = functions.reynolds_number(flow, mean_id, rho, mu)
                            #zeta = functions.get_zeta_in(inlet_id, previous_tube_id, Re=max(reynolds, 1))
                            #zeta += functions.OutletDrag_coefficient(outlet_id, next_tube_id)
                            #zeta += functions.friction_factor(max(reynolds, 1e-12)) * length / mean_id
                            #loss = count * zeta * rho * velocity**2 / 2
                            loss = functions.Lbarb_dp (flow,mu, rho, inlet_id,outlet_id, previous_tube_id, next_tube_id, length)    

            elif kind in {"Valve Cv", "Valve Kv"}:
                coefficient = float(component["valve_coefficient"])
                if coefficient <= 0:
                    raise ValueError("Valve Cv/Kv must be greater than zero.")
                specific_gravity = rho / 997.0
                if kind == "Valve Cv":
                    loss = count * specific_gravity * (flow * 15_850.323 / coefficient) ** 2 * 6_894.757
                else:
                    loss = count * specific_gravity * (flow * 3600 / coefficient) ** 2 * 100_000
            else:
                raise ValueError(f"Unknown component type: {kind}")

            total += loss
            rows.append({"Component": name, "Type": kind, "ΔP (Pa)": loss, "Re": reynolds})
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
        "flow_m3_s": flow, "flow_l_min": flow * 60_000, "total_loss_pa": total,
        "density_kg_m3": rho, "viscosity_pa_s": mu, "breakdown": breakdown,
    }
