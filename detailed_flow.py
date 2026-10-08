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
    if not all(math.isfinite(value) for value in (pressure_inlet_pa, pressure_outlet_pa, pressure_drop_pa)):
        raise ValueError("Inlet and outlet pressures must be finite.")
    if pressure_drop_pa <= 0:
        raise ValueError("Inlet pressure must be greater than outlet pressure.")
    if not path_components:
        raise ValueError("Add at least one path component.")

    chemical = Chemical(fluid, T=temperature_k, P=(pressure_inlet_pa + pressure_outlet_pa) / 2) if isinstance(fluid, str) else fluid
    rho, mu = chemical.rho, chemical.mu
    if not rho or not mu or not math.isfinite(rho) or not math.isfinite(mu) or rho <= 0 or mu <= 0:
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

    restriction_types = {"Barb", "LBarb", "Orifice", "Compression fitting"}
    tube_types = {"Straight tube", "Bend tube"}

    def upstream_connection_id(index: int) -> float | None:
        """Find the nearest known bore immediately upstream of a component."""
        for component in reversed(path_components[:index]):
            kind = component["type"]
            if kind in tube_types:
                return float(component["id_m"])
            if kind in restriction_types:
                return float(component["outlet_id_m"])
            if kind in {"Valve Cv", "Valve Kv"} and component.get("outlet_id_m") is not None:
                return float(component["outlet_id_m"])
            # Cv/Kv describes valve loss; its optional outlet bore may be unknown.
        return None

    def downstream_connection_id(index: int) -> float | None:
        """Find the nearest known bore immediately downstream of a component."""
        for component in path_components[index + 1:]:
            kind = component["type"]
            if kind in tube_types:
                return float(component["id_m"])
            if kind in restriction_types:
                return float(component["inlet_id_m"])
            if kind in {"Valve Cv", "Valve Kv"} and component.get("inlet_id_m") is not None:
                return float(component["inlet_id_m"])
            # Cv/Kv describes valve loss; its optional inlet bore may be unknown.
        return None

    def transition_zeta(from_id: float, to_id: float, reynolds: float) -> float:
        """Loss coefficient for a contraction or expansion between two bores."""
        if to_id <= from_id:
            return functions.get_zeta_in(to_id, from_id, Re=max(reynolds, 1))
        return functions.OutletDrag_coefficient(from_id, to_id)

    def losses(flow: float) -> tuple[float, list[dict[str, float]]]:
        total = 0.0
        rows: list[dict[str, float]] = []

        for index, component in enumerate(path_components, 1):
            kind = component["type"]
            raw_name = component.get("name")
            name = raw_name.strip() if isinstance(raw_name, str) else ""
            name = name or f"{kind} {index}"
            count = max(int(component.get("count", 1)), 1)
            loss, reynolds = 0.0, 0.0

            if kind == "Straight tube":
                diameter = float(component["id_m"])
                length = float(component["length_m"])
                if diameter <= 0 or length < 0:
                    raise ValueError("Straight tube ID must be positive and length cannot be negative.")
                reynolds = functions.reynolds_number(flow, diameter, rho, mu)
                loss = functions.tube_loss (flow,rho, mu,  diameter, length)

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

            elif kind in {"Barb", "Orifice", "LBarb", "Compression fitting"}:
                inlet_id = float(component["inlet_id_m"])
                outlet_id = float(component["outlet_id_m"])
                length = float(component.get("restrictor_length_m", 0))
                if inlet_id <= 0 or outlet_id <= 0 or length < 0:
                    raise ValueError(f"{kind} inlet ID, outlet ID, and length must be valid.")
                mean_id = (inlet_id + outlet_id) / 2
                reynolds = functions.reynolds_number(flow, mean_id, rho, mu)
                upstream_id = upstream_connection_id(index - 1) or inlet_id
                downstream_id = downstream_connection_id(index - 1)
                velocity = flow / (math.pi * mean_id**2 / 4)
                zeta = transition_zeta(upstream_id, inlet_id, reynolds)
                zeta += functions.friction_factor(max(reynolds, 1e-12)) * length / mean_id
                if downstream_id is not None:
                    zeta += transition_zeta(outlet_id, downstream_id, reynolds)
                elif index == len(path_components):
                    # Final restrictor discharges into a large outlet volume.
                    zeta += 1.0
                if kind == "LBarb":
                    zeta += 1.15
                loss = zeta * rho * velocity**2 / 2

            elif kind in {"Valve Cv", "Valve Kv"}:
                coefficient = float(component["valve_coefficient"])
                if coefficient <= 0:
                    raise ValueError("Valve Cv/Kv must be greater than zero.")
                

                reference_id = (component.get("inlet_id_m")
                                or component.get("outlet_id_m")
                                or upstream_connection_id(index - 1)
                                or downstream_connection_id(index - 1)
                            )

                if reference_id is not None:
                    reynolds = functions.reynolds_number(flow, reference_id, rho, mu)
                specific_gravity = rho / 997.0
                if kind == "Valve Cv":
                    loss = count * specific_gravity * (flow * 15_850.323 / coefficient) ** 2 * 6_894.757
                else:
                    loss = count * specific_gravity * (flow * 3600 / coefficient) ** 2 * 100_000
            else:
                raise ValueError(f"Unknown component type: {kind}")

            if not math.isfinite(loss) or loss < 0:
                raise ValueError(f"{name}: calculated pressure loss must be finite and non-negative.")
            total += loss
            rows.append({"Component": name, "Type": kind, "ΔP (Pa)": loss, "Re": reynolds})
        if not math.isfinite(total):
            raise ValueError("Total pressure loss must be finite.")
        return total, rows

    # At zero flow, losses are zero. Do not evaluate the hydraulic helpers at
    # zero: their laminar friction factors contain 64 / Re.
    lower, upper = 0.0, 0.01
    max_flow_m3_s = 10.0
    upper_loss = losses(upper)[0]
    while upper_loss < pressure_drop_pa and upper < max_flow_m3_s:
        upper = min(upper * 2, max_flow_m3_s)
        upper_loss = losses(upper)[0]
    if upper_loss < pressure_drop_pa:
        raise ValueError(
            "Unable to bracket a flow solution within the search limit "
            f"({max_flow_m3_s:g} m³/s). Check path resistance and pressure inputs."
        )
    for _ in range(80):
        midpoint = (lower + upper) / 2
        if losses(midpoint)[0] < pressure_drop_pa:
            lower = midpoint
        else:
            upper = midpoint
    flow = (lower + upper) / 2
    total, breakdown = losses(flow)
    # A narrow flow interval alone does not guarantee pressure balance, e.g.
    # when a friction correlation has a discontinuity at a regime boundary.
    residual_pa = pressure_drop_pa - total
    tolerance_pa = max(1e-6, pressure_drop_pa * 1e-6)
    if abs(residual_pa) > tolerance_pa:
        raise ValueError(
            "Unable to converge to pressure balance: "
            f"available ΔP = {pressure_drop_pa:.6g} Pa, "
            f"calculated loss = {total:.6g} Pa, "
            f"residual = {residual_pa:.6g} Pa "
            f"(tolerance {tolerance_pa:.6g} Pa). "
            "Check component inputs and flow-regime correlations."
        )
    return {
        "flow_m3_s": flow, "flow_l_min": flow * 60_000, "total_loss_pa": total,
        "density_kg_m3": rho, "viscosity_pa_s": mu, "breakdown": breakdown,
    }
