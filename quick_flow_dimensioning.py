import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from thermo import Mixture

from detailed_flow import solve_detailed_flow

PROCESS_PARAMETERS_ROOT = Path(__file__).resolve().parent.parent / "Process_Parameters"
if not PROCESS_PARAMETERS_ROOT.is_dir():
    raise FileNotFoundError(f"Process_Parameters project was not found at: {PROCESS_PARAMETERS_ROOT}")
if str(PROCESS_PARAMETERS_ROOT) not in sys.path:
    sys.path.insert(0, str(PROCESS_PARAMETERS_ROOT))


# ✅ Get the current directory dynamically
BASE_DIR = Path(__file__).resolve().parent
PROJECT_PATH = BASE_DIR / "Process_Parameters"

import general_function as process_functions

st.set_page_config(page_title="Detailed Fluid Path Calculator", layout="wide")
st.title("🌊 Detailed Fluid Path Calculator")
st.write("Solve flow rate from tube friction and a selectable restriction on each tube section.")

st.header("1. Fluid and pressure")
c1, c2, c3, c4 = st.columns(4)
with c1:
    p_inlet_bar = st.number_input("Inlet pressure (bar)", min_value=0.0, value=5.0, step=0.1)
with c2:
    p_outlet_bar = st.number_input("Outlet pressure (bar)", min_value=0.0, value=1.0, step=0.1)
with c3:
    fluid_mode = st.selectbox("Fluid type", ["Pure liquid", "Two-liquid mixture"])
with c4:
    temperature_c = st.number_input("Temperature (°C)", value=25.0, step=1.0)

fluid_options = ["water", "acetone", "ethanol", "methanol", "toluene", "MTBE", "Nitrogen", "Hydrogen"]
if fluid_mode == "Pure liquid":
    fluid = st.selectbox("Fluid (thermo)", fluid_options)
else:
    c1, c2, c3 = st.columns(3)
    with c1:
        component_1 = st.selectbox("Component 1", fluid_options, key="component_1")
    with c2:
        component_2 = st.selectbox("Component 2", fluid_options, index=1, key="component_2")
    with c3:
        component_1_fraction = st.number_input("Component 1 mole fraction", min_value=0.001, max_value=0.999, value=0.5, step=0.01)
    if component_1 == component_2:
        st.error("Choose two different components for a mixture.")
        st.stop()
    fluid = Mixture(
        IDs=[component_1, component_2],
        zs=[component_1_fraction, 1 - component_1_fraction],
        T=temperature_c + 273.15,
        P=((p_inlet_bar + p_outlet_bar) / 2) * 100_000,
    )
st.caption("Density and viscosity come from thermo at the selected temperature and mean pressure.")

st.header("2. Tube sections and restrictions")
default_tubes = pd.DataFrame([
    {"Section": "Feed tube", "Tube ID (mm)": 6.35, "Length (m)": 2.0, "Restriction": "None", "Count": 1, "Restriction ID (mm)": 1.0, "Restriction Length (mm)": 10.0, "Valve Cv/Kv": 1.0, "Bend r/ID": 11.71},
    {"Section": "Process tube", "Tube ID (mm)": 6.00, "Length (m)": 1.0, "Restriction": "None", "Count": 1, "Restriction ID (mm)": 1.0, "Restriction Length (mm)": 10.0, "Valve Cv/Kv": 1.0, "Bend r/ID": 11.71},
])
tubes = st.data_editor(
    default_tubes, num_rows="dynamic", width="stretch",
    column_config={
        "Section": st.column_config.TextColumn(required=True),
        "Tube ID (mm)": st.column_config.NumberColumn(min_value=0.01, format="%.3f"),
        "Length (m)": st.column_config.NumberColumn(min_value=0.0, format="%.3f"),
        "Restriction": st.column_config.SelectboxColumn(options=["None", "Bend", "Barb", "Orifice", "Valve Cv", "Valve Kv"], required=True),
        "Count": st.column_config.NumberColumn(min_value=1, step=1),
        "Restriction ID (mm)": st.column_config.NumberColumn(min_value=0.01, format="%.3f"),
        "Restriction Length (mm)": st.column_config.NumberColumn(min_value=0.0, format="%.3f"),
        "Valve Cv/Kv": st.column_config.NumberColumn(min_value=0.01, format="%.3f"),
        "Bend r/ID": st.column_config.SelectboxColumn(options=[2.26, 3.04, 6.53, 11.71]),
    },
)
st.caption("Use only the relevant fields for the chosen restriction: bend r/ID for bends; restriction ID/length for barbs; restriction ID for orifices; and Valve Cv/Kv for valves.")

pressure_drop_pa = (p_inlet_bar - p_outlet_bar) * 100_000
if pressure_drop_pa <= 0:
    st.warning("Outlet pressure must be lower than inlet pressure.")
else:
    try:
        sections = [
            {
                "name": row["Section"], "id_m": row["Tube ID (mm)"] / 1000,
                "length_m": row["Length (m)"], "restriction": row["Restriction"],
                "restriction_count": row["Count"], "restriction_id_m": row["Restriction ID (mm)"] / 1000,
                "restriction_length_m": row["Restriction Length (mm)"] / 1000,
                "valve_coefficient": row["Valve Cv/Kv"], "bend_radius_ratio": row["Bend r/ID"],
            }
            for _, row in tubes.iterrows()
        ]
        result = solve_detailed_flow(
            pressure_inlet_pa=p_inlet_bar * 100_000,
            pressure_outlet_pa=p_outlet_bar * 100_000,
            temperature_k=temperature_c + 273.15,
            fluid=fluid,
            tube_sections=sections,
            functions=process_functions,
        )

        st.header("3. Results")
        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric("Solved flow rate", f"{result['flow_l_min']:.2f} L/min")
        m2.metric("Pressure drop", f"{pressure_drop_pa / 100_000:.2f} bar")
        m3.metric("Calculated losses", f"{result['total_loss_pa'] / 100_000:.2f} bar")
        m4.metric("Density", f"{result['density_kg_m3']:.1f} kg/m³")
        m5.metric("Viscosity", f"{result['viscosity_pa_s'] * 1000:.3f} mPa·s")
        breakdown = pd.DataFrame(result["breakdown"])
        breakdown["ΔP (bar)"] = breakdown["Total ΔP (Pa)"] / 100_000
        st.subheader("Pressure-loss breakdown")
        st.dataframe(breakdown[["Component", "Restriction", "ΔP (bar)", "Re"]], width="stretch", hide_index=True)
        fig, ax = plt.subplots(figsize=(9, 3.5))
        ax.bar(breakdown["Component"], breakdown["ΔP (bar)"], color="#2b5c8f")
        ax.set_ylabel("Pressure loss (bar)")
        ax.set_title("Loss by tube section")
        ax.tick_params(axis="x", rotation=20)
        st.pyplot(fig)
    except (ValueError, ZeroDivisionError, TypeError) as error:
        st.error(f"Unable to calculate this path: {error}")
