import sys
import importlib
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from thermo import Mixture

import detailed_flow

# Streamlit retains imported modules between reruns; reload this local solver
# so edits to detailed_flow.py are applied without restarting the app.
detailed_flow = importlib.reload(detailed_flow)


# ✅ Get the current directory dynamically
#BASE_DIR = Path(__file__).resolve().parent
#PROJECT_PATH = BASE_DIR / "Process_Parameters"



import os
import urllib.request

RAW_URL = "https://raw.githubusercontent.com/wen-NV/general_fluid_function/main/general_function.py"

# Download and overwrite 'general_function.py' to ensure you always have the latest updates
urllib.request.urlretrieve(RAW_URL, "general_function.py")

# Now import directly
import general_function as  process_functions

st.set_page_config(page_title="Detailed Fluid Path Calculator", layout="wide")
st.title("🌊 Detailed Fluid Path Calculator")
st.write("Solve flow rate from tube friction and a selectable restriction on each tube section.")

st.header("1. Fluid and pressure")
c1, c2, c3, c4 = st.columns(4)
with c1:
    p_inlet_bar = st.number_input("Inlet pressure (bar)", min_value=0.0, value=1.8, step=0.1)
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

st.header("2. Ordered flow path")
default_path = pd.DataFrame([
    {"Name": "Feed tube from vessel", "Type": "Straight tube", "Tube ID (mm)": 4, "Tube Length (mm)": 25.0, "Restrictor Inlet ID (mm)": 2.4, "Restrictor Outlet ID (mm)": 2.4, "Restrictor Length (mm)": 34.0, "Count": 1, "Valve Cv/Kv": 1.0, "Bend Radius (mm)": 50.0, "Bend Angle (°)": 90.0},
    {"Name": "Barb", "Type": "Barb", "Tube ID (mm)": 0.0, "Tube Length (mm)": 0.0, "Restrictor Inlet ID (mm)": 2.40, "Restrictor Outlet ID (mm)": 2.40, "Restrictor Length (mm)": 34.0, "Count": 1, "Valve Cv/Kv": 1.0, "Bend Radius (mm)": 50.0, "Bend Angle (°)": 90.0},
    {"Name": "inter tube", "Type": "Bend tube", "Tube ID (mm)": 3.175, "Tube Length (mm)": 40.0, "Restrictor Inlet ID (mm)": 2.4, "Restrictor Outlet ID (mm)": 2.4, "Restrictor Length (mm)": 34.0, "Count": 1, "Valve Cv/Kv": 1.0, "Bend Radius (mm)": 50.0, "Bend Angle (°)": 90.0},
    {"Name": "Barb", "Type": "Barb", "Tube ID (mm)": 0.0, "Tube Length (mm)": 0.0, "Restrictor Inlet ID (mm)": 2.40, "Restrictor Outlet ID (mm)": 2.40, "Restrictor Length (mm)": 34.0, "Count": 1, "Valve Cv/Kv": 1.0, "Bend Radius (mm)": 50.0, "Bend Angle (°)": 90.0},
    {"Name": "Outlet tube to vessel", "Type": "Straight tube", "Tube ID (mm)": 4.00, "Tube Length (mm)": 25.0, "Restrictor Inlet ID (mm)": 2.4, "Restrictor Outlet ID (mm)": 2.4, "Restrictor Length (mm)": 10.0, "Count": 1, "Valve Cv/Kv": 1.0, "Bend Radius (mm)": 50.0, "Bend Angle (°)": 90.0},

])
path = st.data_editor(
    default_path, num_rows="dynamic", width="stretch",
    column_config={
        "Name": st.column_config.TextColumn(required=True),
        "Type": st.column_config.SelectboxColumn(options=["Straight tube", "Bend tube", "Barb", "Orifice", "LBarb","Valve Cv", "Valve Kv"], required=True),
        "Tube ID (mm)": st.column_config.NumberColumn(min_value=0.01, format="%.3f"),
        "Tube Length (mm)": st.column_config.NumberColumn(min_value=0.0, format="%.3f"),
        "Restrictor Inlet ID (mm)": st.column_config.NumberColumn(min_value=0.01, format="%.3f"),
        "Restrictor Outlet ID (mm)": st.column_config.NumberColumn(min_value=0.01, format="%.3f"),
        "Restrictor Length (mm)": st.column_config.NumberColumn(min_value=0.0, format="%.3f"),
        "Count": st.column_config.NumberColumn(min_value=1, step=1),
        "Valve Cv/Kv": st.column_config.NumberColumn(min_value=0.01, format="%.3f"),
        "Bend Radius (mm)": st.column_config.NumberColumn(min_value=0.01, format="%.3f"),
        "Bend Angle (°)": st.column_config.NumberColumn(min_value=1.0, max_value=360.0, format="%.1f"),
    },
)
st.caption("Rows are evaluated top to bottom. All lengths are entered in mm. A bend is a tube component with its own ID, radius, and now angle is only for 90 degree. Barb/orifice rows must sit between two tube rows; their inlet and outlet IDs are checked against their adjacent tubes.")

pressure_drop_pa = (p_inlet_bar - p_outlet_bar) * 100_000
calculate_clicked = st.button("Calculate flow", type="primary")
if calculate_clicked and pressure_drop_pa <= 0:
    st.warning("Outlet pressure must be lower than inlet pressure.")
elif calculate_clicked:
    try:
        components = [
            {
                "name": row["Name"], "type": row["Type"], "id_m": row["Tube ID (mm)"] / 1000,
                "inlet_id_m": row["Restrictor Inlet ID (mm)"] / 1000,
                "outlet_id_m": row["Restrictor Outlet ID (mm)"] / 1000,
                "length_m": row["Tube Length (mm)"] / 1000,
                "restrictor_length_m": row["Restrictor Length (mm)"] / 1000, "count": row["Count"],
                "valve_coefficient": row["Valve Cv/Kv"],
                "bend_radius_m": row["Bend Radius (mm)"] / 1000,
                "bend_angle_deg": row["Bend Angle (°)"],
            }
            for _, row in path.iterrows()
        ]
        result = detailed_flow.solve_detailed_flow(
            pressure_inlet_pa=p_inlet_bar * 100_000,
            pressure_outlet_pa=p_outlet_bar * 100_000,
            temperature_k=temperature_c + 273.15,
            fluid=fluid,
            path_components=components,
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
        breakdown["ΔP (bar)"] = breakdown["ΔP (Pa)"] / 100_000
        st.subheader("Pressure-loss breakdown")
        st.dataframe(breakdown[["Component", "Type", "ΔP (bar)", "Re"]], width="stretch", hide_index=True)
        fig, ax = plt.subplots(figsize=(9, 3.5))
        ax.bar(breakdown["Component"], breakdown["ΔP (bar)"], color="#2b5c8f")
        ax.set_ylabel("Pressure loss (bar)")
        ax.set_title("Loss by path component")
        ax.tick_params(axis="x", rotation=20)
        st.pyplot(fig)
    except (ValueError, ZeroDivisionError, TypeError) as error:
        st.error(f"Unable to calculate this path: {error}")
