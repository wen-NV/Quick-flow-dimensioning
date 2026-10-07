import sys
import importlib
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from thermo import Chemical, Mixture

import detailed_flow

# Streamlit retains imported modules between reruns; reload this local solver
# so edits to detailed_flow.py are applied without restarting the app.
detailed_flow = importlib.reload(detailed_flow)


# ✅ Get the current directory dynamically
#BASE_DIR = Path(__file__).resolve().parent
#PROJECT_PATH = BASE_DIR / "Process_Parameters"



import importlib
import urllib.request

@st.cache_data(ttl=300)
def sync_fluid_functions():
    url = "https://raw.githubusercontent.com/wen-NV/general_fluid_function/main/general_function.py"
    urllib.request.urlretrieve(url, "general_function.py")

# 1. Download file
sync_fluid_functions()

# 2. Tell Python to refresh its file scanner
importlib.invalidate_caches()

# 3. Import safely
import general_function as process_functions


st.set_page_config(page_title="Detailed Fluid Path Calculator", layout="wide")
st.title("🌊 Detailed Fluid Path Calculator")
st.write("Solve flow rate from tube friction and a selectable restriction on each tube section.")

with st.expander("How to use this app", expanded=False):
    st.markdown(
        """
1. Set inlet/outlet pressure, temperature, and either a pure liquid or a two-liquid mixture with mole fractions (can be changed to mass fraction or volume fraction).
2. Build the flow path in the table from top to bottom. One row is one component.
3. Select the component type first, then fill only the fields that are not marked `—`.
4. The Name field is optional. You can also name each component for later tracking, or it will be named automatically with type and number.
5. Use **Move** `↑` or `↓` to change the component order. 
6. Use the **+** button to add a new row at the bottom, or click the whole role and click "trash bin" at the top right of the table to delete a row.
7. Click **Calculate flow** to see total flow, fluid properties, and pressure loss by component.

    **Path rule:** Fittings and valves may connect directly in series. For a fitting-to-fitting interface, the calculation uses the previous component's outlet ID and the next component's inlet ID. A final fitting may discharge directly to the outlet.

All path-table dimensions are in **mm**. This model is for liquids and gas in low pressure differences, not compressible gas flow for now.
        """
    )

st.header("1. Fluid and pressure")
c1, c2, c3, c4 = st.columns(4)
with c1:
    p_inlet_bar = st.number_input("Inlet absolute pressure (bar(a))", min_value=0.001, value=1.8, step=0.1)
with c2:
    p_outlet_bar = st.number_input("Outlet absolute pressure (bar(a))", min_value=0.001, value=1.0, step=0.1)
with c3:
    fluid_mode = st.selectbox("Fluid type", ["Pure fluid", "Two-component mixture"])
with c4:
    temperature_c = st.number_input("Temperature (°C)", value=25.0, step=1.0)

fluid_options = ["water", "acetone", "ethanol", "methanol", "toluene", "MTBE", "air", "nitrogen", "hydrogen"]
if fluid_mode == "Pure fluid":
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
st.caption("Density and viscosity come from thermo at the selected temperature and mean pressure of inlet and outlet.")

gas_fluids = {"air", "nitrogen", "hydrogen"}
if fluid_mode == "Pure fluid" and fluid.lower() in gas_fluids:
    gas = Chemical(
        fluid,
        T=temperature_c + 273.15,
        P=((p_inlet_bar + p_outlet_bar) / 2) * 100_000,
    )
    gamma = gas.isentropic_exponent or 1.4
    critical_ratio = (2 / (gamma + 1)) ** (gamma / (gamma - 1))
    pressure_ratio = p_outlet_bar / p_inlet_bar
    critical_pressure_bar = p_inlet_bar * critical_ratio
    st.info(
        f"Gas screen — γ = {gamma:.3f}; Pout/Pin = {pressure_ratio:.3f}; "
        f"ideal-gas critical ratio = {critical_ratio:.3f}; "
        f"critical outlet pressure = {critical_pressure_bar:.3f} bar(a)."
    )
    if pressure_ratio <= critical_ratio:
        st.warning(
            "Possible choked flow: the overall outlet pressure is at or below the ideal-gas critical pressure. "
            "A compressible gas model is required; the current liquid-model flow result is not valid."
        )
    else:
        st.warning(
            "Gas selected: the global pressure ratio is above the ideal-gas choking threshold, but a local "
            "restriction may still choke. Use a compressible gas model for a design-quality result."
        )

st.header("2. Ordered flow path")
component_types = [
    "Straight tube", "Bend tube", "Barb", "LBarb", "Orifice",
    "Compression fitting", "Valve Cv", "Valve Kv",
]
default_path = pd.DataFrame([
    {"Name": "Feed tube", "Type": "Straight tube", "Tube ID (mm)": "4.0", "Tube length (mm)": "25.0", "Restrictor inlet ID (mm)": "—", "Restrictor outlet ID (mm)": "—", "Restrictor length (mm)": "—", "Valve Cv/Kv": "—", "Bend radius (mm)": "—", "Move": ""},
    {"Name": "Barb", "Type": "Barb", "Tube ID (mm)": "—", "Tube length (mm)": "—", "Restrictor inlet ID (mm)": "2.4", "Restrictor outlet ID (mm)": "2.4", "Restrictor length (mm)": "34.0", "Valve Cv/Kv": "—", "Bend radius (mm)": "—", "Move": ""},
    {"Name": "Bend", "Type": "Bend tube", "Tube ID (mm)": "3.175", "Tube length (mm)": "40.0", "Restrictor inlet ID (mm)": "—", "Restrictor outlet ID (mm)": "—", "Restrictor length (mm)": "—", "Valve Cv/Kv": "—", "Bend radius (mm)": "50.0", "Move": ""},
    {"Name": "Outlet tube", "Type": "Straight tube", "Tube ID (mm)": "4.0", "Tube length (mm)": "25.0", "Restrictor inlet ID (mm)": "—", "Restrictor outlet ID (mm)": "—", "Restrictor length (mm)": "—", "Valve Cv/Kv": "—", "Bend radius (mm)": "—", "Move": ""},
])
if "path_table" not in st.session_state:
    st.session_state.path_table = default_path
elif "Move" not in st.session_state.path_table.columns:
    st.session_state.path_table = st.session_state.path_table.copy()
    st.session_state.path_table["Move"] = ""
if "Order" not in st.session_state.path_table.columns:
    st.session_state.path_table = st.session_state.path_table.copy()
    st.session_state.path_table.insert(0, "Order", range(1, len(st.session_state.path_table) + 1))
else:
    st.session_state.path_table["Order"] = range(1, len(st.session_state.path_table) + 1)
if "path_editor_revision" not in st.session_state:
    st.session_state.path_editor_revision = 0

type_fields = {
    "Straight tube": {"Tube ID (mm)", "Tube length (mm)"},
    "Bend tube": {"Tube ID (mm)", "Tube length (mm)", "Bend radius (mm)"},
    "Barb": {"Restrictor inlet ID (mm)", "Restrictor outlet ID (mm)", "Restrictor length (mm)"},
    "Orifice": {"Restrictor inlet ID (mm)", "Restrictor outlet ID (mm)", "Restrictor length (mm)"},
    "LBarb": {"Restrictor inlet ID (mm)", "Restrictor outlet ID (mm)", "Restrictor length (mm)"},
    "Compression fitting": {"Restrictor inlet ID (mm)", "Restrictor outlet ID (mm)", "Restrictor length (mm)"},
    "Valve Cv": {"Restrictor inlet ID (mm)", "Restrictor outlet ID (mm)", "Valve Cv/Kv"},
    "Valve Kv": {"Restrictor inlet ID (mm)", "Restrictor outlet ID (mm)", "Valve Cv/Kv"},
}
input_columns = set().union(*type_fields.values())

def normalize_path(table):
    """Set non-applicable cells to — and clear newly applicable cells."""
    normalized = table.copy()
    for index, row in normalized.iterrows():
        kind = row["Type"] if row["Type"] in type_fields else "Straight tube"
        normalized.at[index, "Type"] = kind
        for column in input_columns:
            if column not in type_fields[kind]:
                normalized.at[index, column] = "—"
            else:
                value = str(normalized.at[index, column]).strip()
                is_valve_bore = kind in {"Valve Cv", "Valve Kv"} and column in {
                    "Restrictor inlet ID (mm)", "Restrictor outlet ID (mm)"
                }
                if is_valve_bore and value in {"", "—", "-", "None", "nan", "optional"}:
                    normalized.at[index, column] = "optional"
                elif value in {"—", "-", "None", "nan", "optional"}:
                    normalized.at[index, column] = ""
    return normalized

def sync_path_table(editor_key):
    """Save table edits and set non-applicable component fields to —."""
    changes = st.session_state[editor_key]
    updated = st.session_state.path_table.copy()
    for row_index, values in changes.get("edited_rows", {}).items():
        for column, value in values.items():
            updated.at[row_index, column] = value
    for row_index, values in changes.get("edited_rows", {}).items():
        direction = values.get("Move")
        target_index = row_index - 1 if direction == "↑" else row_index + 1 if direction == "↓" else row_index
        if direction in {"↑", "↓"} and 0 <= target_index < len(updated):
            order = list(range(len(updated)))
            order[row_index], order[target_index] = order[target_index], order[row_index]
            updated = updated.iloc[order].reset_index(drop=True)
    for row_index in sorted(changes.get("deleted_rows", []), reverse=True):
        updated = updated.drop(updated.index[row_index])
    added_rows = changes.get("added_rows", [])
    if added_rows:
        updated = pd.concat([updated, pd.DataFrame(added_rows)], ignore_index=True)
    if "Move" not in updated.columns:
        updated["Move"] = ""
    updated["Move"] = ""
    if "Order" not in updated.columns:
        updated.insert(0, "Order", range(1, len(updated) + 1))
    else:
        updated["Order"] = range(1, len(updated) + 1)
    st.session_state.path_table = normalize_path(updated).reset_index(drop=True)
    st.session_state.path_editor_revision += 1

editor_key = f"path_editor_{st.session_state.path_editor_revision}"
path = st.data_editor(
    st.session_state.path_table,
    num_rows="dynamic",
    width="stretch",
    key=editor_key,
    on_change=sync_path_table,
    args=(editor_key,),
    disabled=["Order"],
    column_config={
        "Order": st.column_config.NumberColumn("Order", format="%d", width="small"),
        "Name": st.column_config.TextColumn(),
        "Type": st.column_config.SelectboxColumn(options=component_types, required=True),
        "Tube ID (mm)": st.column_config.TextColumn(),
        "Tube length (mm)": st.column_config.TextColumn(),
        "Restrictor inlet ID (mm)": st.column_config.TextColumn(),
        "Restrictor outlet ID (mm)": st.column_config.TextColumn(),
        "Restrictor length (mm)": st.column_config.TextColumn(),
        "Valve Cv/Kv": st.column_config.TextColumn(),
        "Bend radius (mm)": st.column_config.TextColumn(),
        "Move": st.column_config.SelectboxColumn(options=["", "↑", "↓"], width="small"),
    },
)
st.caption("Use Move ↑ or ↓ to reorder a row. Fields not used by the selected component type reset to —.")


validation_errors = []

def required_number(row, column):
    value = str(row[column]).strip()
    if value in {"", "-", "—"}:
        validation_errors.append(f"{row['Name']}: enter a value for {column}.")
        return 0.0
    try:
        return float(value)
    except ValueError:
        validation_errors.append(f"{row['Name']}: {column} must be a number.")
        return 0.0

def optional_number(row, column):
    """Return an optional valve-bore value in mm, or None when left optional."""
    value = str(row[column]).strip()
    if value.lower() in {"", "-", "—", "optional", "none", "nan"}:
        return None
    try:
        number = float(value)
    except ValueError:
        validation_errors.append(f"{row['Name']}: {column} must be a number or optional.")
        return None
    if number <= 0:
        validation_errors.append(f"{row['Name']}: {column} must be greater than zero when provided.")
        return None
    return number

components = []
for _, row in path.iterrows():
    kind = row["Type"]
    component = {"name": row["Name"], "type": kind, "count": 1, "bend_angle_deg": 90.0}
    if kind in {"Straight tube", "Bend tube"}:
        component["id_m"] = required_number(row, "Tube ID (mm)") / 1000
        component["length_m"] = required_number(row, "Tube length (mm)") / 1000
        if kind == "Bend tube":
            component["bend_radius_m"] = required_number(row, "Bend radius (mm)") / 1000
    elif kind in {"Barb", "Orifice", "LBarb", "Compression fitting"}:
        component["inlet_id_m"] = required_number(row, "Restrictor inlet ID (mm)") / 1000
        component["outlet_id_m"] = required_number(row, "Restrictor outlet ID (mm)") / 1000
        component["restrictor_length_m"] = required_number(row, "Restrictor length (mm)") / 1000
    else:
        component["valve_coefficient"] = required_number(row, "Valve Cv/Kv")
        inlet_id = optional_number(row, "Restrictor inlet ID (mm)")
        outlet_id = optional_number(row, "Restrictor outlet ID (mm)")
        if inlet_id is not None:
            component["inlet_id_m"] = inlet_id / 1000
        if outlet_id is not None:
            component["outlet_id_m"] = outlet_id / 1000
    components.append(component)

pressure_drop_pa = (p_inlet_bar - p_outlet_bar) * 100_000
calculate_clicked = st.button("Calculate flow", type="primary")
if calculate_clicked and pressure_drop_pa <= 0:
    st.warning("Outlet pressure must be lower than inlet pressure.")
elif calculate_clicked and validation_errors:
    st.error(" ".join(validation_errors))
elif calculate_clicked:
    try:
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
        # Format as display text so all result-table values align consistently.
        breakdown["ΔP (bar)"] = (breakdown["ΔP (Pa)"] / 100_000).map(lambda value: f"{value:.4f}")
        breakdown["Reynolds number"] = breakdown["Re"].map(lambda value: f"{value:,.0f}")
        st.subheader("Pressure-loss breakdown")
        st.dataframe(
            breakdown[["Component", "Type", "ΔP (bar)", "Reynolds number"]],
            width="stretch",
            hide_index=True,
        )
        fig, ax = plt.subplots(figsize=(9, 3.5))
        ax.bar(breakdown["Component"], breakdown["ΔP (bar)"], color="#2b5c8f")
        ax.set_ylabel("Pressure loss (bar)")
        ax.set_title("Loss by path component")
        ax.tick_params(axis="x", rotation=20)
        st.pyplot(fig)
    except (ValueError, ZeroDivisionError, TypeError) as error:
        st.error(f"Unable to calculate this path: {error}")
