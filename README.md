# Detailed Fluid Path Calculator

A Streamlit app for estimating liquid flow through an ordered fluid path. It balances the available pressure drop against tube friction and local losses from bends, barbs, L-barbs, orifices, and valves.

## Start the app

Create and activate a Python environment, then install the dependencies:

```bash
pip install -r requirements.txt
streamlit run quick_flow_dimensioning.py
```

The app downloads the latest `general_function.py` helper from the configured public GitHub repository. Internet access is required when that helper is refreshed.

## 1. Set fluid and pressure conditions

Enter:

- **Inlet pressure** and **outlet pressure** in bar.
- **Temperature** in °C.
- A pure liquid, or a two-liquid mixture and the mole fraction of component 1.

The app uses `thermo` to obtain density and dynamic viscosity at the specified temperature and mean system pressure. The calculated density and viscosity are displayed in the results.

This model is intended for incompressible liquids. Do not use the liquid-flow result for nitrogen, hydrogen, or other gases without a compressible-flow model.

## 2. Build the ordered flow path

The table is read from top to bottom. One row represents one physical component. The **Name** column is optional and is only used as a result label.

Use the built-in table controls to add or delete rows. Use the **Move** column to rearrange a row:

- `↑` moves it up one position.
- `↓` moves it down one position.

Fields that do not apply to the selected component type are automatically shown as `—` and are ignored by the calculation. All lengths and diameters in the table are in **mm**.

| Component type | Required inputs |
| --- | --- |
| Straight tube | Tube ID, Tube length |
| Bend tube | Tube ID, Tube length, Bend radius |
| Barb | Restrictor inlet ID, Restrictor outlet ID, Restrictor length |
| LBarb | Restrictor inlet ID, Restrictor outlet ID, Restrictor length |
| Orifice | Restrictor inlet ID, Restrictor outlet ID, Restrictor length |
| Compression fitting | Restrictor inlet ID, Restrictor outlet ID, Restrictor length |
| Valve Cv | Valve Cv/Kv (enter Cv) |
| Valve Kv | Valve Cv/Kv (enter Kv) |

### Path rules

- A barb, L-barb, or orifice normally sits between two tube rows.
- A compression fitting uses the same inlet/outlet-ID and length inputs as a barb, and normally sits between two tube rows.
- A Cv or Kv valve must also follow a tube and normally sit between two tube rows.
- Its inlet ID must fit the preceding tube ID; its outlet ID must fit the following tube ID.
- A barb, L-barb, or orifice may be the final row. It is then treated as a free outlet discharging into a large volume. A valve may also be the final row.
- A restrictor that is not the final row must be followed by a straight or bend tube. The app warns about invalid ordering and blocks calculation until it is corrected.

Example path:

```text
Straight tube → Barb → Bend tube → LBarb → Straight tube → Orifice
```

The final orifice is valid because it is the outlet.

## 3. Calculate and interpret results

Click **Calculate flow** after configuring the path.

The app displays:

- Solved volumetric flow rate in L/min
- Available pressure drop and calculated losses in bar
- Fluid density and dynamic viscosity
- Per-component pressure-loss and Reynolds-number table
- Pressure-loss bar chart

The solver uses Darcy–Weisbach tube losses and the hydraulic helper functions in `general_function.py` for friction factors and local-loss coefficients.

## Notes and limitations

- Confirm the geometry, units, material properties, and component coefficients independently before using results for design decisions.
- The bend-loss data is based on the available 90° bend correlations; bend tubes use a fixed 90° angle in the current UI.
- Valve Cv/Kv calculations use liquid-flow conventions and fluid specific gravity relative to water.
