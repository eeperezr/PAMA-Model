import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

# --- Model functions (unchanged, from your original code) ---
def model_basic_pama(C, MW, eta7_exp):
    Temp = 298
    eta_in = np.exp(-3.7188 + (578.919 / (-137.546 + Temp)))
    ceta = np.arange(0.1, 100.1, 0.1)
    eta_0 = 1 + ceta + 0.582 * ceta**2.009 + 0.022 * ceta**4
    parc = -0.08212 * ceta
    n = 1 - (0.6187 - 0.5203 * np.exp(parc))
    parc2 = (n - 1) / 2
    l_d = 0.251 + 1.54 * MW * ceta / (C * Temp)
    l = l_d * (0.810 + 0.0230 * ceta**2.438)
    parc1 = (l * 7.3) ** 2
    eta7 = 1 + (eta_0 - 1) * (1 + parc1) ** parc2
    f = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")
    ceta_exp = f(eta7_exp)
    eta_0C = 1 + ceta_exp + 0.582 * ceta_exp**2.009 + 0.022 * ceta_exp**4
    parcC = -0.08212 * ceta_exp
    nC = 1 - (0.6187 - 0.5203 * np.exp(parcC))
    l_dC = 0.251 + 1.54 * MW * ceta_exp / (C * Temp)
    lC = l_dC * (0.810 + 0.0230 * ceta_exp**2.438)
    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 1000, 1)])
    etaC = eta_in + (eta_0C - eta_in) * (1 + (lC * shear) ** 2) ** ((nC - 1) / 2)
    return {"shear": shear.tolist(), "Viscosity (cP)": etaC.tolist()}

def model_pama_temperature(C, MW, eta7_exp, T_wanted_C):
    Temp = 298
    T_wanted = T_wanted_C + 273
    eta_in = np.exp(-3.7188 + (578.919 / (-137.546 + Temp)))
    eta_in_T = np.exp(-3.7188 + (578.919 / (-137.546 + T_wanted)))
    ceta = np.arange(0.1, 100.1, 0.1)
    eta_0 = 1 + ceta + 0.582 * ceta**2.009 + 0.022 * ceta**4
    parc = -0.08212 * ceta
    n = 1 - (0.6187 - 0.5203 * np.exp(parc))
    parc2 = (n - 1) / 2
    l_d = 0.251 + 1.54 * MW * ceta / (C * Temp)
    l = l_d * (0.810 + 0.0230 * ceta**2.438)
    parc1 = (l * 7.3) ** 2
    eta7 = 1 + (eta_0 - 1) * (1 + parc1) ** parc2
    f = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")
    ceta_exp = f(eta7_exp)
    eta_0C = 1 + ceta_exp + 0.582 * ceta_exp**2.009 + 0.022 * ceta_exp**4
    parcC = -0.08212 * ceta_exp
    nC = 1 - (0.6187 - 0.5203 * np.exp(parcC))
    l_dC = 0.251 + 1.54 * MW * ceta_exp / (C * Temp)
    lC = l_dC * (0.810 + 0.0230 * ceta_exp**2.438)
    lC_T = 0.251 + 1.54 * MW * ceta_exp / (C * T_wanted)
    lC_T = lC_T * (0.810 + 0.0230 * ceta_exp**2.438)
    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 1000, 1)])
    etaC_ref = eta_in + (eta_0C - eta_in) * (1 + (lC * shear)**2)**((nC - 1) / 2)
    etaC_temp = eta_in_T + (eta_0C - eta_in_T) * (1 + (lC_T * shear)**2)**((nC - 1) / 2)
    return {"shear": shear.tolist(), "25°C Reference": etaC_ref.tolist(), f"{T_wanted_C}°C Model": etaC_temp.tolist()}

def model_pama_degradation(C, MW, eta7_exp, eta7_exp_D):
    alpha = 0.763
    Temp = 298
    eta_in = np.exp(-3.7188 + (578.919 / (-137.546 + Temp)))
    c_med = C * 0.5
    ceta = np.arange(0.1, 100.1, 0.1)
    eta_0 = 1 + ceta + 0.582 * ceta**2.009 + 0.022 * ceta**4
    etab = ceta / (2 ** alpha)
    cetab = etab * c_med
    eta_0_FD = 1 + cetab + 0.582 * cetab**2.009 + 0.022 * cetab**4
    parc = -0.08212 * ceta
    parcb = -0.08212 * cetab
    n = 1 - (0.6187 - 0.5203 * np.exp(parc))
    n_FD = 1 - (0.6187 - 0.5203 * np.exp(parcb))
    parc2 = (n - 1) / 2
    parc2b = (n_FD - 1) / 2
    l_d = 0.251 + 1.54 * MW * ceta / (C * Temp)
    l_d_b = 0.251 + 1.54 * MW * cetab / (C * Temp)
    l = l_d * (0.810 + 0.0230 * ceta**2.438)
    lb = l_d_b * (0.810 + 0.0230 * cetab**2.438)
    parc1 = (l * 7.3) ** 2
    parc1b = (lb * 7.3) ** 2
    eta7 = 1 + (eta_0 - 1) * (1 + parc1) ** parc2
    eta7_FD = 1 + (eta_0_FD - 1) * (1 + parc1b) ** parc2b
    f = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")
    ceta_exp = f(eta7_exp)
    f_D = interp1d(eta7_FD, cetab, bounds_error=False, fill_value="extrapolate")
    cetab_exp = f_D(eta7_exp_D)
    eta_0C = 1 + ceta_exp + 0.582 * ceta_exp**2.009 + 0.022 * ceta_exp**4
    eta_0C_D = 1 + cetab_exp + 0.582 * cetab_exp**2.009 + 0.022 * cetab_exp**4
    parcC = -0.08212 * ceta_exp
    parcC_D = -0.08212 * cetab_exp
    nC = 1 - (0.6187 - 0.5203 * np.exp(parcC))
    nC_D = 1 - (0.6187 - 0.5203 * np.exp(parcC_D))
    l_dC = 0.251 + 1.54 * MW * ceta_exp / (C * Temp)
    l_dC_D = 0.251 + 1.54 * MW * cetab_exp / (C * Temp)
    lC = l_dC * (0.810 + 0.0230 * ceta_exp**2.438)
    lC_D = l_dC_D * (0.810 + 0.0230 * cetab_exp**2.438)
    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 1000, 1)])
    etaC = eta_in + (eta_0C - eta_in) * (1 + (lC * shear)**2)**((nC - 1) / 2)
    etaC_D = eta_in + (eta_0C_D - eta_in) * (1 + (lC_D * shear)**2)**((nC_D - 1) / 2)
    return {"shear": shear.tolist(), "Polymer UD": etaC.tolist(), "Polymer Degraded": etaC_D.tolist()}

# --- Start of UI setup ---

# Page config
st.set_page_config(page_title="PAMA Rheology Models", page_icon="🧪", layout="wide")

# Inject CSS for background, fonts, sidebar, and header styling
st.markdown("""
<style>
    /* Background color */
    .reportview-container, .main {
        background: #f5f7fa;
    }

    /* Sidebar background */
    .sidebar .sidebar-content {
        background: #dff9fb;
        color: #0a3d62;
    }

    /* Title styling */
    .title {
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        font-size: 48px;
        font-weight: 700;
        color: #0a3d62;
        text-align: center;
        margin-bottom: 0;
        margin-top: 15px;
    }
    .subtitle {
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        font-size: 20px;
        color: #3c6382;
        text-align: center;
        margin-top: 5px;
        margin-bottom: 40px;
    }

    /* Style buttons */
    div.stButton > button {
        background-color: #0a3d62;
        color: white;
        font-weight: 600;
        border-radius: 8px;
        padding: 8px 20px;
        margin-top: 10px;
        margin-bottom: 20px;
        border: none;
    }
    div.stButton > button:hover {
        background-color: #07406a;
        cursor: pointer;
    }

    /* Center the download button */
    .download-button {
        text-align: center;
        margin-top: 15px;
    }
</style>
""", unsafe_allow_html=True)

# Title and subtitle
st.markdown("""
    <h1 class="title">PAMA Rheology Models</h1>
    <p class="subtitle">Explore polymer viscosity behavior with interactive multi-model comparisons</p>
""", unsafe_allow_html=True)

# Main layout with two columns: left for screenshots, right for controls + plots
col_left, col_right = st.columns([1, 2])  # Adjust widths as needed

with col_left:
    st.markdown("### Example Graphs")
    st.image("https://images.unsplash.com/photo-1503023345310-bd7c1de61c7d?auto=format&fit=crop&w=400&q=60", caption="Model 1 - Viscosity Curve", use_column_width=True)
    st.image("https://images.unsplash.com/photo-1498050108023-c5249f4df085?auto=format&fit=crop&w=400&q=60", caption="Model 2 - Temperature Effect", use_column_width=True)
    st.image("https://images.unsplash.com/photo-1530023367847-0e4f16c69c14?auto=format&fit=crop&w=400&q=60", caption="Model 3 - Polymer Degradation", use_column_width=True)

with col_right:
    st.markdown("### Select Models to Compare")

    # Multi-select model choice
    models_selected = st.multiselect(
        "Choose one or more models",
        options=["Basic PAMA", "PAMA with Temperature", "PAMA with Degradation"],
        default=["Basic PAMA"]
    )

    # Shared parameter inputs
    conc = st.number_input("Concentration (g/L)", min_value=0.1, max_value=20.0, value=2.0, step=0.1, format="%.3f")
    mw = st.number_input("Molecular Weight (MDa)", min_value=0.1, max_value=50.0, value=8.0, step=0.1, format="%.3f")
    eta7 = st.number_input("η@7.3 (cP)", min_value=0.01, value=20.0, format="%.3f")

    # Model-specific inputs:
    temp = 35
    eta7d = 7.354
    if "PAMA with Temperature" in models_selected:
        temp = st.slider("Temperature (°C)", min_value=0, max_value=100, value=35, step=1)
    if "PAMA with Degradation" in models_selected:
        eta7d = st.number_input("η@7.3 (Polymer D) (cP)", min_value=0.01, value=7.354, format="%.3f")

    # Run button
    run = st.button("Run Models")

    # Run models and plot results
    if run and models_selected:
        fig = go.Figure()
        combined_df = pd.DataFrame()

        for model_name in models_selected:
            if model_name == "Basic PAMA":
                data = model_basic_pama(conc, mw, eta7)
                # Plot Basic PAMA curve
                fig.add_trace(go.Scatter(x=data["shear"], y=data["Viscosity (cP)"], mode="lines", name="Basic PAMA", line=dict(color="#1e90ff")))
                df_tmp = pd.DataFrame({"shear": data["shear"], "Basic PAMA": data["Viscosity (cP)"]})
            elif model_name == "PAMA with Temperature":
                data = model_pama_temperature(conc, mw, eta7, temp)
                # Plot both reference and temp curves
                fig.add_trace(go.Scatter(x=data["shear"], y=data["25°C Reference"], mode="lines", name="25°C Reference", line=dict(color="#2e86de")))
                # Find the other temp key dynamically
                other_key = [k for k in data.keys() if k not in ("shear", "25°C Reference")][0]
                fig.add_trace(go.Scatter(x=data["shear"], y=data[other_key], mode="lines", name=f"{temp}°C Model", line=dict(color="#1e3799")))
                df_tmp = pd.DataFrame({"shear": data["shear"], "25°C Reference": data["25°C Reference"], f"{temp}°C Model": data[other_key]})
            elif model_name == "PAMA with Degradation":
                data = model_pama_degradation(conc, mw, eta7, eta7d)
                fig.add_trace(go.Scatter(x=data["shear"], y=data["Polymer UD"], mode="lines", name="Polymer UD", line=dict(color="#1e90ff")))
                fig.add_trace(go.Scatter(x=data["shear"], y=data["Polymer Degraded"], mode="lines", name="Polymer Degraded", line=dict(color="#ff4500")))
                df_tmp = pd.DataFrame({"shear": data["shear"], "Polymer UD": data["Polymer UD"], "Polymer Degraded": data["Polymer Degraded"]})
            else:
                continue

            # Merge dfs on shear for combined table (outer join on shear)
            if combined_df.empty:
                combined_df = df_tmp
            else:
                combined_df = pd.merge(combined_df, df_tmp, on="shear", how="outer")

        # Final plot layout
        fig.update_layout(
            title="PAMA Models Comparison",
            xaxis_title="Shear Rate (s⁻¹)",
            yaxis_title="Viscosity (cP)",
            xaxis_type="log",
            yaxis_type="log",
            template="plotly_white",
            legend=dict(title="Models", x=0.8, y=0.95)
        )

        st.plotly_chart(fig, use_container_width=True)

        # Show combined data table with centered text and borders
        st.markdown("### Combined Data Table")
        st.dataframe(combined_df.style.set_properties(**{'text-align': 'center', 'border': '1px solid #e0e0e0'}))

        # CSV download button for combined data
        csv = combined_df.to_csv(index=False)
        st.download_button(
            label="Download Combined Data as CSV",
            data=csv,
            file_name="pama_models_comparison.csv",
            mime="text/csv"
        )
    elif not run:
        st.info("Select one or more models and click 'Run Models' to see the results here.")

# Sidebar content
with st.sidebar:
    st.title("PAMA Controls")
    st.markdown("""
    ### Navigation
    - [Documentation](#)
      - Overview
      - Setup
      - Manual
      - Developers
      - Contributors
      - Changelog
    - [Source Code](https://example.com)
    - Created by: Eduar Perez (University of Buenos Aires)
    """)

