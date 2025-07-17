import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.interpolate import interp1d

# --- CUSTOM CSS for WHITE THEME ---
st.markdown("""
<style>
    /* Set white background and black text */
    .main {
        background-color: #ffffff !important;
        color: #111111 !important;
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    }
    /* Sidebar white with subtle border */
    .css-1d391kg {
        background-color: #fafafa !important;
        border-right: 1px solid #ddd;
    }
    /* Buttons styling */
    div.stButton > button {
        background-color: #0078d4;
        color: white;
        font-weight: 600;
        padding: 0.6rem 1.5rem;
        border-radius: 8px;
        border: none;
        transition: background-color 0.3s ease;
    }
    div.stButton > button:hover {
        background-color: #005a9e;
        cursor: pointer;
    }
    /* Headers styling */
    h1, h2, h3 {
        color: #222222;
    }
    /* Dataframe styling */
    .dataframe th, .dataframe td {
        padding: 0.5rem 0.8rem;
        text-align: center;
    }
    /* Plotly chart container */
    .stPlotlyChart > div {
        border-radius: 12px;
        box-shadow: 0 0 15px rgba(0,0,0,0.05);
    }
</style>
""", unsafe_allow_html=True)

# --- PAGE CONFIG ---
st.set_page_config(
    page_title="PAMA Rheology Models",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- APP HEADER ---
st.title("🧪 PAMA Rheology Modeling App")
st.markdown(
    "<p style='font-size:18px; color:#444;'>Explore different polymer rheology models with interactive parameters and clean visualizations.</p>",
    unsafe_allow_html=True
)
st.markdown("---")

# --- HELPER: PLOT FUNCTION ---
def make_plot(df, title, ycols, labels):
    fig = go.Figure()
    for ycol, label in zip(ycols, labels):
        fig.add_trace(go.Scatter(
            x=df["shear"], y=df[ycol], mode='lines+markers',
            name=label, line=dict(width=2), marker=dict(size=4)
        ))
    fig.update_layout(
        title=title,
        xaxis=dict(title="Shear rate (s⁻¹)", type="log", gridcolor='lightgray'),
        yaxis=dict(title="Viscosity (cP)", type="log", gridcolor='lightgray'),
        template="plotly_white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        margin=dict(t=50, b=40, l=40, r=40),
        hovermode="x unified"
    )
    return fig

# --- MODELS IMPLEMENTATION ---
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
    parc1 = (l * 7.3)**2
    eta7 = 1 + (eta_0 - 1) * (1 + parc1)**parc2

    f = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")
    ceta_exp = f(eta7_exp)

    eta_0C = 1 + ceta_exp + 0.582 * ceta_exp**2.009 + 0.022 * ceta_exp**4
    parcC = -0.08212 * ceta_exp
    nC = 1 - (0.6187 - 0.5203 * np.exp(parcC))
    l_dC = 0.251 + 1.54 * MW * ceta_exp / (C * Temp)
    lC = l_dC * (0.810 + 0.0230 * ceta_exp**2.438)

    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 1000, 1)])
    etaC = eta_in + (eta_0C - eta_in) * (1 + (lC * shear)**2)**((nC - 1) / 2)

    return pd.DataFrame({"shear": shear, "Viscosity (cP)": etaC})

def model_pama_temperature(C, MW, eta7_exp, T_wanted_C):
    T = T_wanted_C + 273
    Temp = 298
    eta_in = np.exp(-3.7188 + (578.919 / (-137.546 + Temp)))
    eta_in_T = np.exp(-3.7188 + (578.919 / (-137.546 + T)))

    ceta = np.arange(0.1, 100.1, 0.1)
    eta_0 = 1 + ceta + 0.582 * ceta**2.009 + 0.022 * ceta**4
    parc = -0.08212 * ceta
    n = 1 - (0.6187 - 0.5203 * np.exp(parc))
    parc2 = (n - 1) / 2
    l_d = 0.251 + 1.54 * MW * ceta / (C * Temp)
    l = l_d * (0.810 + 0.0230 * ceta**2.438)
    eta7 = 1 + (eta_0 - 1) * (1 + (l * 7.3)**2)**parc2

    f = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")
    ceta_exp = f(eta7_exp)

    eta_0C = 1 + ceta_exp + 0.582 * ceta_exp**2.009 + 0.022 * ceta_exp**4
    parcC = -0.08212 * ceta_exp
    nC = 1 - (0.6187 - 0.5203 * np.exp(parcC))
    l_dC = 0.251 + 1.54 * MW * ceta_exp / (C * Temp)
    lC = l_dC * (0.810 + 0.0230 * ceta_exp**2.438)
    lC_T = 0.251 + 1.54 * MW * ceta_exp / (C * T)
    lC_T = lC_T * (0.810 + 0.0230 * ceta_exp**2.438)

    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 1000, 1)])
    etaC_ref = eta_in + (eta_0C - eta_in) * (1 + (lC * shear)**2)**((nC - 1) / 2)
    etaC_temp = eta_in_T + (eta_0C - eta_in_T) * (1 + (lC_T * shear)**2)**((nC - 1) / 2)

    return pd.DataFrame({"shear": shear, "25°C Reference": etaC_ref, f"{T_wanted_C}°C Model": etaC_temp})

def model_pama_degradation(C, MW, eta7_exp, eta7_exp_D):
    Temp = 298
    eta_in = np.exp(-3.7188 + (578.919 / (-137.546 + Temp)))

    ceta = np.arange(0.1, 100.1, 0.1)
    eta_0 = 1 + ceta + 0.582 * ceta**2.009 + 0.022 * ceta**4
    parc = -0.08212 * ceta
    n = 1 - (0.6187 - 0.5203 * np.exp(parc))
    parc2 = (n - 1) / 2
    l_d = 0.251 + 1.54 * MW * ceta / (C * Temp)
    l = l_d * (0.810 + 0.0230 * ceta**2.438)
    parc1 = (l * 7.3)**2
    eta7 = 1 + (eta_0 - 1) * (1 + parc1)**parc2

    f = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")
    ceta_exp = f(eta7_exp)
    ceta_exp_D = f(eta7_exp_D)

    eta_0C = 1 + ceta_exp + 0.582 * ceta_exp**2.009 + 0.022 * ceta_exp**4
    parcC = -0.08212 * ceta_exp
    nC = 1 - (0.6187 - 0.5203 * np.exp(parcC))
    l_dC = 0.251 + 1.54 * MW * ceta_exp / (C * Temp)
    lC = l_dC * (0.810 + 0.0230 * ceta_exp**2.438)

    eta_0D = 1 + ceta_exp_D + 0.582 * ceta_exp_D**2.009 + 0.022 * ceta_exp_D**4
    parcD = -0.08212 * ceta_exp_D
    nD = 1 - (0.6187 - 0.5203 * np.exp(parcD))
    l_dD = 0.251 + 1.54 * MW * ceta_exp_D / (C * Temp)
    lD = l_dD * (0.810 + 0.0230 * ceta_exp_D**2.438)

    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 1000, 1)])
    polymer_UD = eta_in + (eta_0C - eta_in) * (1 + (lC * shear)**2)**((nC - 1) / 2)
    polymer_D = eta_in + (eta_0D - eta_in) * (1 + (lD * shear)**2)**((nD - 1) / 2)

    return pd.DataFrame({"shear": shear, "Polymer UD": polymer_UD, "Polymer Degraded": polymer_D})

# --- APP LAYOUT with TABS ---
tab1, tab2, tab3 = st.tabs(["Basic PAMA", "PAMA with Temperature", "PAMA with Degradation"])

with tab1:
    st.header("Basic PAMA Model")
    st.markdown("Input parameters to calculate viscosity using the basic PAMA model.")
    C = st.slider("Concentration (g/L)", 0.1, 20.0, 2.0, 0.1, help="Polymer concentration in grams per liter")
    MW = st.slider("Molecular Weight (MDa)", 0.1, 50.0, 8.0, 0.1, help="Molecular weight in MegaDaltons")
    eta7_exp = st.number_input("η@7.3 experimental (cP)", min_value=0.01, value=20.0, format="%.3f")
    if st.button("Run Basic Model"):
        df_basic = model_basic_pama(C, MW, eta7_exp)
        st.plotly_chart(make_plot(df_basic, "Basic PAMA Viscosity vs Shear Rate", ["Viscosity (cP)"], ["Viscosity"]), use_container_width=True)
        st.dataframe(df_basic.style.format("{:.3f}"))
        st.download_button("Download Data (CSV)", df_basic.to_csv(index=False), "basic_pama.csv", "text/csv")

with tab2:
    st.header("PAMA with Temperature Effects")
    st.markdown("Explore how temperature influences viscosity in the PAMA model.")
    C = st.slider("Concentration (g/L)", 0.1, 20.0, 2.0, 0.1)
    MW = st.slider("Molecular Weight (MDa)", 0.1, 50.0, 8.0, 0.1)
    eta7_exp = st.number_input("η@7.3 experimental at 25°C (cP)", min_value=0.01, value=15.653, format="%.3f")
    T_wanted_C = st.slider("Target Temperature (°C)", 0, 100, 35, 1)
    if st.button("Run Temperature Model"):
        df_temp = model_pama_temperature(C, MW, eta7_exp, T_wanted_C)
        st.plotly_chart(make_plot(df_temp, "PAMA Temperature Model Viscosity", ["25°C Reference", f"{T_wanted_C}°C Model"], ["25°C", f"{T_wanted_C}°C"]), use_container_width=True)
        st.dataframe(df_temp.style.format("{:.3f}"))
        st.download_button("Download Data (CSV)", df_temp.to_csv(index=False), "pama_temperature.csv", "text/csv")

with tab3:
    st.header("PAMA with Degradation Effects")
    st.markdown("Compare viscosity before and after polymer degradation.")
    C = st.slider("Concentration (g/L)", 0.1, 20.0, 2.0, 0.1)
    MW = st.slider("Molecular Weight (MDa)", 0.1, 50.0, 8.0, 0.1)
    eta7_exp = st.number_input("η@7.3 experimental (Polymer UD) (cP)", min_value=0.01, value=15.653, format="%.3f")
    eta7_exp_D = st.number_input("η@7.3 experimental (Polymer D) (cP)", min_value=0.01, value=7.354, format="%.3f")
    if st.button("Run Degradation Model"):
        df_deg = model_pama_degradation(C, MW, eta7_exp, eta7_exp_D)
        st.plotly_chart(make_plot(df_deg, "PAMA Degradation Model Viscosity", ["Polymer UD", "Polymer Degraded"], ["Polymer UD", "Degraded"]), use_container_width=True)
        st.dataframe(df_deg.style.format("{:.3f}"))
        st.download_button("Download Data (CSV)", df_deg.to_csv(index=False), "pama_degradation.csv", "text/csv")

# --- FOOTER ---
st.markdown("---")
st.markdown("<p style='text-align:center; color:#666;'>Made by Eduar — Powered by Streamlit & Plotly</p>", unsafe_allow_html=True)
