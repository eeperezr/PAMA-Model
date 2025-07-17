import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import uuid

# --- PAGE CONFIG ---
st.set_page_config(
    page_title="PAMA Rheology Models",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- CUSTOM CSS for WHITE THEME ---
st.markdown(
    """
    <style>
    body {
        background-color: #ffffff;
        color: #222222;
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    }
    /* Sidebar styling */
    .css-1d391kg {
        background-color: #fafafa !important;
        border-right: 1px solid #ddd;
    }
    /* Buttons styling */
    .stButton>button {
        background-color: #0a84ff;
        color: white;
        font-weight: 600;
        padding: 0.6rem 1.5rem;
        border-radius: 8px;
        border: none;
        transition: background-color 0.3s ease;
    }
    .stButton>button:hover {
        background-color: #006fd6;
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
    .block-container {
        padding: 2rem 3rem 3rem 3rem;
    }
    /* Plotly chart container */
    .stPlotlyChart > div {
        border-radius: 12px;
        box-shadow: 0 0 15px rgba(0,0,0,0.05);
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# --- APP HEADER ---
st.title("🧪 PAMA Rheology Modeling App")
st.markdown(
    "<p style='font-size:18px; color:#444;'>Explore polymer rheology models with interactive parameters and visualizations.</p>",
    unsafe_allow_html=True,
)
st.markdown("---")

# --- HELPER: PLOT FUNCTION ---
def make_plot(df, title, ycols, labels):
    fig = go.Figure()
    for ycol, label in zip(ycols, labels):
        fig.add_trace(
            go.Scatter(
                x=df["shear"],
                y=df[ycol],
                mode="lines+markers",
                name=label,
                line=dict(width=2),
                marker=dict(size=4),
            )
        )
    fig.update_layout(
        title=title,
        xaxis=dict(title="Shear rate (s⁻¹)", type="log", gridcolor="lightgray"),
        yaxis=dict(title="Viscosity (cP)", type="log", gridcolor="lightgray"),
        template="plotly_white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5),
        margin=dict(t=50, b=40, l=40, r=40),
        hovermode="x unified",
    )
    return fig

# --- MODEL IMPLEMENTATIONS ---
# Model 1: Basic PAMA
def model_basic_pama(C, MW, eta7_exp):
    Temp = 298  # 25°C in Kelvin
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

    return pd.DataFrame({"shear": shear, "Viscosity (cP)": etaC})

# Model 2: PAMA with Temperature
def model_pama_temperature(C, MW, eta7_exp, T_wanted_C):
    Temp = 298  # Reference temp 25°C in Kelvin
    T_wanted = T_wanted_C + 273  # Convert to Kelvin
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

    return pd.DataFrame({"shear": shear, "25°C Reference": etaC_ref, f"{T_wanted_C}°C Model": etaC_temp})

# Model 3: PAMA with Degradation
def model_pama_degradation(C, MW, eta7_exp, eta7_exp_D):
    alpha = 0.763
    Temp = 298  # 25°C in Kelvin
    eta_in = np.exp(-3.7188 + (578.919 / (-137.546 + Temp)))

    c_med = C * 0.5
    ceta = np.arange(0.1, 100.1, 0.1)
    eta_ Ang = eta / (2 ** alpha)
    cetab = eta * c_med
    eta_0 = 1 + ceta + 0.582 * ceta**2.009 + 0.022 * ceta**4
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

    return pd.DataFrame({"shear": shear, "Polymer UD": etaC, "Polymer Degraded": etaC_D})

# --- SIDEBAR for MODEL SELECTION and PARAMETERS ---
st.sidebar.header("Model Selection and Parameters")
model_choice = st.sidebar.radio("Select Model", ("Basic PAMA", "PAMA with Temperature", "PAMA with Degradation"))

# --- MAIN CONTENT ---
if model_choice == "Basic PAMA":
    st.header("Basic PAMA Model")
    st.markdown("Input parameters to calculate viscosity using the basic PAMA model.")
    C = st.sidebar.slider("Concentration (g/L)", 0.1, 20.0, 2.0, 0.1, key="C_basic")
    MW = st.sidebar.slider("Molecular Weight (MDa)", 0.1, 50.0, 8.0, 0.1, key="MW_basic")
    eta7_exp = st.sidebar.number_input("η@7.3 experimental (cP)", min_value=0.01, value=20.0, format="%.3f", key="eta7_basic")
    
    if st.sidebar.button("Run Basic Model", key="btn_basic"):
        df_basic = model_basic_pama(C, MW, eta7_exp)
        st.plotly_chart(
            make_plot(df_basic, "Basic PAMA Viscosity vs Shear Rate", ["Viscosity (cP)"], ["Viscosity"]),
            use_container_width=True,
        )
        st.dataframe(df_basic.style.format("{:.3f}"))
        st.download_button(
            label="Download Data (CSV)",
            data=df_basic.to_csv(index=False),
            file_name="basic_pama.csv",
            mime="text/csv",
        )

elif model_choice == "PAMA with Temperature":
    st.header("PAMA with Temperature Effects")
    st.markdown("Explore how temperature influences viscosity in the PAMA model.")
    C = st.sidebar.slider("Concentration (g/L)", 0.1, 20.0, 2.0, 0.1, key="C_temp")
    MW = st.sidebar.slider("Molecular Weight (MDa)", 0.1, 50.0, 8.0, 0.1, key="MW_temp")
    eta7_exp = st.sidebar.number_input("η@7.3 experimental at 25°C (cP)", min_value=0.01, value=15.653, format="%.3f", key="eta7_temp")
    T_wanted_C = st.sidebar.slider("Target Temperature (°C)", 0, 100, 35, 1, key="T_temp")
    
    if st.sidebar.button("Run Temperature Model", key="btn_temp"):
        df_temp = model_pama_temperature(C, MW, eta7_exp, T_wanted_C)
        st.plotly_chart(
            make_plot(
                df_temp,
                "PAMA Temperature Model Viscosity",
                ["25°C Reference", f"{T_wanted_C}°C Model"],
                ["25°C", f"{T_wanted_C}°C"],
            ),
            use_container_width=True,
        )
        st.dataframe(df_temp.style.format("{:.3f}"))
        st.download_button(
            label="Download Data (CSV)",
            data=df_temp.to_csv(index=False),
            file_name="pama_temperature.csv",
            mime="text/csv",
        )

else:
    st.header("PAMA with Degradation Effects")
    st.markdown("Compare viscosity before and after polymer degradation.")
    C = st.sidebar.slider("Concentration (g/L)", 0.1, 20.0, 2.0, 0.1, key="C_deg")
    MW = st.sidebar.slider("Molecular Weight (MDa)", 0.1, 50.0, 8.0, 0.1, key="MW_deg")
    eta7_exp = st.sidebar.number_input("η@7.3 experimental (Polymer UD) (cP)", min_value=0.01, value=15.653, format="%.3f", key="eta7_deg")
    eta7_exp_D = st.sidebar.number_input("η@7.3 experimental (Polymer D) (cP)", min_value=0.01, value=7.354, format="%.3f", key="eta7_deg_D")
    
    if st.sidebar.button("Run Degradation Model", key="btn_deg"):
        df_deg = model_pama_degradation(C, MW, eta7_exp, eta7_exp_D)
        st.plotly_chart(
            make_plot(
                df_deg,
                "PAMA Degradation Model Viscosity",
                ["Polymer UD", "Polymer Degraded"],
                ["Polymer UD", "Degraded"],
            ),
            use_container_width=True,
        )
        st.dataframe(df_deg.style.format("{:.3f}"))
        st.download_button(
            label="Download Data (CSV)",
            data=df_deg.to_csv(index=False),
            file_name="pama_degradation.csv",
            mime="text/csv",
        )

# --- FOOTER ---
st.markdown("---")
st.markdown(
    "<p style='text-align:center; color:#666;'>Made by Eduar — Powered by Streamlit & Plotly</p>",
    unsafe_allow_html=True,
)
