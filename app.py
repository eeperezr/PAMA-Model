import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

# --- PAGE CONFIG ---
st.set_page_config(
    page_title="PAMA Rheology Models",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="collapsed",  # Collapse sidebar by default
)

# --- CUSTOM CSS for RepTate-like UI ---
st.markdown(
    """
    <style>
    body {
        background-color: #ffffff;
        color: #222222;
        font-family: 'Arial', sans-serif;
    }
    /* Sidebar styling (Table of Contents) */
    .css-1d391kg {
        width: 200px;
        background-color: #f0f0f0;
        border-right: 1px solid #ccc;
        padding: 10px;
        position: fixed;
        top: 60px;
        bottom: 0;
        overflow-y: auto;
    }
    /* Header styling */
    .css-1aumxhk {
        background-color: #1e90ff;
        color: white;
        padding: 10px;
        text-align: center;
    }
    /* Main content area */
    .block-container {
        margin-left: 220px;
        padding: 20px;
    }
    /* Buttons styling */
    .stButton>button {
        background-color: #1e90ff;
        color: white;
        font-weight: 600;
        padding: 0.5rem 1rem;
        border-radius: 5px;
        border: none;
        transition: background-color 0.3s ease;
    }
    .stButton>button:hover {
        background-color: #104e8b;
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
        border-radius: 5px;
        box-shadow: 0 0 10px rgba(0,0,0,0.1);
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# --- HEADER ---
st.markdown("<h1 style='background-color: #1e90ff; color: white; padding: 10px; text-align: center;'>PAMA Rheology Modeling App 🧪</h1>", unsafe_allow_html=True)
st.markdown("<p style='font-size: 16px; color: #444; text-align: center;'>Tool for analyzing polymer rheology models with interactive visualizations.</p>", unsafe_allow_html=True)

# --- SIDEBAR (Table of Contents) ---
with st.sidebar:
    st.markdown("<h3>Table of Contents</h3>", unsafe_allow_html=True)
    st.markdown("<ul><li><a href='#'>PAMA Documentation</a></li>", unsafe_allow_html=True)
    st.markdown("<ul><li><a href='#'>Contents</a></li>", unsafe_allow_html=True)
    st.markdown("<ul><li><a href='#'>About PAMA</a></li>", unsafe_allow_html=True)
    st.markdown("<li><a href='#'>Installation</a></li>", unsafe_allow_html=True)
    st.markdown("<li><a href='#'>User Manual</a></li>", unsafe_allow_html=True)
    st.markdown("<li><a href='#'>PAMA for Developers</a></li>", unsafe_allow_helper=True)
    st.markdown("<li><a href='#'>PAMA Contributors</a></li>", unsafe_allow_html=True)
    st.markdown("<li><a href='#'>Version History</a></li></ul>", unsafe_allow_html=True)
    st.markdown("<p>More Info:</p>", unsafe_allow_html=True)
    st.markdown("<ul><li><a href='https://example.com'>Source Code</a></li></ul>", unsafe_allow_html=True)
    st.markdown("<p>Authors:</p>", unsafe_allow_html=True)
    st.markdown("<ul><li>Jorge Ramirez (jorge.ramirez@upm.es)</li>", unsafe_allow_html=True)
    st.markdown("<li>Victor Boudara (victor.bc@gmail.com)</li></ul>", unsafe_allow_html=True)

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

    return pd.DataFrame({"shear": shear, "Viscosity (cP)": etaC})

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

    return pd.DataFrame({"shear": shear, "25°C Reference": etaC_ref, f"{T_wanted_C}°C Model": etaC_temp})

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

    return pd.DataFrame({"shear": shear, "Polymer UD": etaC, "Polymer Degraded": etaC_D})

# --- MAIN CONTENT ---
st.markdown("<h2>PAMA Rheology Models</h2>", unsafe_allow_html=True)
st.markdown("<p>PAMA (Polymer Analysis and Modeling Application) is a tool for analyzing polymer rheology models with interactive visualizations.</p>", unsafe_allow_html=True)

# Model selection and parameters
model_choice = st.selectbox("Select Model", ("Basic PAMA", "PAMA with Temperature", "PAMA with Degradation"))

if model_choice == "Basic PAMA":
    st.markdown("<h3>Basic PAMA Model</h3>", unsafe_allow_html=True)
    st.markdown("<p>Input parameters to calculate viscosity using the basic PAMA model.</p>", unsafe_allow_html=True)
    col1, col2 = st.columns([1, 2])
    with col1:
        C = st.slider("Concentration (g/L)", 0.1, 20.0, 2.0, 0.1)
        MW = st.slider("Molecular Weight (MDa)", 0.1, 50.0, 8.0, 0.1)
        eta7_exp = st.number_input("η@7.3 experimental (cP)", min_value=0.01, value=20.0, format="%.3f")
        if st.button("Run Basic Model"):
            df_basic = model_basic_pama(C, MW, eta7_exp)
            with col2:
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
    st.markdown("<h3>PAMA with Temperature Effects</h3>", unsafe_allow_html=True)
    st.markdown("<p>Explore how temperature influences viscosity in the PAMA model.</p>", unsafe_allow_html=True)
    col1, col2 = st.columns([1, 2])
    with col1:
        C = st.slider("Concentration (g/L)", 0.1, 20.0, 2.0, 0.1)
        MW = st.slider("Molecular Weight (MDa)", 0.1, 50.0, 8.0, 0.1)
        eta7_exp = st.number_input("η@7.3 experimental at 25°C (cP)", min_value=0.01, value=15.653, format="%.3f")
        T_wanted_C = st.slider("Target Temperature (°C)", 0, 100, 35, 1)
        if st.button("Run Temperature Model"):
            df_temp = model_pama_temperature(C, MW, eta7_exp, T_wanted_C)
            with col2:
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
    st.markdown("<h3>PAMA with Degradation Effects</h3>", unsafe_allow_html=True)
    st.markdown("<p>Compare viscosity before and after polymer degradation.</p>", unsafe_allow_html=True)
    col1, col2 = st.columns([1, 2])
    with col1:
        C = st.slider("Concentration (g/L)", 0.1, 20.0, 2.0, 0.1)
        MW = st.slider("Molecular Weight (MDa)", 0.1, 50.0, 8.0, 0.1)
        eta7_exp = st.number_input("η@7.3 experimental (Polymer UD) (cP)", min_value=0.01, value=15.653, format="%.3f")
        eta7_exp_D = st.number_input("η@7.3 experimental (Polymer D) (cP)", min_value=0.01, value=7.354, format="%.3f")
        if st.button("Run Degradation Model"):
            df_deg = model_pama_degradation(C, MW, eta7_exp, eta7_exp_D)
            with col2:
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
st.markdown("<hr>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center; color: #666;'>Made by Eduar — Powered by Streamlit & Plotly</p>", unsafe_allow_html=True)
