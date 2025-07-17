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
)

# --- SIDEBAR (Parameters and Model Selection) ---
with st.sidebar:
    st.title("PAMA Controls")
    st.markdown("### Navigation")
    st.markdown("""
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

    model = st.selectbox("Choose Model", ["Basic PAMA", "PAMA with Temperature", "PAMA with Degradation"])

    if 'conc' not in st.session_state:
        st.session_state.conc = 2.0
    if 'mw' not in st.session_state:
        st.session_state.mw = 8.0
    if 'eta7' not in st.session_state:
        st.session_state.eta7 = 20.0
    if 'temp' not in st.session_state:
        st.session_state.temp = 35
    if 'eta7d' not in st.session_state:
        st.session_state.eta7d = 7.354

    st.session_state.conc = st.slider("Concentration (g/L)", 0.1, 20.0, 2.0, 0.1)
    st.session_state.mw = st.slider("Molecular Weight (MDa)", 0.1, 50.0, 8.0, 0.1)
    st.session_state.eta7 = st.number_input("η@7.3 (cP)", min_value=0.01, value=20.0, format="%.3f")

    if model == "PAMA with Temperature":
        st.session_state.temp = st.slider("Temperature (°C)", 0, 100, 35, 1)
    elif model == "PAMA with Degradation":
        st.session_state.eta7d = st.number_input("η@7.3 (Polymer D) (cP)", min_value=0.01, value=7.354, format="%.3f")

    if st.button("Run Model"):
        if model == "Basic PAMA":
            data = model_basic_pama(st.session_state.conc, st.session_state.mw, st.session_state.eta7)
        elif model == "PAMA with Temperature":
            data = model_pama_temperature(st.session_state.conc, st.session_state.mw, st.session_state.eta7, st.session_state.temp)
        elif model == "PAMA with Degradation":
            data = model_pama_degradation(st.session_state.conc, st.session_state.mw, st.session_state.eta7, st.session_state.eta7d)
        
        st.session_state.data = data
        st.session_state.model = model

# --- MAIN CONTENT (Graphs and Fancy Stuff) ---
st.title("PAMA Rheology Models")
st.markdown("### Visualize Polymer Behavior")
st.markdown("Explore detailed rheology insights with interactive plots and data.")

with st.container():
    if 'data' in st.session_state and 'model' in st.session_state:
        data = st.session_state.data
        model = st.session_state.model

        # --- CHART ---
        traces = []
        columns = []
        if model == "Basic PAMA":
            traces.append(go.Scatter(x=data["shear"], y=data["Viscosity (cP)"], name="Viscosity", mode="lines+markers", line=dict(color="#1e90ff")))
            columns = ["shear", "Viscosity (cP)"]
        elif model == "PAMA with Temperature":
            traces.append(go.Scatter(x=data["shear"], y=data["25°C Reference"], name="25°C", mode="lines+markers", line=dict(color="#1e90ff")))
            traces.append(go.Scatter(x=data["shear"], y=data[list(data.keys())[2]], name=list(data.keys())[2], mode="lines+markers", line=dict(color="#4169e1")))
            columns = ["shear", "25°C Reference", list(data.keys())[2]]
        elif model == "PAMA with Degradation":
            traces.append(go.Scatter(x=data["shear"], y=data["Polymer UD"], name="Polymer UD", mode="lines+markers", line=dict(color="#1e90ff")))
            traces.append(go.Scatter(x=data["shear"], y=data["Polymer Degraded"], name="Polymer Degraded", mode="lines+markers", line=dict(color="#ff4500")))
            columns = ["shear", "Polymer UD", "Polymer Degraded"]

        fig = go.Figure(
            data=traces,
            layout=go.Layout(
                title=f"{model} Analysis",
                xaxis={"title": "Shear Rate (s⁻¹)", "type": "log", "gridcolor": "#e0e0e0"},
                yaxis={"title": "Viscosity (cP)", "type": "log", "gridcolor": "#e0e0e0"},
                template="plotly_white",
                legend={"orientation": "h", "yanchor": "bottom", "y": 1.1, "xanchor": "center", "x": 0.5"},
                margin={"t": 60, "b": 50, "l": 50, "r": 50"},
                plot_bgcolor="#f9f9f9",
                paper_bgcolor="#f9f9f9",
                hovermode="x unified"
            )
        )
        st.plotly_chart(fig, use_container_width=True)

        # --- TABLE ---
        df = pd.DataFrame({col: data[col] for col in columns})
        st.table(df.style.set_properties(**{'text-align': 'center', 'border': '1px solid #e0e0e0'}))

        # --- DOWNLOAD BUTTON ---
        csv = df.to_csv(index=False)
        st.download_button(
            label="Download Data",
            data=csv,
            file_name=f"{model.replace(' ', '_').lower()}_data.csv",
            mime="text/csv",
            help="Export the current data as a CSV file"
        )

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
