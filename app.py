import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.interpolate import interp1d

st.set_page_config(page_title="PAMA Method Models", layout="wide")

st.title("PAMA Method – Rheological Models")
st.markdown("""
This app simulates and visualizes **3 different rheological models** for HPAM solutions:

1. **Basic Model**
2. **Model with Temperature Adjustment**
3. **Model with Polymer Degradation**

Developed from the methodology by *E. Pérez, D. Alviso, E. Manrique, G. Artana*.
""")

# Select Model
model_choice = st.sidebar.radio(
    "Choose a model:",
    ["1. Basic", "2. With Temperature", "3. With Degradation"]
)

# Shared Inputs
C = st.sidebar.number_input("Concentration (g/L)", value=2.0)
MW = st.sidebar.number_input("Molecular Weight (MDa)", value=8.0)
eta7_exp = st.sidebar.number_input("Experimental η @ 7.3 s⁻¹ (cP)", value=20.0)

if model_choice == "3. With Degradation":
    eta7_exp_D = st.sidebar.number_input("Degraded ηD @ 7.3 s⁻¹ (cP)", value=10.0)

# Constants
Temp = 298
eta_in = np.exp((-3.7188 + (578.919 / (-137.546 + Temp)) + (0.0608 * (0 + 0) ** 1.3533)))
shear = np.concatenate((np.arange(0.01, 1, 0.01), np.arange(1, 10000, 1)))

# === Model Calculations ===

def base_model():
    ceta = np.arange(0.1, 100, 0.1)
    eta_0 = 1 + ceta + 0.582 * ceta**2.009 + 0.022 * ceta**4
    parc = -0.08212 * ceta
    n = 1 - (0.6187 - 0.5203 * np.exp(parc))
    parc2 = (n - 1) / 2
    l_d = 0.251 + 1.54 * MW * ceta / (C * Temp)
    l = l_d * (0.810 + 0.0230 * ceta**2.438)
    eta7 = 1 + (eta_0 - 1) * (1 + (l * 7.3))**parc2
    interp = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")
    ceta_exp = interp(eta7_exp)
    eta_0C = 1 + ceta_exp + 0.582 * ceta_exp**2.009 + 0.022 * ceta_exp**4
    parcC = -0.08212 * ceta_exp
    nC = 1 - (0.6187 - 0.5203 * np.exp(parcC))
    parc2C = (nC - 1) / 2
    l_dC = 0.251 + 1.54 * MW * ceta_exp / (C * Temp)
    lC = l_dC * (0.810 + 0.0230 * ceta_exp**2.438)
    etaC = eta_in + (eta_0C - eta_in) * (1 + (lC * shear)**2)**((nC - 1)/2)
    return shear, etaC, None

def model_with_temperature():
    alpha = 0.763
    c_med = C * 0.5
    ceta = np.arange(0.1, 100, 0.1)
    eta = ceta / c_med
    eta_0 = 1 + ceta + 0.582 * ceta**2.009 + 0.022 * ceta**4
    etab = eta / (2 ** alpha)
    cetab = etab * c_med
    eta_0_FD = 1 + cetab + 0.582 * cetab**2.009 + 0.022 * cetab**4
    beta = np.log(eta_0_FD) - np.log(eta_0)
    parc = -0.08212 * ceta
    parcb = -0.08212 * cetab
    n = 1 - (0.6187 - 0.5203 * np.exp(parc))
    n_FD = 1 - (0.6187 - 0.5203 * np.exp(parcb))
    gama = (n / n_FD) - 1
    parc2 = (n - 1) / 2
    parc2b = (n_FD - 1) / 2
    l_d = 0.251 + 1.54 * MW * ceta / (C * Temp)
    l_db = 0.251 + 1.54 * MW * cetab / (c_med * Temp)
    l = l_d * (0.810 + 0.0230 * ceta**2.438)
    lb = l_db * (0.810 + 0.0230 * cetab**2.438)
    omega = np.log(lb) - np.log(l)
    eta7 = 1 + (eta_0 - 1) * (1 + l * 7.3)**parc2
    eta7b = 1 + (eta_0_FD - 1) * (1 + lb * 7.3)**parc2b
    rho = np.log(eta7b) - np.log(eta7)
    interp = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")
    ceta_find = interp(eta7_exp)
    eta_0C = np.exp(interp1d(ceta, np.log(eta_0))(ceta_find))
    beta_find = interp1d(ceta, beta)(ceta_find)
    gama_find = interp1d(ceta, gama)(ceta_find)
    omega_find = interp1d(ceta, omega)(ceta_find)
    rho_find = interp1d(ceta, rho)(ceta_find)
    FD = 1.0  # assume temperature-induced shift
    eta_0_D = np.exp(beta_find * FD) * eta_0C
    nC = 1 - (0.6187 - 0.5203 * np.exp(-0.08212 * ceta_find))
    n_D = np.exp(gama_find * FD) * nC
    lC = (0.251 + 1.54 * MW * ceta_find / (C * Temp)) * (0.810 + 0.0230 * ceta_find**2.438)
    la_D = np.exp(omega_find * FD) * lC
    etaC = eta_in + (eta_0C - eta_in) * (1 + (lC * shear)**2)**((nC - 1)/2)
    etaD = eta_in + (eta_0_D - eta_in) * (1 + (la_D * shear)**2)**((n_D - 1)/2)
    return shear, etaC, etaD

def model_with_degradation():
    shear, etaC, etaD = model_with_temperature()
    FD = np.log(eta7_exp_D / eta7_exp) / np.log(etaD[0] / etaC[0])
    etaD = etaC * np.exp(FD * np.log(etaD / etaC))
    return shear, etaC, etaD

# Run selected model
if model_choice == "1. Basic":
    shear, etaC, _ = base_model()
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=shear, y=etaC, mode='lines', name="Viscosity η"))
elif model_choice == "2. With Temperature":
    shear, etaC, etaD = model_with_temperature()
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=shear, y=etaC, mode='lines', name="Undegraded"))
    fig.add_trace(go.Scatter(x=shear, y=etaD, mode='lines', name="Temp. Shifted"))
elif model_choice == "3. With Degradation":
    shear, etaC, etaD = model_with_degradation()
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=shear, y=etaC, mode='lines', name="Undegraded"))
    fig.add_trace(go.Scatter(x=shear, y=etaD, mode='lines', name="Degraded"))

# Plot
fig.update_layout(
    title="Carreau-Yasuda Rheological Curve",
    xaxis=dict(title="Shear Rate (1/s)", type="log"),
    yaxis=dict(title="Viscosity (cP)", type="log"),
    template="plotly_white",
    height=600
)
st.plotly_chart(fig, use_container_width=True)

# Download
df = pd.DataFrame({'shear_rate (1/s)': shear, 'etaC (cP)': etaC})
if etaD is not None:
    df['etaD (cP)'] = etaD

csv = df.to_csv(index=False)
st.download_button("Download CSV", csv, "pama_model_output.csv")
