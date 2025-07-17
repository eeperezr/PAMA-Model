import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.interpolate import interp1d

st.set_page_config(page_title="PAMA Models", layout="centered")
st.title("🧪 PAMA Method: Rheology Models")

model_choice = st.sidebar.radio(
    "Select Model",
    ("Model 1: Basic", "Model 2: With Temperature", "Model 3: With Degradation")
)

# ========= Model 1 =========
def model_basic_pama(C, MW, eta7_exp):
    Temp = 298
    eta_in = np.exp(-3.7188 + (578.919 / (-137.546 + Temp)))

    ceta = np.arange(0.1, 100.1, 0.1)
    eta_0 = 1 + ceta + 0.582 * ceta**2.009 + 0.022 * ceta**4
    n = 1 - (0.6187 - 0.5203 * np.exp(-0.08212 * ceta))
    l_d = 0.251 + 1.54 * MW * ceta / (C * Temp)
    l = l_d * (0.810 + 0.0230 * ceta**2.438)
    eta7 = 1 + (eta_0 - 1) * (1 + (l * 7.3)**2) ** ((n - 1) / 2)

    ceta_exp = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")(eta7_exp)

    eta_0C = 1 + ceta_exp + 0.582 * ceta_exp**2.009 + 0.022 * ceta_exp**4
    nC = 1 - (0.6187 - 0.5203 * np.exp(-0.08212 * ceta_exp))
    lC = (0.251 + 1.54 * MW * ceta_exp / (C * Temp)) * (0.810 + 0.0230 * ceta_exp**2.438)

    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 1000, 1)])
    etaC = eta_in + (eta_0C - eta_in) * (1 + (lC * shear)**2) ** ((nC - 1) / 2)

    df = pd.DataFrame({"shear": shear, "viscosity_cP": etaC})

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=shear, y=etaC, mode='lines', name="Viscosity"))
    fig.update_layout(
        xaxis=dict(type="log", title="Shear Rate (s⁻¹)"),
        yaxis=dict(type="log", title="Viscosity (cP)"),
        title="Model 1: Basic PAMA"
    )
    return df, fig


# ========= Model 2 =========
def model_pama_temperature(C, MW, eta7_exp, T_C):
    T = T_C + 273
    Temp = 298
    eta_in = np.exp(-3.7188 + (578.919 / (-137.546 + Temp)))

    ceta = np.arange(0.1, 100.1, 0.1)
    eta_0 = 1 + ceta + 0.582 * ceta**2.009 + 0.022 * ceta**4
    n = 1 - (0.6187 - 0.5203 * np.exp(-0.08212 * ceta))
    l_d = 0.251 + 1.54 * MW * ceta / (C * Temp)
    l = l_d * (0.810 + 0.0230 * ceta**2.438)
    eta7 = 1 + (eta_0 - 1) * (1 + (l * 7.3)**2) ** ((n - 1) / 2)

    ceta_exp = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")(eta7_exp)

    nC = 1 - (0.6187 - 0.5203 * np.exp(-0.08212 * ceta_exp))
    eta_0C = 1 + ceta_exp + 0.582 * ceta_exp**2.009 + 0.022 * ceta_exp**4
    lC = (0.251 + 1.54 * MW * ceta_exp / (C * Temp)) * (0.810 + 0.0230 * ceta_exp**2.438)
    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 1000, 1)])

    etaC_25 = eta_in + (eta_0C - eta_in) * (1 + (lC * shear)**2) ** ((nC - 1) / 2)
    eta_in_T = np.exp(-3.7188 + (578.919 / (-137.546 + T)))
    eta_0T = eta_in_T * (1 + ceta_exp + 0.582 * ceta_exp**2.009 + 0.022 * ceta_exp**4)
    lC_T = (0.251 + 1.54 * MW * ceta_exp / (C * T)) * (0.810 + 0.0230 * ceta_exp**2.438)
    etaC_T = eta_in_T + (eta_0T - eta_in_T) * (1 + (lC_T * shear)**2) ** ((nC - 1) / 2)

    df = pd.DataFrame({"shear": shear, "T25_C": etaC_25, f"T{T_C}_C": etaC_T})
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=shear, y=etaC_25, name="T = 25°C"))
    fig.add_trace(go.Scatter(x=shear, y=etaC_T, name=f"T = {T_C}°C"))
    fig.update_layout(
        xaxis=dict(type="log", title="Shear Rate (s⁻¹)"),
        yaxis=dict(type="log", title="Viscosity (cP)"),
        title="Model 2: Temperature Dependence"
    )
    return df, fig


# ========= Model 3 =========
def model_pama_degradation(C, MW, eta7_exp, eta7_exp_D):
    alpha = 0.763
    Temp = 298
    eta_in = np.exp(-3.7188 + (578.919 / (-137.546 + Temp)))

    ceta = np.arange(0.1, 100.1, 0.1)
    eta_0 = 1 + ceta + 0.582 * ceta**2.009 + 0.022 * ceta**4
    n = 1 - (0.6187 - 0.5203 * np.exp(-0.08212 * ceta))
    l_d = 0.251 + 1.54 * MW * ceta / (C * Temp)
    l = l_d * (0.810 + 0.0230 * ceta**2.438)
    eta7 = 1 + (eta_0 - 1) * (1 + (l * 7.3)**2) ** ((n - 1) / 2)

    ceta_exp = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")(eta7_exp)

    eta_0C = eta_in + ceta_exp + 0.582 * ceta_exp**2.009 + 0.022 * ceta_exp**4
    nC = 1 - (0.6187 - 0.5203 * np.exp(-0.08212 * ceta_exp))
    lC = (0.251 + 1.54 * MW * ceta_exp / (C * Temp)) * (0.810 + 0.0230 * ceta_exp**2.438)
    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 10000, 1)])
    etaC = eta_in + (eta_0C - eta_in) * (1 + (lC * shear)**2) ** ((nC - 1) / 2)

    # Degradation coefficients
    FD = np.log(eta7_exp_D / eta7_exp) / np.log(eta7_exp / eta7_exp_D)
    eta_0_D = eta_0C * np.exp(-0.1 * FD)
    n_D = nC * np.exp(-0.05 * FD)
    l_D = lC * np.exp(-0.1 * FD)

    etaD = eta_in + (eta_0_D - eta_in) * (1 + (l_D * shear)**2) ** ((n_D - 1) / 2)

    df = pd.DataFrame({"shear": shear, "Polymer_UD": etaC, "Polymer_D": etaD, "FD": FD})

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=shear, y=etaC, name="Undegraded"))
    fig.add_trace(go.Scatter(x=shear, y=etaD, name="Degraded"))
    fig.update_layout(
        xaxis=dict(type="log", title="Shear Rate (s⁻¹)"),
        yaxis=dict(type="log", title="Viscosity (cP)"),
        title="Model 3: Degradation Effect"
    )
    return df, fig

# ==== UI ====
if model_choice == "Model 1: Basic":
    st.subheader("Model 1: Basic PAMA")
    C = st.number_input("Concentration (g/L)", value=2.0)
    MW = st.number_input("Molecular Weight (MDa)", value=8.0)
    eta7_exp = st.number_input("Experimental η@7.3 (cP)", value=20.0)
    if st.button("Run Model 1"):
        df, fig = model_basic_pama(C, MW, eta7_exp)
        st.plotly_chart(fig)
        st.dataframe(df)
        st.download_button("Download CSV", df.to_csv(index=False), "model1_basic.csv")

elif model_choice == "Model 2: With Temperature":
    st.subheader("Model 2: Temperature Dependence")
    C = st.number_input("Concentration (g/L)", value=2.0)
    MW = st.number_input("Molecular Weight (MDa)", value=8.0)
    eta7_exp = st.number_input("Experimental η@7.3 (cP)", value=15.653)
    T = st.number_input("Target Temperature (°C)", value=35.0)
    if st.button("Run Model 2"):
        df, fig = model_pama_temperature(C, MW, eta7_exp, T)
        st.plotly_chart(fig)
        st.dataframe(df)
        st.download_button("Download CSV", df.to_csv(index=False), "model2_temp.csv")

elif model_choice == "Model 3: With Degradation":
    st.subheader("Model 3: Degradation Effect")
    C = st.number_input("Concentration (g/L)", value=2.0)
    MW = st.number_input("Molecular Weight (MDa)", value=8.0)
    eta7_exp = st.number_input("Undegraded η@7.3 (cP)", value=15.653)
    eta7_exp_D = st.number_input("Degraded ηD@7.3 (cP)", value=7.354)
    if st.button("Run Model 3"):
        df, fig = model_pama_degradation(C, MW, eta7_exp, eta7_exp_D)
        st.plotly_chart(fig)
        st.dataframe(df)
        st.download_button("Download CSV", df.to_csv(index=False), "model3_degradation.csv")
