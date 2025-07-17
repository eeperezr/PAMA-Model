import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.interpolate import interp1d

st.title("PAMA Method Models")

model_choice = st.sidebar.radio(
    "Select Model",
    ("Basic PAMA", "PAMA with Temperature", "PAMA with Degradation")
)

# ========== Model 1: Basic PAMA ==========
def model_basic_pama(C, MW, eta7_exp):
    Temp = 298  # constant 25 C in Kelvin
    eta_in = np.exp(-3.7188 + (578.919 / (-137.546 + Temp)) + (0.0608 * 0 ** 1.3533))

    ceta = np.arange(0.1, 100.1, 0.1)
    eta_0 = 1 + ceta + 0.582 * ceta ** 2.009 + 0.022 * ceta ** 4
    parc = -0.08212 * ceta
    n = 1 - (0.6187 - 0.5203 * np.exp(parc))
    parc2 = (n - 1) / 2
    l_d = 0.251 + 1.54 * MW * ceta / (C * Temp)
    l = l_d * (0.810 + 0.0230 * ceta ** 2.438)
    parc1 = (l * 7.3) ** 2
    eta7 = 1 + (eta_0 - 1) * (1 + parc1) ** parc2

    # interpolate ceta at experimental eta7
    f = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")
    ceta_exp = f(eta7_exp)

    # Carreau-Yassuda coefficients
    eta_0C = 1 + ceta_exp + 0.582 * ceta_exp ** 2.009 + 0.022 * ceta_exp ** 4
    parcC = -0.08212 * ceta_exp
    nC = 1 - (0.6187 - 0.5203 * np.exp(parcC))
    l_dC = 0.251 + 1.54 * MW * ceta_exp / (C * Temp)
    lC = l_dC * (0.810 + 0.0230 * ceta_exp ** 2.438)

    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 1000, 1)])
    etaC = eta_in + (eta_0C - eta_in) * (1 + (lC * shear) ** 2) ** ((nC - 1) / 2)

    df = pd.DataFrame({"shear": shear, "eta": etaC})

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=shear, y=etaC, mode='lines', name="Viscosity"))
    fig.update_layout(
        xaxis=dict(type="log", title="Shear rate (s^-1)"),
        yaxis=dict(type="log", title="Viscosity (cP)"),
        title="Basic PAMA Model"
    )

    return df, fig


# ========== Model 2: PAMA with Temperature ==========
def model_pama_temperature(C, MW, eta7_exp, T_wanted_C):
    T_wanted = T_wanted_C + 273  # Convert to Kelvin
    Temp = 298  # Reference temp 25 C in Kelvin
    eta_in = np.exp(-3.7188 + (578.919 / (-137.546 + Temp)) + (0.0608 * 0 ** 1.3533))

    ceta = np.arange(0.1, 100.1, 0.1)
    eta_0 = 1 + ceta + 0.582 * ceta ** 2.009 + 0.022 * ceta ** 4
    parc = -0.08212 * ceta
    n = 1 - (0.6187 - 0.5203 * np.exp(parc))
    parc2 = (n - 1) / 2
    l_d = 0.251 + 1.54 * MW * ceta / (C * Temp)
    l = l_d * (0.810 + 0.0230 * ceta ** 2.438)
    parc1 = (l * 7.3) ** 2
    eta7 = 1 + (eta_0 - 1) * (1 + parc1) ** parc2

    f = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")
    ceta_exp = f(eta7_exp)

    # Carreau coefficients at reference Temp
    eta_0C = 1 + ceta_exp + 0.582 * ceta_exp ** 2.009 + 0.022 * ceta_exp ** 4
    parcC = -0.08212 * ceta_exp
    nC = 1 - (0.6187 - 0.5203 * np.exp(parcC))
    parc2C = (nC - 1) / 2
    l_dC = 0.251 + 1.54 * MW * ceta_exp / (C * Temp)
    lC = l_dC * (0.810 + 0.0230 * ceta_exp ** 2.438)

    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 1000, 1)])

    etaC_ref = eta_in + (eta_0C - eta_in) * (1 + (lC * shear) ** 2) ** ((nC - 1) / 2)

    # Temperature dependent eta_in
    eta_in_T = np.exp(-3.7188 + (578.919 / (-137.546 + T_wanted)) + (0.0608 * 0 ** 1.3533))

    # Calculate temperature dependent eta_0 and lC (simplified, detailed from R can be added)
    eta_0T = eta_in_T * (1 + ceta_exp + 0.582 * ceta_exp ** 2.009 + 0.022 * ceta_exp ** 4)
    l_dC_T = 0.251 + 1.54 * MW * ceta_exp / (C * T_wanted)
    lC_T = l_dC_T * (0.810 + 0.0230 * ceta_exp ** 2.438)

    etaC_temp = eta_in_T + (eta_0T - eta_in_T) * (1 + (lC_T * shear) ** 2) ** ((nC - 1) / 2)

    df = pd.DataFrame({
        "shear": shear,
        "eta_ref": etaC_ref,
        "eta_temp": etaC_temp
    })

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=shear, y=etaC_ref, mode='lines', name="T = 25.0°C"))
    fig.add_trace(go.Scatter(x=shear, y=etaC_temp, mode='lines', name=f"T = {T_wanted_C}°C"))
    fig.update_layout(
        xaxis=dict(type="log", title="Shear rate (s^-1)"),
        yaxis=dict(type="log", title="Viscosity (cP)"),
        title="PAMA Model with Temperature Dependence"
    )

    return df, fig


# ========== Model 3: PAMA with Degradation ==========
def model_pama_degradation(C, MW, eta7_exp, eta7_exp_D):
    alpha = 0.763
    Temp = 298
    eta_in = np.exp(-3.7188 + (578.919 / (-137.546 + Temp)) + (0.0608 * 0 ** 1.3533))

    c_med = C * 0.5
    ceta = np.arange(0.1, 100.1, 0.1)
    eta = ceta / c_med

    eta_0 = 1 + ceta + 0.582 * ceta ** 2.009 + 0.022 * ceta ** 4
    etab = eta / (2 ** alpha)
    cetab = etab * c_med
    eta_0_FD = 1 + cetab + 0.582 * cetab ** 2.009 + 0.022 * cetab ** 4

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

    l = l_d * (0.810 + 0.0230 * ceta ** 2.438)
    lb = l_db * (0.810 + 0.0230 * cetab ** 2.438)

    omega = np.log(lb) - np.log(l)

    parc1 = l * 7.3
    parc1b = lb * 7.3

    eta7 = 1 + (eta_0 - 1) * (1 + parc1) ** parc2
    eta7b = 1 + (eta_0_FD - 1) * (1 + parc1b) ** parc2b

    rho = np.log(eta7b) - np.log(eta7)

    # interpolate ceta at experimental eta7
    f_ceta = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")
    ceta_exp = f_ceta(eta7_exp)

    y_beta = interp1d(ceta, beta, bounds_error=False, fill_value="extrapolate")
    beta_find = y_beta(ceta_exp)

    y_gamma = interp1d(ceta, gama, bounds_error=False, fill_value="extrapolate")
    gama_find = y_gamma(ceta_exp)

    y_omega = interp1d(ceta, omega, bounds_error=False, fill_value="extrapolate")
    omega_find = y_omega(ceta_exp)

    y_rho = interp1d(ceta, rho, bounds_error=False, fill_value="extrapolate")
    rho_find = y_rho(ceta_exp)

    eta_0C = eta_in + ceta_exp + 0.582 * ceta_exp ** 2.009 + 0.022 * ceta_exp ** 4
    parcC = -0.08212 * ceta_exp
    nC = 1 - (0.6187 - 0.5203 * np.exp(parcC))
    parc2C = (nC - 1) / 2
    l_dC = 0.251 + 1.54 * MW * ceta_exp / (C * Temp)
    lC = l_dC * (0.810 + 0.0230 * ceta_exp ** 2.438)

    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 10000, 1)])

    etaC = eta_in + (eta_0C - eta_in) * (1 + (lC * shear) ** 2) ** ((nC - 1) / 2)

    FD = np.log(eta7_exp_D / eta7_exp) / rho_find

    eta_0_D = np.exp(beta_find * FD) * eta_0C
    n_D = np.exp(gama_find * FD) * nC
    la_D = np.exp(omega_find * FD) * lC

    eta_D = eta_in + (eta_0_D - eta_in) * (1 + (la_D * shear) ** 2) ** ((n_D - 1) / 2)

    df = pd.DataFrame({
        "shear": shear,
        "Polymer_UD_cP": etaC,
        "Polymer_D_cP": eta_D,
        "FD": FD
    })

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=shear, y=etaC, mode='lines', name="Polymer UD, cP"))
    fig.add_trace(go.Scatter(x=shear, y=eta_D, mode='lines', name="Polymer D, cP"))
    fig.update_layout(
        xaxis=dict(type="log", title="Shear, s^-1"),
        yaxis=dict(type="log", title="η, cP"),
        title="PAMA Model with Degradation"
    )

    return df, fig


# ==== UI ====

if model_choice == "Basic PAMA":
    st.header("Basic PAMA Model")
    C = st.number_input("Concentration, g/l", value=2.0)
    MW = st.number_input("Molecular Weight, MDa", value=8.0)
    eta7_exp = st.number_input("Experimental η@7.3, cP", value=20.0)

    if st.button("Run Basic Model"):
        df, fig = model_basic_pama(C, MW, eta7_exp)
        st.plotly_chart(fig)
        st.dataframe(df)
        st.download_button("Download Data", df.to_csv(index=False), "basic_pama.csv")

elif model_choice == "PAMA with Temperature":
    st.header("PAMA Model with Temperature")
    C = st.number_input("Concentration, g/l", value=2.0)
    MW = st.number_input("Molecular Weight, MDa", value=8.0)
    eta7_exp = st.number_input("Experimental η@7.3, cP", value=15.653)
    T_wanted_C = st.number_input("Temperature, °C", value=35.0)

    if st.button("Run Temperature Model"):
        df, fig = model_pama_temperature(C, MW, eta7_exp, T_wanted_C)
        st.plotly_chart(fig)
        st.dataframe(df)
        st.download_button("Download Data", df.to_csv(index=False), "pama_temperature.csv")

elif model_choice == "PAMA with Degradation":
    st.header("PAMA Model with Degradation")
    C = st.number_input("Concentration, g/l", value=2.0)
    MW = st.number_input("Molecular Weight, MDa", value=8.0)
    eta7_exp = st.number_input("Experimental η@7.3, cP", value=15.653)
    eta7_exp_D = st.number_input("Degraded solution ηD@7.3, cP", value=7.354)

    if st.button("Run Degradation Model"):
        df, fig = model_pama_degradation(C, MW, eta7_exp, eta7_exp_D)
        st.plotly_chart(fig)
        st.dataframe(df)
        st.download_button("Download Data", df.to_csv(index=False), "pama_degradation.csv")
