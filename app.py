import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.interpolate import interp1d

# Page config and style
st.set_page_config(page_title="PAMA Method Models", layout="wide")

st.markdown(
    """
    <style>
    body {
        background-color: #ffffff;
        color: #222222;
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    }
    .stButton>button {
        background-color: #0a84ff;
        color: white;
        font-weight: bold;
        border-radius: 8px;
        padding: 0.5em 1.5em;
        transition: background-color 0.3s ease;
    }
    .stButton>button:hover {
        background-color: #006fd6;
    }
    .block-container {
        padding: 2rem 3rem 3rem 3rem;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("PAMA Method Models")
st.sidebar.header("Choose your model")

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
        xaxis=dict(type="log", title="Shear rate (s⁻¹)"),
        yaxis=dict(type="log", title="Viscosity (cP)"),
        title="Basic PAMA Model",
        template="plotly_white"
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

    # Calculate temperature dependent eta_0 and lC
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
    fig.add_trace(go.Scatter(x=shear, y=etaC_ref, mode='lines', name="T = 25°C"))
    fig.add_trace(go.Scatter(x=shear, y=etaC_temp, mode='lines', name=f"T = {T_wanted_C}°C"))
    fig.update_layout(
        xaxis=dict(type="log", title="Shear rate (s⁻¹)"),
        yaxis=dict(type="log", title="Viscosity (cP)"),
        title="PAMA Model with Temperature Dependence",
        template="plotly_white"
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

    parc2 = (n - 1) / 2
    parc2b = (n_FD - 1) / 2

    l_d = 0.251 + 1.54 * MW * ceta / (C * Temp)
    l_d_b = 0.251 + 1.54 * MW * cetab / (C * Temp)

    l = l_d * (0.810 + 0.0230 * ceta ** 2.438)
    lb = l_d_b * (0.810 + 0.0230 * cetab ** 2.438)

    parc1 = (l * 7.3) ** 2
    parc1b = (lb * 7.3) ** 2

    eta7 = 1 + (eta_0 - 1) * (1 + parc1) ** parc2
    eta7_FD = 1 + (eta_0_FD - 1) * (1 + parc1b) ** parc2b

    f = interp1d(eta7, ceta, bounds_error=False, fill_value="extrapolate")
    ceta_exp = f(eta7_exp)

    f_D = interp1d(eta7_FD, cetab, bounds_error=False, fill_value="extrapolate")
    cetab_exp = f_D(eta7_exp_D)

    eta_0C = 1 + ceta_exp + 0.582 * ceta_exp ** 2.009 + 0.022 * ceta_exp ** 4
    eta_0C_D = 1 + cetab_exp + 0.582 * cetab_exp ** 2.009 + 0.022 * cetab_exp ** 4

    parcC = -0.08212 * ceta_exp
    parcC_D = -0.08212 * cetab_exp

    nC = 1 - (0.6187 - 0.5203 * np.exp(parcC))
    nC_D = 1 - (0.6187 - 0.5203 * np.exp(parcC_D))

    parc2C = (nC - 1) / 2
    parc2C_D = (nC_D - 1) / 2

    l_dC = 0.251 + 1.54 * MW * ceta_exp / (C * Temp)
    l_dC_D = 0.251 + 1.54 * MW * cetab_exp / (C * Temp)

    lC = l_dC * (0.810 + 0.0230 * ceta_exp ** 2.438)
    lC_D = l_dC_D * (0.810 + 0.0230 * cetab_exp ** 2.438)

    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 1000, 1)])

    etaC = eta_in + (eta_0C - eta_in) * (1 + (lC * shear) ** 2) ** ((nC - 1) / 2)
    etaC_D = eta_in + (eta_0C_D - eta_in) * (1 + (lC_D * shear) ** 2) ** ((nC_D - 1) / 2)

    df = pd.DataFrame({"shear": shear, "eta": etaC, "eta_degraded": etaC_D})

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=shear, y=etaC, mode='lines', name="Undegraded"))
    fig.add_trace(go.Scatter(x=shear, y=etaC_D, mode='lines', name="Degraded"))
    fig.update_layout(
        xaxis=dict(type="log", title="Shear rate (s⁻¹)"),
        yaxis=dict(type="log", title="Viscosity (cP)"),
        title="PAMA Model with Degradation",
        template="plotly_white"
    )

    return df, fig


# -------- UI Inputs and Logic --------
if model_choice == "Basic PAMA":
    st.header("Basic PAMA Model Parameters")

    C = st.slider("Concentration (g/L)", 0.1, 20.0, 2.0, 0.1, key="C_basic")
    MW = st.slider("Molecular Weight (MDa)", 0.1, 50.0, 8.0, 0.1, key="MW_basic")
    eta7_exp = st.slider("Intrinsic Viscosity η₇ Experimental (cP)", 1.0, 100.0, 15.0, 0.5, key="eta7_basic")

    if st.button("Calculate Basic PAMA", key="btn_basic"):
        df_basic, fig_basic = model_basic_pama(C, MW, eta7_exp)
        st.plotly_chart(fig_basic, use_container_width=True)

elif model_choice == "PAMA with Temperature":
    st.header("PAMA Model with Temperature Dependence")

    C = st.slider("Concentration (g/L)", 0.1, 20.0, 2.0, 0.1, key="C_temp")
    MW = st.slider("Molecular Weight (MDa)", 0.1, 50.0, 8.0, 0.1, key="MW_temp")
    eta7_exp = st.slider("Intrinsic Viscosity η₇ Experimental (cP)", 1.0, 100.0, 15.0, 0.5, key="eta7_temp")
    T_wanted = st.slider("Temperature (°C)", 20, 80, 25, 1, key="T_temp")

    if st.button("Calculate PAMA with Temperature", key="btn_temp"):
        df_temp, fig_temp = model_pama_temperature(C, MW, eta7_exp, T_wanted)
        st.plotly_chart(fig_temp, use_container_width=True)

else:
    st.header("PAMA Model with Degradation")

    C = st.slider("Concentration (g/L)", 0.1, 20.0, 2.0, 0.1, key="C_deg")
    MW = st.slider("Molecular Weight (MDa)", 0.1, 50.0, 8.0, 0.1, key="MW_deg")
    eta7_exp = st.slider("Intrinsic Viscosity η₇ Experimental (cP)", 1.0, 100.0, 15.0, 0.5, key="eta7_deg")
    eta7_exp_D = st.slider("Degraded Intrinsic Viscosity η₇ Experimental (cP)", 1.0, 100.0, 10.0, 0.5, key="eta7_deg_D")

    if st.button("Calculate PAMA with Degradation", key="btn_deg"):
        df_deg, fig_deg = model_pama_degradation(C, MW, eta7_exp, eta7_exp_D)
        st.plotly_chart(fig_deg, use_container_width=True)
