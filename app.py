import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.interpolate import interp1d

# === PAGE CONFIG ===
st.set_page_config(page_title="PAMA Rheology Models", layout="wide")

# === TITLE ===
st.markdown("<h1 style='text-align: center;'>🧪 PAMA Rheology Modeling App</h1>", unsafe_allow_html=True)
st.markdown("<p style='text-align: center;'>Select and run different models to analyze polymer rheology under various conditions.</p>", unsafe_allow_html=True)
st.markdown("---")

# === SIDEBAR ===
with st.sidebar:
    st.header("Model Settings")
    model_choice = st.radio("Choose a Model", [
        "Basic PAMA",
        "PAMA with Temperature",
        "PAMA with Degradation"
    ])

    st.markdown("Made with ❤️ using Streamlit + Plotly")

# === COMMON FUNCTION FOR PLOTLY ===
def make_plot(df, title, ycols, labels):
    fig = go.Figure()
    for ycol, label in zip(ycols, labels):
        fig.add_trace(go.Scatter(x=df["shear"], y=df[ycol], mode='lines', name=label))
    fig.update_layout(
        xaxis=dict(type="log", title="Shear rate (s⁻¹)"),
        yaxis=dict(type="log", title="Viscosity (cP)"),
        template="plotly_white",
        title=title,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5)
    )
    return fig

# === MODEL 1: Basic PAMA ===
def model_basic_pama(C, MW, eta7_exp):
    Temp = 298
    eta_in = np.exp(-3.7188 + (578.919 / (-137.546 + Temp)))

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

    eta_0C = 1 + ceta_exp + 0.582 * ceta_exp ** 2.009 + 0.022 * ceta_exp ** 4
    parcC = -0.08212 * ceta_exp
    nC = 1 - (0.6187 - 0.5203 * np.exp(parcC))
    l_dC = 0.251 + 1.54 * MW * ceta_exp / (C * Temp)
    lC = l_dC * (0.810 + 0.0230 * ceta_exp ** 2.438)

    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 1000, 1)])
    etaC = eta_in + (eta_0C - eta_in) * (1 + (lC * shear) ** 2) ** ((nC - 1) / 2)

    df = pd.DataFrame({"shear": shear, "eta": etaC})
    return df

# === MODEL 2: PAMA with Temperature ===
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
    l = l_d * (0.810 + 0.0230 * ceta ** 2.438)
    eta7 = 1 + (eta_0 - 1) * (1 + (l * 7.3)**2) ** parc2

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
    etaC_ref = eta_in + (eta_0C - eta_in) * (1 + (lC * shear)**2) ** ((nC - 1) / 2)
    etaC_temp = eta_in_T + (eta_0C - eta_in_T) * (1 + (lC_T * shear)**2) ** ((nC - 1) / 2)

    df = pd.DataFrame({"shear": shear, "eta_ref": etaC_ref, "eta_temp": etaC_temp})
    return df

# === MODEL 3: PAMA with Degradation ===
def model_pama_degradation(C, MW, eta7_exp, eta7_exp_D):
    Temp = 298
    eta_in = np.exp(-3.7188 + (578.919 / (-137.546 + Temp)))

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
    ceta_exp_D = f(eta7_exp_D)

    # Polymer UD parameters
    eta_0C = 1 + ceta_exp + 0.582 * ceta_exp ** 2.009 + 0.022 * ceta_exp ** 4
    parcC = -0.08212 * ceta_exp
    nC = 1 - (0.6187 - 0.5203 * np.exp(parcC))
    l_dC = 0.251 + 1.54 * MW * ceta_exp / (C * Temp)
    lC = l_dC * (0.810 + 0.0230 * ceta_exp ** 2.438)

    # Polymer D parameters
    eta_0D = 1 + ceta_exp_D + 0.582 * ceta_exp_D ** 2.009 + 0.022 * ceta_exp_D ** 4
    parcD = -0.08212 * ceta_exp_D
    nD = 1 - (0.6187 - 0.5203 * np.exp(parcD))
    l_dD = 0.251 + 1.54 * MW * ceta_exp_D / (C * Temp)
    lD = l_dD * (0.810 + 0.0230 * ceta_exp_D ** 2.438)

    shear = np.concatenate([np.arange(0.01, 1, 0.01), np.arange(1, 1000, 1)])
    Polymer_UD_cP = eta_in + (eta_0C - eta_in) * (1 + (lC * shear) ** 2) ** ((nC - 1) / 2)
    Polymer_D_cP = eta_in + (eta_0D - eta_in) * (1 + (lD * shear) ** 2) ** ((nD - 1) / 2)

    df = pd.DataFrame({
        "shear": shear,
        "Polymer_UD_cP": Polymer_UD_cP,
        "Polymer_D_cP": Polymer_D_cP
    })
    return df

# === MAIN UI ===
if model_choice == "Basic PAMA":
    col1, col2 = st.columns([1, 2])
    with col1:
        C = st.number_input("Concentration (g/L)", value=2.0, min_value=0.01, format="%.3f")
        MW = st.number_input("Molecular Weight (MDa)", value=8.0, min_value=0.01, format="%.3f")
        eta7_exp = st.number_input("η@7.3 experimental (cP)", value=20.0, min_value=0.01, format="%.3f")
        if st.button("Run Model"):
            df = model_basic_pama(C, MW, eta7_exp)
            fig = make_plot(df, "Basic PAMA Model", ["eta"], ["Viscosity"])
            st.plotly_chart(fig, use_container_width=True)
            st.dataframe(df)
            st.download_button("📥 Download CSV", df.to_csv(index=False), "basic_pama.csv")

elif model_choice == "PAMA with Temperature":
    col1, col2 = st.columns([1, 2])
    with col1:
        C = st.number_input("Concentration (g/L)", value=2.0, min_value=0.01, format="%.3f")
        MW = st.number_input("Molecular Weight (MDa)", value=8.0, min_value=0.01, format="%.3f")
        eta7_exp = st.number_input("η@7.3 experimental (cP)", value=15.653, min_value=0.01, format="%.3f")
        T_wanted_C = st.number_input("Target Temp (°C)", value=35.0, format="%.1f")
        if st.button("Run Temperature Model"):
            df = model_pama_temperature(C, MW, eta7_exp, T_wanted_C)
            fig = make_plot(df, "PAMA Temperature Model", ["eta_ref", "eta_temp"], ["T = 25°C", f"T = {T_wanted_C}°C"])
            st.plotly_chart(fig, use_container_width=True)
            st.dataframe(df)
            st.download_button("📥 Download CSV", df.to_csv(index=False), "pama_temperature.csv")

elif model_choice == "PAMA with Degradation":
    col1, col2 = st.columns([1, 2])
    with col1:
        C = st.number_input("Concentration (g/L)", value=2.0, min_value=0.01, format="%.3f")
        MW = st.number_input("Molecular Weight (MDa)", value=8.0, min_value=0.01, format="%.3f")
        eta7_exp = st.number_input("η@7.3 experimental (cP)", value=15.653, min_value=0.01, format="%.3f")
        eta7_exp_D = st.number_input("η@7.3 degraded (cP)", value=7.354, min_value=0.01, format="%.3f")
        if st.button("Run Degradation Model"):
            df = model_pama_degradation(C, MW, eta7_exp, eta7_exp_D)
            fig = make_plot(df, "PAMA Degradation Model", ["Polymer_UD_cP", "Polymer_D_cP"], ["Polymer UD", "Polymer D"])
            st.plotly_chart(fig, use_container_width=True)
            st.dataframe(df)
            st.download_button("📥 Download CSV", df.to_csv(index=False), "pama_degradation.csv")
