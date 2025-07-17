import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy.interpolate import interp1d

# Page config
st.set_page_config(page_title="PAMA Viscosity Models", layout="wide")

# Theme switcher (manual override for colors)
theme = st.sidebar.selectbox("🌙 Theme", ["Light", "Dark"])
if theme == "Dark":
    bg_color = "#1e1e1e"
    text_color = "#ffffff"
else:
    bg_color = "#f5f8fc"
    text_color = "#222222"

st.markdown(f"""
    <style>
    body {{
        background-color: {bg_color};
        color: {text_color};
        font-family: 'Segoe UI', sans-serif;
    }}
    .stButton>button {{
        background-color: #005bbb;
        color: white;
        font-weight: bold;
        border-radius: 6px;
        padding: 0.5em 1.2em;
    }}
    .stButton>button:hover {{
        background-color: #003e87;
    }}
    </style>
""", unsafe_allow_html=True)


# ---------- MODEL FUNCTIONS (paste unchanged) ----------
# Paste the functions:
# model_basic_pama()
# model_pama_temperature()
# model_pama_degradation()
# (from your original script—they stay the same)

# --- UI Tabs ---
st.title("🧪 PAMA Viscosity Modeling Platform")

tabs = st.tabs(["Basic PAMA", "PAMA with Temperature", "PAMA with Degradation"])
multi_plot = st.sidebar.checkbox("📊 Show All Curves on Same Graph", value=False)

# ---- Tab 1: Basic PAMA ----
with tabs[0]:
    st.subheader("Basic PAMA")
    C1 = st.number_input("Concentration (g/L)", key="C1", value=2.0, min_value=0.1, max_value=50.0)
    MW1 = st.number_input("Molecular Weight (MDa)", key="MW1", value=8.0, min_value=0.1, max_value=100.0)
    eta7_1 = st.number_input("Experimental η₇ (cP)", key="eta7_1", value=15.0, min_value=0.1)
    run1 = st.button("Generate", key="btn1")

    if run1:
        df1, fig1 = model_basic_pama(C1, MW1, eta7_1)
        if multi_plot:
            st.session_state.fig_base = fig1
            st.session_state.df_base = df1
        else:
            st.plotly_chart(fig1, use_container_width=True)
            st.download_button("Download CSV", df1.to_csv(index=False), "basic_pama.csv", "text/csv")


# ---- Tab 2: PAMA with Temperature ----
with tabs[1]:
    st.subheader("PAMA with Temperature")
    C2 = st.number_input("Concentration (g/L)", key="C2", value=2.0)
    MW2 = st.number_input("Molecular Weight (MDa)", key="MW2", value=8.0)
    eta7_2 = st.number_input("Experimental η₇ (cP)", key="eta7_2", value=15.0)
    T2 = st.number_input("Temperature (°C)", key="T2", value=25)
    run2 = st.button("Generate", key="btn2")

    if run2:
        df2, fig2 = model_pama_temperature(C2, MW2, eta7_2, T2)
        if multi_plot:
            st.session_state.fig_temp = fig2
            st.session_state.df_temp = df2
        else:
            st.plotly_chart(fig2, use_container_width=True)
            st.download_button("Download CSV", df2.to_csv(index=False), "pama_temperature.csv", "text/csv")


# ---- Tab 3: PAMA with Degradation ----
with tabs[2]:
    st.subheader("PAMA with Degradation")
    C3 = st.number_input("Concentration (g/L)", key="C3", value=2.0)
    MW3 = st.number_input("Molecular Weight (MDa)", key="MW3", value=8.0)
    eta7_3 = st.number_input("Undegraded η₇ (cP)", key="eta7_3", value=15.0)
    eta7_3D = st.number_input("Degraded η₇ (cP)", key="eta7_3D", value=10.0)
    run3 = st.button("Generate", key="btn3")

    if run3:
        df3, fig3 = model_pama_degradation(C3, MW3, eta7_3, eta7_3D)
        if multi_plot:
            st.session_state.fig_deg = fig3
            st.session_state.df_deg = df3
        else:
            st.plotly_chart(fig3, use_container_width=True)
            st.download_button("Download CSV", df3.to_csv(index=False), "pama_degradation.csv", "text/csv")


# ---- Overlay All Curves ----
if multi_plot and any(k in st.session_state for k in ["fig_base", "fig_temp", "fig_deg"]):
    st.subheader("📊 Combined Comparison")
    combined_fig = go.Figure()

    if "fig_base" in st.session_state:
        df = st.session_state.df_base
        combined_fig.add_trace(go.Scatter(x=df["shear"], y=df["eta"], name="Basic PAMA", mode='lines'))

    if "fig_temp" in st.session_state:
        df = st.session_state.df_temp
        combined_fig.add_trace(go.Scatter(x=df["shear"], y=df["eta_temp"], name="PAMA Temp", mode='lines'))

    if "fig_deg" in st.session_state:
        df = st.session_state.df_deg
        combined_fig.add_trace(go.Scatter(x=df["shear"], y=df["eta"], name="Undegraded", mode='lines'))
        combined_fig.add_trace(go.Scatter(x=df["shear"], y=df["eta_degraded"], name="Degraded", mode='lines'))

    combined_fig.update_layout(
        xaxis=dict(type="log", title="Shear Rate (s⁻¹)"),
        yaxis=dict(type="log", title="Viscosity (cP)"),
        title="Overlay of All Selected PAMA Models",
        template="plotly_white"
    )
    st.plotly_chart(combined_fig, use_container_width=True)
