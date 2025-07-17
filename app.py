import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import json
import io

# --- PAGE CONFIG ---
st.set_page_config(
    page_title="PAMA Rheology Models",
    page_icon="🧪",
    layout="wide",
)

# --- CSS STYLING (from styles.css, inlined) ---
st.markdown(
    """
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Arial&display=swap');
    body {
        font-family: 'Arial', sans-serif;
        margin: 0;
        padding: 0;
        background-color: #f3f4f6;
        color: #1f2937;
    }
    header {
        background-color: #1e90ff;
        color: white;
        text-align: center;
        padding: 1rem;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    #sidebar {
        width: 16rem;
        position: fixed;
        top: 4rem;
        bottom: 0;
        background-color: #e5e7eb;
        padding: 1rem;
        border-right: 1px solid #d1d5db;
        overflow-y: auto;
        transition: all 0.3s ease;
    }
    #content {
        margin-left: 16rem;
        padding: 1.5rem;
    }
    .model-section {
        background-color: white;
        padding: 1rem;
        border-radius: 0.5rem;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        margin-bottom: 1.5rem;
    }
    select, input {
        width: 100%;
        padding: 0.5rem;
        border: 1px solid #d1d5db;
        border-radius: 0.375rem;
        margin-bottom: 0.5rem;
        outline: none;
        transition: border-color 0.2s ease;
    }
    select:focus, input:focus {
        border-color: #3b82f6;
        ring: 2px solid rgba(59,130,246,0.5);
    }
    button {
        background-color: #1e90ff;
        color: white;
        padding: 0.5rem 1rem;
        border-radius: 0.375rem;
        border: none;
        cursor: pointer;
        transition: background-color 0.2s ease;
    }
    button:hover {
        background-color: #104e8b;
    }
    #chart {
        width: 100%;
        height: 24rem;
        background-color: white;
        border-radius: 0.5rem;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
    }
    table {
        width: 100%;
        background-color: white;
        border-radius: 0.5rem;
        box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        margin-top: 1rem;
    }
    th {
        background-color: #e5e7eb;
        color: #374151;
        padding: 0.5rem;
    }
    td {
        padding: 0.5rem;
        border-top: 1px solid #d1d5db;
    }
    #download-link {
        color: #1e90ff;
        margin-top: 0.5rem;
        display: inline-block;
    }
    #download-link:hover {
        text-decoration: underline;
    }
    h1, h2, h3 {
        color: #1f2937;
    }
    ul {
        list-style-type: disc;
        padding-left: 1.5rem;
    }
    a {
        color: #1e90ff;
    }
    a:hover {
        color: #104e8b;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# --- HTML STRUCTURE (from index.html, adapted with download button) ---
html_content = """
<header>
    <h1>PAMA Rheology Modeling App 🧪</h1>
    <p>Tool for analyzing polymer rheology models with interactive visualizations.</p>
</header>
<div id="sidebar">
    <h3>Table of Contents</h3>
    <ul>
        <li><a href="#">PAMA Documentation</a></li>
        <ul>
            <li><a href="#">Contents</a></li>
            <li><a href="#">About PAMA</a></li>
            <li><a href="#">Installation</a></li>
            <li><a href="#">User Manual</a></li>
            <li><a href="#">PAMA for Developers</a></li>
            <li><a href="#">PAMA Contributors</a></li>
            <li><a href="#">Version History</a></li>
        </ul>
    </ul>
    <p>More Info:</p>
    <ul>
        <li><a href="https://example.com">Source Code</a></li>
    </ul>
    <p>Authors:</p>
    <ul>
        <li>Eduar Perez (University of Buenos Aires)</li>
    </ul>
</div>
<div id="content">
    <h2>PAMA Rheology Models</h2>
    <p>PAMA (Polymer Analysis and Modeling Application) is a tool for analyzing polymer rheology models with interactive visualizations.</p>
    <div class="model-section">
        <label for="model-select">Select Model:</label>
        <select id="model-select">
            <option value="basic">Basic PAMA</option>
            <option value="temperature">PAMA with Temperature</option>
            <option value="degradation">PAMA with Degradation</option>
        </select>
        <div id="params-basic" class="params">
            <label>Concentration (g/L):</label><span id="conc-basic"></span><br>
            <label>Molecular Weight (MDa):</label><span id="mw-basic"></span><br>
            <label>η@7.3 experimental (cP):</label><span id="eta7-basic"></span><br>
            <button id="run-basic">Run Basic Model</button>
        </div>
        <div id="params-temperature" class="params" style="display:none;">
            <label>Concentration (g/L):</label><span id="conc-temp"></span><br>
            <label>Molecular Weight (MDa):</label><span id="mw-temp"></span><br>
            <label>η@7.3 at 25°C (cP):</label><span id="eta7-temp"></span><br>
            <label>Target Temperature (°C):</label><span id="temp-target"></span><br>
            <button id="run-temp">Run Temperature Model</button>
        </div>
        <div id="params-degradation" class="params" style="display:none;">
            <label>Concentration (g/L):</label><span id="conc-deg"></span><br>
            <label>Molecular Weight (MDa):</label><span id="mw-deg"></span><br>
            <label>η@7.3 (Polymer UD) (cP):</label><span id="eta7-ud"></span><br>
            <label>η@7.3 (Polymer D) (cP):</label><span id="eta7-d"></span><br>
            <button id="run-deg">Run Degradation Model</button>
        </div>
        <div id="chart"></div>
        <table id="data-table" style="display:none;">
            <thead><tr><th>Shear</th><th>Viscosity</th></tr></thead>
            <tbody></tbody>
        </table>
        <div id="download-container"></div>
    </div>
</div>
<script>
    const modelSelect = document.getElementById('model-select');
    const params = document.querySelectorAll('.params');
    const chart = document.getElementById('chart');
    const dataTable = document.getElementById('data-table');
    const downloadContainer = document.getElementById('download-container');

    modelSelect.addEventListener('change', (e) => {
        params.forEach(p => p.style.display = 'none');
        document.getElementById(`params-${e.target.value}`).style.display = 'block';
    });

    function showLoading() {
        chart.innerHTML = '<div class="flex items-center justify-center h-full"><span class="text-gray-500">Loading...</span></div>';
        dataTable.style.display = 'none';
        downloadContainer.innerHTML = '';
    }

    function updateParams() {
        const params = {
            'basic': { conc: 'conc-basic', mw: 'mw-basic', eta7: 'eta7-basic' },
            'temperature': { conc: 'conc-temp', mw: 'mw-temp', eta7: 'eta7-temp', temp: 'temp-target' },
            'degradation': { conc: 'conc-deg', mw: 'mw-deg', eta7ud: 'eta7-ud', eta7d: 'eta7-d' }
        };
        const model = modelSelect.value;
        const p = params[model];
        if (p.conc) document.getElementById(p.conc).textContent = ` ${window.st_params.conc || 2.0} (Slider: ${window.st_params.conc || 2.0})`;
        if (p.mw) document.getElementById(p.mw).textContent = ` ${window.st_params.mw || 8.0} (Slider: ${window.st_params.mw || 8.0})`;
        if (p.eta7) document.getElementById(p.eta7).textContent = ` ${window.st_params.eta7 || 20.0} (Input: ${window.st_params.eta7 || 20.0})`;
        if (p.temp) document.getElementById(p.temp).textContent = ` ${window.st_params.temp || 35} (Slider: ${window.st_params.temp || 35})`;
        if (p.eta7d) document.getElementById(p.eta7d).textContent = ` ${window.st_params.eta7d || 7.354} (Input: ${window.st_params.eta7d || 7.354})`;
    }

    document.getElementById('run-basic').addEventListener('click', () => {
        showLoading();
        const data = window.model_basic_pama(
            window.st_params.conc || 2.0,
            window.st_params.mw || 8.0,
            window.st_params.eta7 || 20.0
        );
        updateChart(data, 'basic');
        window.setData(data, 'basic');
    });

    document.getElementById('run-temp').addEventListener('click', () => {
        showLoading();
        const data = window.model_pama_temperature(
            window.st_params.conc || 2.0,
            window.st_params.mw || 8.0,
            window.st_params.eta7 || 15.653,
            window.st_params.temp || 35
        );
        updateChart(data, 'temperature');
        window.setData(data, 'temperature');
    });

    document.getElementById('run-deg').addEventListener('click', () => {
        showLoading();
        const data = window.model_pama_degradation(
            window.st_params.conc || 2.0,
            window.st_params.mw || 8.0,
            window.st_params.eta7 || 15.653,
            window.st_params.eta7d || 7.354
        );
        updateChart(data, 'degradation');
        window.setData(data, 'degradation');
    });

    function updateChart(data, model) {
        let traces = [];
        let columns = [];
        if (model === 'basic') {
            traces.push({ x: data.shear, y: data['Viscosity (cP)'], name: 'Viscosity', mode: 'lines+markers' });
            columns = ['shear', 'Viscosity (cP)'];
        } else if (model === 'temperature') {
            traces.push({ x: data.shear, y: data['25°C Reference'], name: '25°C', mode: 'lines+markers' });
            traces.push({ x: data.shear, y: data[Object.keys(data)[2]], name: Object.keys(data)[2], mode: 'lines+markers' });
            columns = ['shear', '25°C Reference', Object.keys(data)[2]];
        } else if (model === 'degradation') {
            traces.push({ x: data.shear, y: data['Polymer UD'], name: 'Polymer UD', mode: 'lines+markers' });
            traces.push({ x: data.shear, y: data['Polymer Degraded'], name: 'Polymer Degraded', mode: 'lines+markers' });
            columns = ['shear', 'Polymer UD', 'Polymer Degraded'];
        }

        Plotly.newPlot(chart, traces, {
            title: `${model.charAt(0).toUpperCase() + model.slice(1)} PAMA Viscosity vs Shear Rate`,
            xaxis: { title: 'Shear rate (s⁻¹)', type: 'log', gridcolor: 'lightgray' },
            yaxis: { title: 'Viscosity (cP)', type: 'log', gridcolor: 'lightgray' },
            template: 'plotly_white',
            legend: { orientation: 'h', yanchor: 'bottom', y: 1.02, xanchor: 'center', x: 0.5 },
            margin: { t: 50, b: 40, l: 40, r: 40 },
            hovermode: 'x unified'
        });

        const tbody = dataTable.querySelector('tbody');
        tbody.innerHTML = '';
        data.shear.forEach((s, i) => {
            const row = tbody.insertRow();
            columns.forEach(col => row.insertCell().textContent = data[col][i].toFixed(3));
        });
        dataTable.style.display = 'table';
    }

    // Expose setData function to Python
    window.setData = (data, model) => {
        window.parent.postMessage({ type: 'setData', data: data, model: model }, '*');
    };
</script>
"""

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

# --- EXPOSE MODELS AND PARAMETERS TO JAVASCRIPT ---
st.components.v1.html(f"""
<script>
    window.model_basic_pama = {json.dumps(model_basic_pama.__code__.co_consts[0])};
    window.model_pama_temperature = {json.dumps(model_pama_temperature.__code__.co_consts[0])};
    window.model_pama_degradation = {json.dumps(model_pama_degradation.__code__.co_consts[0])};
    window.st_params = {json.dumps(st.session_state)};
    // Receive data from JS and set in session_state
    window.addEventListener('message', (event) => {
        if (event.data.type === 'setData') {
            st.session_state.data = event.data.data;
            st.session_state.model = event.data.model;
        }
    });
</script>
""" + html_content, height=800)

# --- DOWNLOAD BUTTON ---
if 'data' in st.session_state and 'model' in st.session_state:
    data = st.session_state.data
    model = st.session_state.model
    columns = []
    if model == 'basic':
        columns = ['shear', 'Viscosity (cP)']
    elif model == 'temperature':
        columns = ['shear', '25°C Reference', list(data.keys())[2]]
    elif model == 'degradation':
        columns = ['shear', 'Polymer UD', 'Polymer Degraded']
    
    df = pd.DataFrame({col: data[col] for col in columns})
    csv = df.to_csv(index=False)
    st.download_button(
        label="Download CSV",
        data=csv,
        file_name=f"{model}_pama.csv",
        mime="text/csv",
        key="download-button"
    )

# --- INTERACTIVE WIDGETS ---
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

st.session_state.conc = st.sidebar.slider("Concentration (g/L)", 0.1, 20.0, 2.0, 0.1)
st.session_state.mw = st.sidebar.slider("Molecular Weight (MDa)", 0.1, 50.0, 8.0, 0.1)
st.session_state.eta7 = st.sidebar.number_input("η@7.3 experimental (cP)", min_value=0.01, value=20.0, format="%.3f")
st.session_state.temp = st.sidebar.slider("Target Temperature (°C)", 0, 100, 35, 1)
st.session_state.eta7d = st.sidebar.number_input("η@7.3 experimental (Polymer D) (cP)", min_value=0.01, value=7.354, format="%.3f")
