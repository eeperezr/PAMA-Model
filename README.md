# 🧪 PAMA Rheology Modeling App

A Streamlit web application for simulating polymer viscosity behavior using the **PAMA method**. It provides three modeling approaches:
- 📘 Basic PAMA Model
- 🌡️ PAMA Model with Temperature Dependence
- ⚠️ PAMA Model with Degradation

---

## 🔗 Live App

> 🚀 [Try it Live on Streamlit Cloud](https://your-app-url.streamlit.app)  
*(Replace with your actual deployed link)*

---

## 📌 Features

- Interactive Streamlit UI with sliders and inputs
- Real-time dynamic plots using **Plotly**
- Export results as CSV
- Switch between 3 scientifically-backed models
- Log–log shear rate vs. viscosity plots
- Easy deployment-ready Python code

---

## 🧠 About the PAMA Models

### 1. Basic PAMA Model
- Predicts polymer viscosity from concentration and molecular weight
- Assumes standard conditions (25°C)
- Implements the Carreau-Yasuda equation

### 2. PAMA Model with Temperature
- Adjusts viscosity predictions based on user-specified temperature (in °C)
- Calculates temperature-adjusted η₀ and λ (relaxation time)
- Allows simulation of polymer performance in varying thermal environments

### 3. PAMA Model with Degradation
- Models viscosity reduction due to polymer degradation
- Uses two experimental η@7.3 values (undegraded and degraded)
- Computes degradation factor (FD) and predicts degraded rheological profile

---

## 📸 Screenshots

> 📈 Example plots will appear here.

---

## 🔧 How to Run Locally

1. **Clone this repository**
   ```bash
   git clone https://github.com/your-username/pama-streamlit.git
   cd pama-streamlit
