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


Background & References

This method and model are based on polymer rheology theory, Carreau-Yassuda viscosity models, and the specific needs of measuring polymer solutions with limited instrumentation. For further reading, please consult:

    E. Pérez, D. Alviso, E. Manrique, G. Artana, Laboratorio de Fluidodinámica, Facultad de Ingeniería, Universidad de Buenos Aires.

    Pérez, E., Alviso, D., Manrique, E., & Artana, G. (2022). Estimation of the rheological curve of HPAM solutions from measurements using        the Brookfield viscometer. Journal of Petroleum Science and Engineering, 216, 110793.
    
    Pérez, E., Alviso, D., Carmona, M., Romero, J., Manrique, E., & Artana, G. (2023, October). Estimation of the concentration of HPAM            polymer required to achieve a desired viscosity in EOR projects. In IOR+ 2023 (Vol. 2023, No. 1, pp. 1-10). European Association of Geoscientists & Engineers.

    Pérez, E., Alviso, D., Carmona, M., Manrique, E., & Artana, G. (2024, October). Estimation of the Molecular Weight of HPAM Solutions from      Brookfield Viscometer Measurements at 298K. In Third EAGE Workshop on EOR (Vol. 2024, No. 1, pp. 1-6). European Association of                 Geoscientists & Engineers.

    Pérez, E., Alviso, D., Carmona, M., Manrique, E., & Artana, G. (2024). A simple model of the rheological curve of HPAM solutions at            different temperatures. Scientific Reports, 14(1), 31601.

    

License

MIT License — feel free to use, modify, and distribute!
Contact

For questions or collaboration, contact:

    E. Pérez: eperez@fi.uba.ar
