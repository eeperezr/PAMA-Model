# PAMA-Model

PAMA Method for HPAM Rheology — Multi-Implementation

This repository contains implementations of the PAMA (Polymer Apparent Modified Apparatus) method to estimate the full rheological curve (viscosity vs shear rate) of HPAM polymer solutions based on limited experimental data from Brookfield viscometers.
Overview

The PAMA method is designed to derive a full viscosity vs. shear rate curve for HPAM polymers when only Brookfield-type viscometers are available, which measure viscosity at a single shear rate (7.3 s⁻¹). This method calculates rheological properties considering polymer concentration, molecular weight, and experimental viscosities, also accounting for polymer degradation.
Implementations
1. Original R Shiny App Version

    Technology: R with Shiny, Shinythemes, Plotly

    Features:

        User interface for input parameters and results visualization.

        Outputs Carreau-Yassuda viscosity model plots.

        Data download functionality.

    Use case: Interactive web app for easy parameter tweaking and instant visualization.

Run:

# Install dependencies
install.packages(c("shiny", "shinythemes", "plotly"))

# Run app
shiny::runApp("app_folder_path")

2. Modified R Shiny App with Degradation Parameter

    Technology: R with Shiny, Shinythemes, Plotly

    Improvements:

        Added handling for degraded polymer solutions.

        Inputs include degraded viscosity at 7.3 s⁻¹.

        Shows viscosity curves for both undegraded and degraded polymers.

    Use case: Web app that models polymer degradation effects on viscosity.

Run:

# Install dependencies
install.packages(c("shiny", "shinythemes", "plotly"))

# Run app
shiny::runApp("modified_app_folder_path")

3. Python Standalone Script (No App)

    Technology: Python 3.x, NumPy, Pandas, Plotly

    Features:

        Command-line script for PAMA method calculations.

        Inputs are defined directly in the script (polymer concentration, molecular weight, viscosities).

        Outputs data tables and saves to carreau_data.txt.

        Generates interactive Plotly plots with:

            Viscosity curve for undegraded polymer

            Viscosity curve for degraded polymer

        Fully reproducible outside R/Shiny environment.

    Use case: Lightweight, flexible script for batch processing or integration into larger Python workflows.

Run:

pip install numpy pandas plotly
python pama_method.py

Inputs Explained
Parameter	Description	Units	Default Example
concentration	Polymer concentration in solution	g/L	2
molecular_weight	Polymer molecular weight	MDa	8
eta7_exp	Experimental viscosity at 7.3 s⁻¹ (undegraded)	cP	15.653
eta7_exp_D	Experimental viscosity at 7.3 s⁻¹ (degraded)	cP	7.354


Output

    Data file: carreau_data.txt containing shear rates and viscosities.

    Plots: Interactive log-log plots showing the rheological curves of HPAM polymers, visualizing the effect of degradation.

Background & References

This method and model are based on polymer rheology theory, Carreau-Yassuda viscosity models, and the specific needs of measuring polymer solutions with limited instrumentation. For further reading, please consult:

    E. Pérez, D. Alviso, E. Manrique, G. Artana, Laboratorio de Fluidodinámica, Facultad de Ingeniería, Universidad de Buenos Aires.

    Relevant rheology textbooks and Carreau-Yassuda modeling papers.

License

MIT License — feel free to use, modify, and distribute!
Contact

For questions or collaboration, contact:

    E. Pérez: eperez@fi.uba.ar

    D. Alviso: dalviso@fi.uba.ar
