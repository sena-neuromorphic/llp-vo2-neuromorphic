# llp-vo2-neuromorphic

This repository contains the code for the paper:

*Explicit Electrothermal LLP VO₂ Model Reproducing Preisach-Like Hysteresis for Memristive and Neuromorphic Devices*  
(accepted for publication in *Scientific Reports*).

## Abstract

Vanadium dioxide (VO₂) devices exhibit a thermally driven metal–insulator transition with pronounced hysteresis, making them promising candidates for memristive and neuromorphic applications.

This repository implements a fully explicit electrothermal LLP–VO₂ model, designed to enable practical use of LLP hysteresis operators in VO₂-based systems. The model couples a thermal balance equation with a discrete LLP operator governing the metallic fraction, capturing key Preisach-like features such as nested hysteresis loops and return-point memory.

The formulation combines Arrhenius-type conduction with a metallic contribution and is evaluated using both explicit Euler and fourth-order Runge–Kutta (RK4) integration schemes to analyze numerical accuracy and stiffness effects.

The results demonstrate that LLP-based hysteresis can be simulated efficiently without additional differential states, providing a lightweight, reproducible, and implementation-ready framework compatible with SPICE-class solvers for future integration and benchmarking.

## Notes

The code was originally developed and tested in a Google Colab environment. Running it in Colab is recommended for faster and smoother execution.

## Main Scripts

The following scripts reproduce the results presented in the paper:

- `Fig2and3_main_article.py`  
  Reproduces Figures 2 and 3 from the main article, as well as Figures 1 and 2 from the Supplementary Information.
  
- `figure_4_triangular_excitation_main_article.py`  
  Reproduces Figure 4 from the main article (damped triangular excitation). Additional figures are also generated but are not included in the manuscript.

- `figure_5_forc_main_article.py`  
  Reproduces Figure 5 from the main article, including the FORC protocol and the corresponding hysteretic response.

- `Figure6_LLP_Experimental.py`
   Reproduces Figure 6 from the main article.
   If you are using Google Colab, make sure to upload the file `Experimental_Data_VO2.xlsx` to the working directory before running the scripts.
   The file `Experimental_Data_VO2.xlsx` contains the experimental data and is only required for the script `Figure6_LLP_Experimental.py`.
  The file `Chapter2_2003.pdf` contains a detailed description of the experimental platform, including the setup and data acquisition procedures used to obtain the experimental data.

  All provided scripts can be used to model VO$_2$ films using the explicit electrothermal LLP formulation. For questions, please contact the authors. If you use this code, please cite the associated publication.
  If you use the LLP model to study other metal–insulator transition materials, feel free to share your results with us — we would be happy to hear about your work.

  ## Requirements

The code was developed and tested in a Google Colab environment.

All required dependencies are listed in the `requirements.txt` file. These packages are available by default in Google Colab, enabling straightforward and reproducible execution.
