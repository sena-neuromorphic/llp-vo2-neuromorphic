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
