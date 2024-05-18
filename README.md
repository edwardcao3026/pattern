# Universal deterministic patterns in stochastic count data

This repository contains the Julia code and generated plots for the paper [1], where the corresponding data are deposited at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11213863.svg)](https://doi.org/10.5281/zenodo.11213863)

**Requirements:**

- Julia v1.9.1
- DelimitedFiles v0.10.12
- Plots v1.39.0
- Dates
- StatsBase v0.33.21
- HDF5 v0.17.1
- SparseArrays v1.10.0
- XLSX v0.10.1

All the codes have been tested on a MacBook Pro with Apple M3 Pro chip (11 cores) and 18 GB RAM.

All the curves are computed using the equation

![equation](https://latex.codecogs.com/svg.image?\langle&space;n\rangle=\frac{1}{2}\left(1-\text{FF}&plus;\sqrt{\frac{8k&plus;n_c(1-\text{FF})^2}{n_c}}\right))

of which the details are described in [1].

**File description:**

- `fig_plot.jl` is the main code for data processing and generating figures.
- `FigS1.svg` and `Fig4.svg` are the results generated by the code `fig_plot.jl` and presented in [1].

**Reference:**

- [1] Z. Cao, Y. Wang, R. Grima. *Universal deterministic patterns in stochastic count data*. 
