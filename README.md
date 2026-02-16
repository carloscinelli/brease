# brease: Causally Sound Priors for Binary Experiments

<!-- badges: start -->
[![R-CMD-check](https://github.com/carloscinelli/brease/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/carloscinelli/brease/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `brease` R package implements the **BREASE** (**B**aseline **R**isk, **E**fficacy, and **A**dverse **S**ide **E**ffects) framework for Bayesian analysis of randomized controlled trials with binary treatment and binary outcome, as described in:

> Irons, N.J. and Cinelli, C. (2025). [Causally Sound Priors for Binary Experiments](https://doi.org/10.1214/25-BA1506). *Bayesian Analysis*.

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("carloscinelli/brease")
```

## Quick Example

```r
library(brease)

# Aspirin and Myocardial Infarction (Physicians' Health Study)
result <- brease(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)
result

# COVID-19 Vaccine Trial (Pfizer)
result <- brease(y0 = 169, y1 = 9, N0 = 20172, N1 = 19965)
result

# Sensitivity analysis
sens <- bf_seq(y0 = 26, y1 = 10, N0 = 11034, N1 = 11037)
contour_bf(sens)
```

## Overview

The BREASE framework parameterizes the likelihood of a binary experiment in terms of three clinically meaningful causal quantities:

- **Baseline risk** (theta_0): probability of the event without treatment
- **Efficacy** (eta_e): fraction of would-be cases prevented by treatment
- **Side-effect risk** (eta_s): fraction of would-be non-cases harmed by treatment

The package provides:

- **Exact posterior sampling** and analytic marginal likelihoods/Bayes factors
- **Optional MCMC backends** via Stan or JAGS
- **Sensitivity analysis** tools (contour plots, side-effect sensitivity)
- **Comparison methods**: Independent Beta and Logit-Transform priors
- **Bounds** on probabilities of causation

## Documentation

For a detailed introduction, see the [package vignette](https://carloscinelli.com/brease/articles/brease.html).

## Citation

If you use this package, please cite:

```
Irons, N.J. and Cinelli, C. (2025). Causally Sound Priors for Binary
Experiments. Bayesian Analysis. doi:10.1214/25-BA1506
```
