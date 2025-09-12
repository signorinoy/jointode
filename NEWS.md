# JointODE Development Roadmap

## Version 0.1.0 (Current Release)

* Joint modeling of longitudinal biomarkers and survival outcomes using ODEs
* Second-order differential equation formulation for biomarker dynamics
* EM algorithm with parallel processing support

## Version 0.2.0 (Planned - Q2 2025)

### Subgroup Heterogeneity

* `JointODE_group()`: Latent subgroup modeling
* Group-specific ODE parameters (κ_g, γ_g)
* Model selection via ICL and entropy criteria
* K-means initialization with stability analysis

## Version 0.3.0 (Planned - Q4 2025)

### Multiple Biomarkers

* `JointODE_multi()`: Multi-marker joint modeling
* `select_biomarkers()`: SIP-based variable selection
* Adaptive LASSO and group penalties
* Rcpp integration for performance

## Version 0.4.0 (Planned - Q4 2026)

### Machine Learning

* `JointODE_nn()`: Neural ODE backend
* Python bridge via reticulate
* Cloud deployment tools
