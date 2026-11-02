# MargRegVarPro

MATLAB implementations and numerical experiments for the paper:

> **Marginalization-Based Regularization for Separable Nonlinear Least Squares Problems with Applications to Semi-Blind Deblurring**  
> J. T. Dworaczyk and M. I. Español.

The routines reproduce the experiments, figures, and numerical tables presented in the paper.

---

## Project Description

MargRegVarPro implements regularized variable projection for semi-blind image deblurring, where the image is a linear variable and the blur parameters are nonlinear variables. A log-determinant penalty motivated by Gaussian marginalization regularizes the blur parameters independently of the image prior.

The implementation includes:

- Fourier-domain evaluation of the reduced objective and its derivatives for Gaussian, defocus, and motion blur.
- A quasi-Newton solver with diagonal damping and Armijo backtracking.
- Comparisons with variable projection without the log-determinant penalty and with Gaussian marginalization.
- Objective landscapes, auxiliary-scale sensitivity studies, and an illustrative experiment using conditional generalized cross-validation for parameter selection.

---

## Contents of the Repository

- **Solver and PSF routines:** the `solve_*.m` functions, `reconstruction_objective.m`, and `*_psf_symbol_derivs.m` implement the optimization and blur models.
- **Experiment drivers:** `run_all.m` runs the complete workflow; the other `run_*.m` functions run individual studies and parameter searches.
- **Figure and table exporters:** `plot_*.m` and `write_manuscript_tables.m` generate the paper's experimental figures and seven numerical tables.
- **Input images:** `cameraman.png` and `cameraman_reconstruction.png` supply the image data used in the experiments.
- **Verification:** `verify_code.m` checks derivatives, objective identities, and solver interfaces; `verify_tables.m` compares generated tables with the manuscript references in `reference_tables/`.

---

## Software Requirements

- MATLAB; tested with R2024b Update 1 on Windows.
- Image Processing Toolbox.

Optimization Toolbox is not required.

---

## Reproducibility

Download or clone the repository, open its folder in MATLAB, and run:

```matlab
run_all
```

This runs the numerical experiments and parameter searches, exports the figures and tables, and performs the verification checks. The random seeds and experiment settings are fixed in the code. Generated files are saved under `output/`; precomputed results are not needed.

After a complete run, regenerate the figures and tables from the saved results with:

```matlab
run_all('figures')
```

This mode reevaluates the objective grids but reuses the saved optimization results. Individual experiment functions can also be called by name.

The main experiment functions are:

| Experiment | Function |
|---|---|
| Objective landscapes | `run_objective_landscapes` |
| Image and blur reconstruction | `run_reconstruction_gaussian`, `run_reconstruction_defocus`, `run_reconstruction_motion` |
| Reference-image weight selection | `run_alpha_sweep_gaussian`, `run_alpha_sweep_defocus`, `run_alpha_sweep_motion` |
| Auxiliary-scale sensitivity | `run_auxiliary_scale_sensitivity` |
| Conditional GCV parameter selection | `run_gcv_single_sample` |

The reconstruction weight searches use the known reference image to compare candidate weights. The conditional GCV experiment selects parameters from the observation.

Solver stop codes are `0` for the iteration limit, `1` for the gradient tolerance, and `2` for a small step, including line-search stagnation. In the GCV results, `alternation_converged` records the relative-change tests for the alternating updates; the inner solver stop codes are stored in column 7 of the fit histories.

Runtime and floating-point stopping behavior may vary across MATLAB environments.

---

## Citation

When using MargRegVarPro in your research, please cite:

> J. T. Dworaczyk and M. I. Español. *Marginalization-Based Regularization for Separable Nonlinear Least Squares Problems with Applications to Semi-Blind Deblurring*. Manuscript, 2026.
