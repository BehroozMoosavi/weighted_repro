# weighted_repro

Code and archived results for the paper

> **Finite-Sample Valid and Asymptotically Efficient Simulation-Based Inference**
> Behrooz Moosavi Ramezanzadeh and Jordan Awan, University of Pittsburgh.

The paper studies how to choose the test statistic in **SimBaRepro**, a simulation-based method that gives confidence sets and tests valid in finite samples whenever the data can be simulated from a parametric model, as is the case for differentially private releases. It shows that Mahalanobis depth is asymptotically optimal for joint inference on the whole parameter, and it proposes the **Penalized Wald** statistic for inference on a subset of the parameters in the presence of nuisance parameters. This repository contains the simulation code behind every numerical result in the paper and the archived outputs used to make its tables and figures.

The release tagged **`aos-submission-v1`** corresponds to the submitted manuscript.

## Repository layout

| Folder | Contents |
|---|---|
| `codes/` | Simulation scripts, one folder per experiment. |
| `results/` | Archived outputs: data files, LaTeX table fragments read by the paper, and figures. `results/MANIFEST.sha256` lists a SHA-256 checksum for every archived file. |
| `tabulation/` | Scripts that turn the archived outputs into the LaTeX tables. |

## Scripts and the results they produce

| Script | Location | Result in the paper |
|---|---|---|
| `exp1_coverage.R` | `codes/location_scale_normal/exp1_coverage.R` | Table 1: Experiment 1, private location–scale model |
| `exp1_lambda_grid_optimal.R` | `codes/sensitivity_analysis/exp1_normal_sensitivity/exp1_lambda_grid_optimal.R` | Sensitivity to the penalty weight (supplementary simulations) |
| `exp1_R_point_cloud_sensitivity.R` | `codes/sensitivity_analysis/exp1_normal_sensitivity/exp1_R_point_cloud_sensitivity.R` | Sensitivity to the inference-cloud size |
| `exp1_R_aux_sensitivity.R` | `codes/sensitivity_analysis/exp1_normal_sensitivity/exp1_R_aux_sensitivity.R` | Sensitivity to the auxiliary-cloud size |
| `exp2_linearreg_eff.R` | `codes/linear_regression/exp2_linearreg_eff.R` | Figure 2: Experiment 2, test of a regression slope (simulation) |
| `heat map.R` | `codes/linear_regression/heat map.R` | Figure 2 (plotting only) |
| `two_target_dim_mh_pw_pbadi.R` | `codes/linear_regression/two_target_dim_mh_pw_pbadi.R` | Table 2: joint target in private linear regression |
| `clamp_sensitivity.R` | `codes/sensitivity_analysis/exp2_linear_regression_sensitivity/clamp_sensitivity.R` | Sensitivity to the clamping range (simulation) |
| `clamp_plots_only.R` | `codes/sensitivity_analysis/exp2_linear_regression_sensitivity/clamp_plots_only.R` | Sensitivity to the clamping range (plotting only) |
| `exp2_sigma_penalized_efficient.R` | `codes/linear_regression/exp2_sigma_penalized_efficient.R` | Scale target in private linear regression |
| `eff.R` | `codes/objective_perturbation/eff.R` | Figure 3: Experiment 3, logistic regression by objective perturbation (simulation) |
| `heatmap.R` | `codes/objective_perturbation/heatmap.R` | Figure 3 (plotting only) |
| `failure_v1.R` | `codes/pb_ADI_failure/failure_v1.R` | Boundary study, hard clamp |
| `failure_v2.R` | `codes/pb_ADI_failure/failure_v2.R` | Boundary study, hard and soft clamp at the boundary |

Plotting and tabulation scripts read the saved results and do not rerun any simulation.

## Requirements

All experiments are written in R. The R version and the versions of every package loaded by the scripts are recorded in `SESSION_INFO.md`. Parallel runs use PSOCK clusters through `foreach` and `doSNOW`, with BLAS and OpenMP threads limited to one inside each worker; the number of workers changes the run time but not the results.

## Reproducing the results

1. Clone the repository and check out the release used for the paper:
   ```sh
   git clone https://github.com/BehroozMoosavi/weighted_repro.git
   cd weighted_repro
   git checkout aos-submission-v1
   ```
2. Run a simulation script from its own folder, for example
   ```sh
   cd codes/location_scale_normal
   Rscript exp1_coverage.R
   ```
3. Run the matching plotting or tabulation script to rebuild the figure or table.

Every replication is seeded explicitly, and within a replication the observed release, the inference cloud and the auxiliary cloud are drawn as consecutive blocks of one random-number stream. Results are therefore reproducible regardless of how many workers are used. The seed of each study is given in the Implementation Details appendix of the paper.

The full Monte Carlo studies use 1000 replications per design point and take many hours on a multicore machine. For a quick check, lower the number of replications defined near the top of a script.

## Comparator results

The comparator values for the slope test (Figure 2) and the logistic-regression study (Figure 3) are the archived results of the public implementation of Wang, Chang and Awan, [`JordanAwan/debiased_dp_parametric_bootstrap`](https://github.com/JordanAwan/debiased_dp_parametric_bootstrap), at commit `232ec0b`: the files `Figure5/results/*.csv` and the values committed in `Figure6/plot_objpert.py`. The comparators in Experiment 1 and in the joint-target study are recomputed with the code in this repository.

## Numerical conventions

The confidence sets in the paper are defined through a supremum over the nuisance parameters, and the code approximates that supremum with a finite search. Finding an accepted value certifies acceptance, but failing to find one does not certify rejection. Following the implementation of Awan and Wang, the search in Experiment 1 and in the slope test of Experiment 2 includes the data-generating value among its starting points; the other studies use starting values computed from the data. The studies also differ in how they record a search that does not resolve. The Implementation Details appendix of the paper describes the search and the failure convention used in each study.

## Known limitations of the archive

- The archived Experiment 1 table does not retain the interval-widening and boundary-hit flags.
- The Experiment 3 plotting archive reports coverage over resolved replications and does not retain the complete counts of unresolved replications.
- The grid-window diagnostics of the joint-target study are not archived.
- The clamping-range study and the slope study share eight design points; their archived rejection rates differ at five of them.

The paper discloses each of these points.

## Table captions

The LaTeX fragments in `results/` carry the captions of the submitted manuscript, which were revised after the tables were generated. Rerunning the tabulation scripts reproduces the same numbers with the earlier caption text.

## License and citation

The code is released under the GPL-3.0 license (see `LICENSE`). If you use it, please cite the paper; `CITATION.cff` gives the citation metadata.
