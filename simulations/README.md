# Simulation scripts

R scripts reproducing the simulations and real-data analyses of the papers
built on `picreg`. Each folder holds the scripts and, in `results/`, the
outputs they produced (CSV, RDS, PDF).

| Folder | Paper |
| --- | --- |
| `jasa_paper/` | Sardy, van Cutsem & van de Geer, *The Pivotal Information Criterion* ([arXiv:2603.04172](https://arxiv.org/abs/2603.04172)) |
| `cox_paper/` | van Cutsem & Sardy, *Survival of the fittest Cox model: Pivotal variable selection for time-to-event data* ([arXiv:2510.19374](https://arxiv.org/abs/2510.19374)) |
| `lassoVSgumbel/` | van de Geer, Sardy & van Cutsem, *A pivotal transform for the high-dimensional location-scale model* ([arXiv:2512.18705](https://arxiv.org/abs/2512.18705)) |

## `jasa_paper/`

| Script | Paper section |
| --- | --- |
| `simu1_phase_transition.R`, `simu1_plot.R` | 4.2, phase transition (Figures 3 and 4) |
| `simu2_weak_signals.R` | 4.3, weak signals and correlated designs (Table 2) |
| `simu3_real_data.R` | 4.4, real data (Table 4) |
| `simu4_calibration.R`, `simu4_plot.R` | 3.1, sensitivity of the calibration (Figure 1) |
| `simu5_illustration.R`, `simu5_plot.R` | 4.1, pivotal versus non-pivotal calibration (Figure 2) |

`common.R` holds the helpers shared by the simulation scripts.

## Data

The UCI datasets used by `jasa_paper/simu3_real_data.R` (`communities.data`,
`ionosphere.data`, `wdbc.data`, `sonar.csv`) are included in
`jasa_paper/data/`; the first three are downloaded again from the UCI
repository if missing. The Diabetes data come from the R package `lars`, and
the PBC, GBSG, lung and NKI70 data of `cox_paper/real_data_simu.R` from the R
packages `survival` and `penalized`.

The following datasets are not redistributed here. Place them at the paths
below to rerun the corresponding analyses:

| File | Source |
| --- | --- |
| `jasa_paper/data/riboflavin.csv` | `riboflavin` data of the R package `hdi`, exported with the response in a column `y` |
| `jasa_paper/data/ProstateCancer.csv` | Prostate cancer data of Stamey et al. (1989), as distributed on Kaggle (with columns `lpsa` and `train`) |
| `cox_paper/ChopData.csv` | CHOP data, see the references of the Cox paper |
| `cox_paper/mclData.csv` | MCL data, see the references of the Cox paper |
| `cox_paper/LymphomaData.rda` | Lymphoma data, see the references of the Cox paper |
