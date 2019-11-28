# Best Orthogonalized Subset Selection (BOSS)
Best Orthogonalized Subset Selection (BOSS) is a least-squares (LS) based subset selection method, that performs best subset selection upon an orthogonalized basis of ordered predictors, with the computational effort of a single ordinary LS fit.

This repository contains an R package **BOSSreg** that provides a highly optimized implementation of BOSS and an heuristic degrees of freedom, which can be further plugged into an information criterion such as AICc in order to selection the subset from candidates. Various choices of information criteria are provided including AICc, Cp, GCV, AIC and BIC. The R package also implements forward stepwise regression (FS) with no additional computational cost, where the subset of FS is selected via cross-validation (CV). CV is also an option for BOSS.

It also contains all the code to reproduce the results in the paper,
[Tian, S., Hurvich, C. and Simonoff, J. (2019): "On the Use of Information Criteria for Subset Selection in Least Squares Regression"](https://arxiv.org/abs/1911.10191).

### Install the R package
To install the latest version of the package, run the following:
```
library(devtools)
install_github(repo="sentian/BOSSreg", subdir="r-package")
```
Alternatively, a stable version can be installed from CRAN
```
install.packages("BOSSreg", repos = "http://cran.rstudio.com")
```

### Use the R package
For a simple guide of the functionalities of the package, please refer to the [BOSSreg's Vignette](https://github.com/sentian/BOSSreg/blob/master/r-package/vignettes/BOSSreg.pdf).

For a complete documentation of the package, please refer to the [BOSSreg's Documentation](https://github.com/sentian/BOSSreg/blob/master/BOSSreg_reference.pdf).

### Reproduce the results in the paper
The structure of the '*code*' directory is shown below.

'*plots.R*' and '*tables.R*' generate the figures and tables that can be found in the *paper/figures* and *paper/tables* directories, respectively. Note that '*plots.R*' is self-contained and does not rely on the simulation or real data results. '*tables.R*' requires simulation and real data results as inputs, which can be found in '*code/run_model/simulation/results*' and '*code/run_model/realdata/results*' directories, respectively.

'*.sh*' are bash files that submit the corresponding R code to a Linux server for running in parallels.

```bash
code
├── plots_tables              ## generate figures and tables in the paper
│   ├── plots.R               # code that reproduces all the figures
│   └── tables.R              # code that reproduces all the tables
├── run_model                 ## fit the models
│   ├── realdata
│   │   ├── results           # realdata results
│   │   │   └── *.rds
│   │   ├── forestfires.csv   # forest fire dataset
│   │   ├── run.R             # code that fits and evaluates methods
│   │   ├── run_forestfire.R  # code that fits MIO BS on forest fire dataset
│   │   └── run_forestfire.sh
│   └── simulation
│   │   ├── para_forhpc       # parameters for each configuration
│   │   │   └── *.txt
│   │   ├── results           # simulation results
│   │   │   ├── generalx      # results for a general X
│   │   │       └── *.rds
│   │   │   └── orthx         # results for an orthogonal X
│   │   │       └── *.rds
│   │   ├── run.R             # code to fit and evaluate all methods
│   │   ├── run_generalx.sh
│   │   └── run_orthx.sh
└── utils.R                   ## code that contains functions shared by other R codes
```
