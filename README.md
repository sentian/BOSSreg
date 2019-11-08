# Best Orthogonalized Subset Selection (BOSS)
This repository contains the R package **BOSSreg** that implements the Best Orthogonalized Subset Selection (BOSS) and Forward Stepwise regression (FS), with feasible selection rules that include information criteria and cross-validation. Various choices of information criteria are provided that include AICc, Cp, GCV, AIC and BIC.

It also contains all the code to reproduce the results in the paper,
[Tian, Hurvich and Simonoff (2019): "On the use of information criterion for least squares based subset selection methods"](https://github.com/sentian/BOSSreg/blob/master/paper/Tian2019.pdf)

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
For a simple guide of the functionalities of the package, please refer to the [BOSSreg's Vigette](https://github.com/sentian/BOSSreg/blob/master/r-package/vignettes/BOSSreg.pdf).

For a complete documentation of the package, please refer to the [BOSSreg's Documentation](https://github.com/sentian/BOSSreg/blob/master/BOSSreg_0.1.0.pdf).

### Reproduce the results in the paper
The structure of the '*code*' directory is as follows. '*plots.R*' and '*tables.R*' are the codes to reproduce the figures and tables in the paper, respectively.
```bash
code
├── plots_tables         ## generate figures and tables in the paper
│   ├── plots            # the figures
│   │   ├── *
│   ├── tables           # the tables
│   │   ├── *
│   ├── plots.R          # code that reproduces all the figures
│   └── tables.R         # code that reproduces all the tables
├── run_model            ## fit the models and evaluate their performances
│   ├── para_forhpc      # parameters for each configuration
│   │   ├── *
│   ├── results          # simulation results
│   │   ├── generalx     # results for a general X
│   │   │   ├── *
│   │   └── orthx        # results for an orthogonal X
│   │       ├── *
│   ├── run.R            # code to fit and evaluate all methods
│   ├── run_generalx.sh  # bash file that submits parallel jobs to a Linux server
│   └── run_orthx.sh     # bash file that submits parallel jobs to a Linux server
└── utils.R              ## code that contains functions shared by all the other R codes
```
