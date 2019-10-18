# Best orthogonalized subset selection (BOSS)
This repository contains the **R package** that implements the best orthogonalized subset selection (BOSS), a least squares based subset selection algorithm. The details of the methodology can be found in our paper: [Tian, Hurvich and Simonoff (2019): "On the use of information criterion for least squares based subset selection methods"](https://github.com/sentian/boss/blob/master/r-package/vignettes/Tian2019.pdf). We also provide a [Vigette](https://github.com/sentian/boss/blob/master/r-package/vignettes/boss.pdf) for the package. Finally, the code for all the simulation results in our paper can be found in the **simulation** directory. 

## Install the R package
To install the **boss** package, run the following:
```
library(devtools)
install_github(repo="sentian/boss", subdir="r-package")
```
