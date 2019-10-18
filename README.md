# Best orthogonalized subset selection (BOSS)
This repository contains the **R package** that implements the best orthogonalized subset selection (BOSS), a least squares based subset selection algorithm. The details of the methodology can be found in our paper: Tian, Hurvich and Simonoff (2019): "On the use of information criterion for least squares based subset selection methods", which can be found in the **paper** directory. We also provide a [Vigettes](https://github.com/sentian/boss/blob/master/r-package/vignettes/boss.html) that is inside the **rpackage** directory. Finally, the code for all the simulation results in our paper can be found in the **simulation** directory. 

## Install the R package
To install the **boss** package, run the following:
```
library(devtools)
install_github(repo="sentian/boss", subdir="r-package")
```
