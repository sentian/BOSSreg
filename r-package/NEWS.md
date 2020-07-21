# BOSSreg v0.2.0.9000
In this development version (available on github but not on CRAN yet), I
* added argument `maxstep` to stop FS and BOSS at a specified step size
* extended the estimation of hdf to the scenario of n<=p
* modified function `boss` to account for p>n
* modified function `cv.boss` to account for p>n (only validates subset with sizes up to min(n - n/n.folds, maxstep))
