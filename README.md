# emdr

Empirical mode decomposition regression

## Description

The `emdr` library provides functions to perform EMD-regression as described in Masselot et al. (2018). It includes functions to perform EMD, interpret the resulting IMFs and use them in a suitable regression analysis. See the help (by typing `?emdr` in the console) for details about the functions and their use.

## Installation

1. In R, install the package directly from github using the command (the package `devtools` is required):
```r
> install_github("PierreMasselot/emdr")
```
2. The package can then be loaded as usual: `library(emdr)`.

## References
Huang, N.E., Shen, Z., Long, S.R., Wu, M.C., Shih, H.H., Zheng, Q., Yen, N.-C., Tung, C.C., Liu, H.H., 1998. The empirical mode decomposition and the Hilbert spectrum for nonlinear and non-stationary time series analysis. *Proceedings of the Royal Society of London. Series A: Mathematical, Physical and Engineering Sciences* **454**, 903–995. https://doi.org/10.1098/rspa.1998.0193

Masselot, P., Chebana, F., Bélanger, D., St-Hilaire, A., Abdous, B., Gosselin, P., Ouarda, T.B.M.J., 2018. EMD-regression for modelling multi-scale relationships, and application to weather-related cardiovascular mortality. *Science of The Total Environment* **612**, 1018–1029. https://doi.org/10.1016/j.scitotenv.2017.08.276

