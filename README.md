# PLFD

## Dependencies

* `R > 3.5.0`
* `Rcpp >= 1.0.2`
* `RcppArmadillo >= 0.9.800`

## Installation

On win-platform, ensure that the [Rtools](https://cran.r-project.org/) is deployed well to compile `PLFD` from source.

One can install it from source within R via
```R
library(remotes)
install_git("https://gitee.com/xu-zc/PLFD.git")
```
or install with a vignette 
```R
library(remotes)
install_git("https://gitee.com/xu-zc/PLFD.git", build_vignettes=TRUE)
```

Pre-compiled version of `PLFD` is not yet released on [CRAN](https://cran.r-project.org/) currently. 

## Usage

See the [usage documentation].

## Reference

A Portmanteau Local Feature Discrimination Approach to the High-Dimensional Matrix-Variate Data
