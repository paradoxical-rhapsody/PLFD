## Dependencies

* `R > 3.5.0`
* If you would like to install `PLFD` from source on win-platform, ensure that the [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is deployed well. 

## Installation

Compiled version of `PLFD` is not yet released on [CRAN](https://cran.r-project.org/). You can install it from source within R via
```R
library(remotes)
install_git("https://gitee.com/xu-zc/PLFD.git")
```
or install with a vignette 
```R
library(remotes)
install_git("https://gitee.com/xu-zc/PLFD.git", build_vignettes=TRUE)
```

## Usage

See the vignettes after installation.

## Reference

A Portmanteau Local Feature Discrimination Approach to the High-Dimensional Matrix-Variate Data
