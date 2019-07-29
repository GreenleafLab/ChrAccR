# ChrAccR <img src="man/figures/chraccr_logo.png" align="right" height="96"/>

[![Build Status](https://travis-ci.org/GreenleafLab/ChrAccR.svg?branch=master)](https://travis-ci.org/GreenleafLab/ChrAccR)

* __Package name:__ ChrAccR
* __Title:__ Analyzing chromatin accessibility data in R
* __Description:__ Tools for analyzing chromatin accessibility data in R. Currently supports ATAC-seq and NOMe-seq data analysis.
* __Author/Maintainer:__ Fabian Mueller (<muellerf@stanford.edu>)
* __Version:__ 0.9.3
* __Date:__ 2019-07-29


## Installation

To install `ChrAccR` and its dependencies, use the `devtools` installation routine:

```r
# install devtools if not previously installed
if (!is.element('devtools', installed.packages()[,"Package"])) install.packages('devtools')

# install dependencies
devtools::install_github("demuellae/muLogR")
devtools::install_github("demuellae/muRtools")
devtools::install_github("demuellae/muReportR")

# install ChrAccR
devtools::install_github("GreenleafLab/ChrAccR")
```

## Getting started

The `ChrAccR` [vignette](https://GreenleafLab.github.io/ChrAccR/articles/overview.html) provides a most excellent starting point to get familiar with the package.
