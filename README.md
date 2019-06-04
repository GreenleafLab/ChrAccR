# ChrAccR <img src="man/figures/chraccr_logo.png" align="right" style="position:absolute; top:0; right:0; padding:10px; height:96px;"/>
* __Package name:__ ChrAccR
* __Title:__ Analyzing chromatin accessibility data in R
* __Description:__ Tools for analyzing chromatin accessibility data in R. Currently supports ATAC-seq and NOMe-seq data analysis.
* __Author/Maintainer:__ Fabian Mueller (<muellerf@stanford.edu>)
* __Version:__ 0.9.1
* __Date:__ 2019-06-04


## Installation

To install `ChrAccR` and its dependencies, use the `devtools` installation routine:

```r
# install devtools if not previously installed
if (!is.element('devtools', installed.packages()[,"Package"])) install.packages('devtools')

# install dependencies
devtools::install_github("demuellae/muLogR")
devtools::install_github("demuellae/muRtools")

# install ChrAccR
devtools::install_github("demuellae/ChrAccR")
```

## Getting started

The `ChrAccR` [vignette](https://muellerf.s3.amazonaws.com/data/ChrAccR/vignette/overview.html) provides a most excellent starting point to get familiar with the package.
