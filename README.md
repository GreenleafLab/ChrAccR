# ChrAccR
* __Package name:__ ChrAccR
* __Title:__ Analyzing chromatin accessibility data in R
* __Description:__ Tools for analyzing chromatin accessibility data in R. Currently supports ATAC-seq and NOMe-seq data analysis.
* __Author:__ Fabian Mueller
* __Maintainer:__ Fabian Mueller <muellerf@stanford.edu>
* __Version:__ 0.9
* __Date:__ 2019-05-24


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
