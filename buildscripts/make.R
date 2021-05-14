## Required packages
library("devtools")
library(roxygen2)

## Document and check
setwd("~/code")
setwd("robustarchitecture")
document()
check()

## Install package locally
setwd("~/code")
install("robustarchitecture")
setwd("robustarchitecture")

## Remember that users should install via:
remotes::install_github("danjlawson/robustarchitecture")
## Unless for some reason you wish to modify the code.
