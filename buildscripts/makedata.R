## This document creates the dataset included in the package

setwd("~/code/robustarchitecture")

library("robustarchitecture")
library("devtools")
devtools::load_all(".")

set.seed(4)
n=100000
s=-0.5
testdata=make_test_data(n,s=s)

topdata=top_data(testdata,10000)
save(topdata, file="data/topdata.RData")
