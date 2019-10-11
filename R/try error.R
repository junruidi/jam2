rm(list = ls())
library(jam2)
load("~/Dropbox/Junrui Di/fosr study/data/processed.rda")
id = ids
Y = act[,1:20]
rm(list = setdiff(ls(),c("Y","id")))

mf_pca = multilevel_pca(Y = Y, id = id, twoway = FALSE, cov.method = "m1", pve = 0.9,smoothing = T, smooth.method = "sc")


