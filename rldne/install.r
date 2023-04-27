#! /usr/bin/env -S Rscript --slave --vanilla

install.packages("devtools", repos='https://cran.rstudio.com/', Ncpus=4)

library("devtools")

devtools::install_github(repo="zakrobinson/RLDNe")

library(RLDNe)

fix_permissions_RLDNe(execute=T)
