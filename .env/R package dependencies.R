# System
options("repos" = c(CRAN = "https://cran.revolutionanalytics.com"))
library("devtools")
library("BiocManager")

# PCA analysis
install_github("willtownes/glmpca", ref="edc04cc", upgrade="never")
