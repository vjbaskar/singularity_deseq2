
chooseCRANmirror(graphics = FALSE, ind = 1, local.only = FALSE)
install.packages("devtools")
devtools::install_version("latticeExtra", version="0.6-28")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")


install.packages(c("getopt", "WriteXLS", "dplyr"))
BiocManager::install(c("DESeq2", "rtracklayer", "GenomicFeatures"))
