#!/usr/bin/env Rscript
# Title: Supp_Fig5.R
# Version: 0.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-06-06
# Modified in: 2021-08-10



#==========#
# Comments #
#==========#

# Generate a figure of global gene and isoform differential expression at different stages for the manuscript



#==========#
# Versions #
#==========#

# v0.1 - 2021-08-10: rename script / update graph path / redesign graph layout
# v0.0 - 2021-06-06: creation



#===========#
# Functions #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

source("RNA-seq_analysis_func.R")
source("functions/line2user.R")



#===========#
# Variables #
#===========#

cat("Setting variables...\n")

# Folders
data_fd    <- "../data/"
graph_fd   <- "../graphs/"
result_fd  <- "../results/"
result2_fd <- "../results/3-quantification/"


mytypes <- c("gene", "isoform")
mygenes <- c("Smp_246790")

# My genes of interest
mygoi3_file <- paste0(result_fd, "/2-QTL/QTL_genes_chr3")
mygenes <- c("Smp_246790")

myann_file <- paste0(data_fd, "/genome/Sm_v7.1_transcript_table_gff-hhpred.tsv")

# Multi-comparisons
mycomp <- matrix(c("SmLE-PZQ-ER-RNA-adu-m", "SmLE-PZQ-ES-RNA-adu-f",
                   "SmLE-PZQ-ES-RNA-adu-m", "SmLE-PZQ-ES-RNA-adu-f",
                   "SmLE-PZQ-ES-RNA-juv-m", "SmLE-PZQ-ES-RNA-adu-m",
                   "SmLE-PZQ-ES-RNA-juv-f", "SmLE-PZQ-ES-RNA-adu-f",
                   "SmLE-PZQ-ER-RNA-juv-m", "SmLE-PZQ-ES-RNA-juv-m",
                   "SmLE-PZQ-ER-RNA-juv-f", "SmLE-PZQ-ES-RNA-juv-f",
                   "SmLE-PZQ-ER-RNA-adu-m", "SmLE-PZQ-ES-RNA-juv-m"),
                  ncol=2, byrow=TRUE)

# Graphic variables

## Gene color vector
clr.vec <- c("grey85", "black", "blue", "red")

## Line color for volcano plots
clr.ln <- "brown"

myconditions <- matrix(c("SmLE-PZQ-ES-RNA-adu-m", "ES adult male",
                         "SmLE-PZQ-ES-RNA-adu-f", "ES adult female",
                         "SmLE-PZQ-ES-RNA-juv-m", "ES juvenile male",
                         "SmLE-PZQ-ES-RNA-juv-f", "ES juvenile female",
                         "SmLE-PZQ-ER-RNA-juv-m", "ER juvenile male",
                         "SmLE-PZQ-ER-RNA-juv-f", "ER juvenile female",
                         "SmLE-PZQ-ER-RNA-adu-m", "ER adult male",
                         "SmLE-PZQ-ER-RNA-adu-f", "ER adult female"),
                       ncol=2, byrow=TRUE)
myconditions2 <- make.names(myconditions[,1]) # To match column names



#=================#
# Data processing #
#=================#

cat("Processing data...\n")

# Loading accessory data
mygoi3 <- readLines(mygoi3_file)

mygenes <- c(list(mygoi3), as.list(mygenes))

load(paste0(result2_fd, "/myDE.tb.RData"))
myDE_all <- myres

load(paste0(result2_fd, "/Normalized_count.RData"))
mycount <- myres

# Reorder vector
reo_count <- sapply(myconditions2, function(x) grep(x, names(mycount[[1]][[1]]))) %>% unname()



#=========#
# Figures #
#=========#

cat("Generating graphs...\n")

if (! dir.exists(graph_fd)) { dir.create(graph_fd) }

pdf(paste0(graph_fd, "Supp. Fig. 5.pdf"), width = 8, heigh = 3 * nrow(mycomp), useDingbats=TRUE)

mycex <- 1.5

# Layout matrix
## Core matrix
mymat <- matrix(c(1:(nrow(mycomp) * 2)), ncol=2, byrow = FALSE)

## Add rows
mymat <- rbind(max(mymat) + 1:2, mymat)
mymat <- rbind(mymat, rep(max(mymat) + 1, 2))
mymat <- rbind(mymat, rep(max(mymat) + 1, 2))

## Add columns
mymat <- cbind(c(0, rep(max(mymat) + 1, nrow(mycomp)), 0, 0), mymat)
mymat <- cbind(mymat, c(0, max(mymat) + 1:nrow(mycomp), 0, 0))

## Create layout
layout(mymat, widths = c(0.1, 1, 1, 0.5), heights = c(0.1, rep(1, nrow(mycomp)), 0.1, 0.2))

# Define margins
par(mar = c(2, 2, 1, 1))

for (t in 1:length(mytypes)) {

    mytype.tmp <- mytypes[t]

    for (j in 1:nrow(mycomp)) {
        
        mynm    <- paste0(c(mytype.tmp, mycomp[j,]), collapse="-")
        
        # Loading data
        myDE.tb <- myDE_all[[t]][[mynm]]
        myresults <- myDE.tb[[1]]

        #--------------#
        # Volcano plot #
        #--------------#

        myclr <- vector("list", length(mygenes)+1)

        myclr <- lapply(myclr, function(x) rep(FALSE, nrow(myresults)))

        myclr[[1]][ myresults$PPDE > 0.95 & abs(log2(myresults$PostFC)) > 1 ] <- TRUE  # !! make variables !

        for (i in 2:(length(mygenes)+1)) {
            myclr[[i]][ grepl(paste0(mygenes[[i-1]], collapse="|"), rownames(myresults)) ] <- TRUE 
        }

        # Axis
        if (t == 2) { myyax <- "n" } else { myyax <- "s" }
        if (j != nrow(mycomp)) { myxax <- "n" } else { myxax <- "s" }

        plot(log2(myresults$PostFC), myresults$PPDE, xlab="", ylab="", pch=19, col=clr.vec[1], xaxt=myxax, yaxt=myyax)
        for (i in 1:length(myclr)) {
            points(log2(myresults$PostFC)[myclr[[i]]], myresults$PPDE[myclr[[i]]], pch=19, col=clr.vec[i+1])
        }

        abline(h=0.95, col=clr.ln, lty=2)
        abline(v=1,    col=clr.ln, lty=2)
        abline(v=-1,   col=clr.ln, lty=2)

    }
    
}

# Redefine margins
par(mar = c(0, 2, 0, 1))

# Titles
for (i in mytypes) {
    plot.new()
    text(0.5, 0.5, labels = paste(tools::toTitleCase(i), "expression"), font = 2, cex = mycex * 1.5)
}

# X axis
plot.new()
text(0.5, 0.5, labels = "Fold change (log2)", cex = mycex)

# Legend
plot.new()
myorder <- matrix(1:length(clr.vec), ncol = 2, byrow = TRUE)
legend("center", c("Non significant differential expression", "Significant differential expression", "Gene / isoform expressed under QTL chr. 3", expression(italic(Sm.TRPM[PZQ]) ~ expression))[myorder], col = clr.vec[myorder], pch = 19, bty = "n", ncol = 2)

# Redefine margins
par(mar = rep(0, 4))

# Y axis
plot.new()
text(0.5, 0.5, labels = "Posterior probability", cex = mycex, srt = 90)

# Comparison titles
for(j in 1:nrow(mycomp)) {
    plot.new()
    mytitle <- sapply(mycomp[j,], function(x) myconditions[ myconditions[,1] == x, 2]) %>% paste0(., collapse="\nvs.\n")
    text(0.5, 0.5, labels = mytitle, cex = mycex, font = 2)
}

dev.off()
