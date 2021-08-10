#!/usr/bin/env Rscript
# Title: Supp. Fig. X2
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-06-06
# Modified in:



#==========#
# Comments #
#==========#

# Generate a figure of global gene and isoform differential expression at different stages for the manuscript



#==========#
# Versions #
#==========#

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


pdf(paste0(graph_fd, "Supp. Fig. X2.pdf"), width = 8, heigh = 4 * nrow(mycomp), useDingbats=TRUE)

layout(matrix(c(1:(nrow(mycomp) * 2)), ncol=2, byrow = FALSE))

for (t in 1:length(mytypes)) {

    mytype.tmp <- mytypes[t]

    for (j in 1:nrow(mycomp)) {
        
        mynm    <- paste0(c(mytype.tmp, mycomp[j,]), collapse="-")
        mytitle <- sapply(mycomp[j,], function(x) myconditions[ myconditions[,1] == x, 2]) %>% paste0(., collapse=" vs. ")
        
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

        plot(log2(myresults$PostFC), myresults$PPDE, xlab="Fold change (log2)", ylab="Posterior probability", pch=19, col=clr.vec[1], main=paste(tools::toTitleCase(mytype.tmp), "expression -", mytitle))
        for (i in 1:length(myclr)) {
            points(log2(myresults$PostFC)[myclr[[i]]], myresults$PPDE[myclr[[i]]], pch=19, col=clr.vec[i+1])
        }

        abline(h=0.95, col=clr.ln, lty=2)
        abline(v=1,    col=clr.ln, lty=2)
        abline(v=-1,   col=clr.ln, lty=2)

        # legend("bottomleft", c("Genes under\nQTL chr. 3", expression(italic("Sm.TRPM"["PZQ"]))), col=clr.vec[3:length(clr.vec)], pch=19, bty="n")
        legend("bottomleft", c("QTL chr. 3", expression(italic("Sm.TRPM"["PZQ"]))), col=clr.vec[3:length(clr.vec)], pch=19, bty="n")
            
    }
    
    mtext(paste0(LETTERS[t], "."), side=3, line=2, at=line2user(par("mar")[2],2), cex=par("cex")*2, adj=0)
}



dev.off()
