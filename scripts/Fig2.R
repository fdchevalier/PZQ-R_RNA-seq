#!/usr/bin/env Rscript
# Title: Fig2.R
# Version: 1.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-05-12
# Modified in: 2023-07-15



#==========#
# Comments #
#==========#

# Generate a figure of global gene and isoform differential expression and specific gene and isoform Sm.TRPM_PZQ expression at different stages for the manuscript



#==========#
# Versions #
#==========#

# v1.0 - 2023-07-15: rename the script / adapt the code to the v10 genome / improve the figure
# v0.3 - 2021-08-15: adapt the figure to the journal guidelines
# v0.2 - 2021-07-27: rename script / update graph path
# v0.1 - 2021-06-06: update code and file path / correct sample order
# v0.0 - 2021-05-12: creation



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
mycomp <- matrix(c("SmLE-PZQ-ER-RNA-adu-m", "SmLE-PZQ-ES-RNA-adu-m"),
                  ncol=2, byrow=TRUE)

# Graphic variables

## Gene color vector
clr.vec <- c("grey85", "black", "blue", "red")

## Line color for volcano plots
clr.ln <- "brown"

## Population colors
mycolor <- c("#d8c99c", "#8cd6cf") %>% setNames(., c("ER", "ES"))

myconditions <- matrix(c("SmLE-PZQ-ES-RNA-adu-m", "ES adult male",
                         "SmLE-PZQ-ES-RNA-adu-f", "ES adult female",
                         "SmLE-PZQ-ES-RNA-juv-m", "ES juvenile male",
                         "SmLE-PZQ-ES-RNA-juv-f", "ES juvenile female",
                         "SmLE-PZQ-ER-RNA-adu-m", "ER adult male",
                         "SmLE-PZQ-ER-RNA-adu-f", "ER adult female",
                         "SmLE-PZQ-ER-RNA-juv-m", "ER juvenile male",
                         "SmLE-PZQ-ER-RNA-juv-f", "ER juvenile female"),
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


pdf(paste0(graph_fd, "Fig. 5 v2.pdf"), width = 7.5, heigh = 7.5, useDingbats = TRUE)

#layout(matrix(c(1, 2, 3, 4), ncol=2, byrow=TRUE))
layout(matrix(c(1, 2, 3, 4, 5, 6), ncol = 2, byrow = TRUE), heights = c(0.48, 0.48, 0.04))

myltr <- 0

for (t in 1:length(mytypes)) {

    par(mar=c(5,4,4,2)+0.1)
    
    mytype.tmp <- mytypes[t]
    myltr <- myltr + 1

    for (j in 1:nrow(mycomp)) {
        
        mynm <- paste0(c(mytype.tmp, mycomp[j,]), collapse="-")
       
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

        plot(log2(myresults$PostFC), myresults$PPDE, xlab=expression("Fold change (log"[2] * ")"), ylab="Posterior probability", pch=19, col=clr.vec[1], main=paste("Global", mytype.tmp, "expression"))
        for (i in 1:length(myclr)) {
            points(log2(myresults$PostFC)[myclr[[i]]], myresults$PPDE[myclr[[i]]], pch=19, col=clr.vec[i+1])
        }

        abline(h=0.95, col=clr.ln, lty=2)
        abline(v=1,    col=clr.ln, lty=2)
        abline(v=-1,   col=clr.ln, lty=2)
    }
    
    
    mtext(LETTERS[myltr], side=3, line=2, at=line2user(par("mar")[2],2), cex=par("cex")*2, adj=0)

    par(mar=c(5,4,4,1)+0.1)

    myltr <- myltr + 1

    myclr.bx <- myconditions[,2]
    myclr.bx <- gsub(".*ER.*", mycolor["ER"], myclr.bx)
    myclr.bx <- gsub(".*ES.*", mycolor["ES"], myclr.bx)

    myarg <- myconditions[,2]
    myarg <- gsub("ER |ES ", "", myarg)
    
    mydata.i.mean <- mycount[[t]][[1]]
    mydata.i.sd   <- mycount[[t]][[2]]

    if (mytype.tmp == "gene") {
        mydata.i.mean <- mydata.i.mean[reo_count]
        mydata.i.sd   <- mydata.i.sd[reo_count]
        mytitle.tmp   <- mytype.tmp
    }

    if (mytype.tmp == "isoform") {
        myiso <- 1
        mydata.i.mean <- mycount[[t]][[1]][myiso,][reo_count]
        mydata.i.sd   <- mycount[[t]][[2]][myiso,][reo_count]
        mytitle.tmp   <- paste(mytype.tmp, myiso)
    }
    barplot2(t(mydata.i.mean), space=c(0, 0.5), ylab="Expression level (TPM)", beside=TRUE, plot.ci=TRUE, ci.u=t(mydata.i.mean+mydata.i.sd), ci.l=t(mydata.i.mean-mydata.i.sd), col=myclr.bx, names.arg=NULL, main=bquote(italic("Sm.TRPM"["PZQ"]) ~ .(mytitle.tmp) ~ "expression"))
    mypos <- seq(1, par("xaxp")[2], par("xaxp")[2] / length(myarg))
    text(mypos, par("usr")[3] * 5, srt = 45, adj = 1, labels = myarg, xpd = TRUE)

    mtext(LETTERS[myltr], side=3, line=2, at=line2user(par("mar")[2],2), cex=par("cex")*2, adj=0)
}

# Legend
## Legend volcano plot
par(mar = rep(0, 4))
plot(NULL, xlim = c(0,1), ylim = c(0,1), axes = FALSE, bty = "n", xlab = "", ylab = "")
legend("center", c("QTL chr. 3", expression(italic("Sm.TRPM"["PZQ"]))), col=clr.vec[3:length(clr.vec)], pch=19, bty="n", horiz = TRUE)

## Legend barplot
plot(NULL, xlim = c(0,1), ylim = c(0,1), axes = FALSE, bty = "n", xlab = "", ylab = "")
legend("center", legend = rev(names(mycolor)), fill = rev(mycolor), bty = "n", horiz = TRUE)


invisible(dev.off())
