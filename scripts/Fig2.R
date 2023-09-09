#!/usr/bin/env Rscript
# Title: Fig2.R
# Version: 1.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-05-12
# Modified in: 2023-09-08



#==========#
# Comments #
#==========#

# Generate a figure of global gene and isoform differential expression and specific gene and isoform Sm.TRPM_PZQ expression at different stages for the manuscript



#==========#
# Versions #
#==========#

# v1.1 - 2023-09-08: add gene and isoform models to the figure
# v1.0 - 2023-07-15: rename the script / adapt the code to the v10 genome / improve the figure
# v0.3 - 2021-08-15: adapt the figure to the journal guidelines
# v0.2 - 2021-07-27: rename script / update graph path
# v0.1 - 2021-06-06: update code and file path / correct sample order
# v0.0 - 2021-05-12: creation



#==========#
# Packages #
#==========#

cat("Loading packages...\n")

suppressMessages({
    library("rtracklayer")
    library("magrittr")
    library("gplots")
})

# Turn off factors
options(stringsAsFactors = FALSE)



#===========#
# Functions #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

source("RNA-seq_analysis_func.R")
source("functions/line2user.R")
source("functions/genemodel.plot.R")

# Format coordinates for genemodel.plot
format_coord <- function(x, intron_size = NULL) {
    # Remove exons
    x <- x[ x[, 1] != "exon", ]

    # Create intron
    x_cds <- x[ x[, 1] == "CDS", ]
    n_intron <- nrow(x_cds) - 1
    x_intron <- data.frame(rep("intron", n_intron), rep(NA, n_intron), rep(NA, n_intron))
    colnames(x_intron) <- colnames(x)
    if (is.null(intron_size)) {
        for (i in 2:nrow(x_cds)) {
            x_intron[i - 1, 2:3] <- c(x_cds[i - 1, 3] + 1, x_cds[i, 2] - 1)
        }
    } else {
        for (i in 2:nrow(x_cds)) {
            x_intron[i - 1, 2:3] <- c(x_cds[i - 1, 3] + 1, x_cds[i - 1, 3] + intron_size + 1)
            x_cds[i, 2:3] <- x_cds[i, 2:3] - (x_cds[i, 2] - x_intron[i - 1, 3]) + 1
        }
        x[ x[, 1] == "three_prime_UTR", 2:3] <- x[ x[, 1] == "three_prime_UTR", 2:3] - (x[x[, 1] == "three_prime_UTR", 2] - x_cds[i, 3]) + 1
    }
    x <- rbind(x[ x[, 1] != "CDS", ], x_cds, x_intron) %>% .[order(.[, 2]), ]

    # Rename items
    x[, 1] <- gsub("five_prime_UTR", "5' utr", x[, 1])
    x[, 1] <- gsub("three_prime_UTR", "3' utr", x[, 1])
    x[, 1] <- gsub("CDS", "coding_region", x[, 1])

    # Reformat coordinates
    x <- apply(x, 1, function(y) c(y[1], paste0(y[2:3], collapse = "-"))) %>% t() %>% as.data.frame()
    colnames(x) <- c("type", "coordinates")

    return(x)
}



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

mygff_fl_v7  <- paste0(data_fd, "genome/schistosoma_mansoni.PRJEA36577.WBPS14.annotations.gff3")
mygff_fl_v10 <- paste0(data_fd, "genome/schistosoma_mansoni.PRJEA36577.WBPS18.annotations.gff3")

trpm_domains_fl <- paste0(data_fd, "genome/Smp_246790_domains.tsv")

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

# Load GFF
mygff_v7  <- readGFF(mygff_fl_v7)
mygff_v10 <- readGFF(mygff_fl_v10)


#--------------------#
# V7 SmTRPM_PZQ gene #
#--------------------#

# Gene model
mygene_v7  <- mygff_v7[ grepl(mygenes, mygff_v7[, "Parent"]), ] %>% as.data.frame()

mygene_v7_cd <- mygene_v7[, 3:5] %>% as.data.frame()
mygene_v7_cd <- mygene_v7_cd[ mygene_v7_cd[, 1] == "exon", ] %>% unique() %>% .[order(.[,2]), ]
myshift <- mygene_v7_cd[1, 2] - 1
mygene_v7_cd[, 2:3] <- mygene_v7_cd[, 2:3] - myshift

# Isoform models
myiso_v7 <- (mygene_v7[,"Parent"] %>% unlist() %>% unique())[-1]

myiso_v7_cd <- lapply(myiso_v7, function(x) mygene_v7[mygene_v7[,"Parent"] == x, 3:5]) %>% setNames(., myiso_v7)
myiso_v7_cd <- lapply(myiso_v7_cd, function(x) cbind(x[, 1], x[, 2:3] - myshift))
myiso_v7_cd <- lapply(myiso_v7_cd, format_coord)

myiso_v7_bd <- lapply(myiso_v7_cd, function(x) as.character(x[,2]) %>% strsplit(., "-") %>% unlist() %>% as.numeric() %>% range())

# Domains
trpm_domains <- read.delim(trpm_domains_fl, header = TRUE, stringsAsFactor = FALSE)

## Remove minor domains
trpm_domains <- trpm_domains[ ! grepl("Ankyrin|Pre-|Coil|Ribe|Pore|S5|S6", trpm_domains[,1] ), ]

## Merge S1-S4
trpm_domains[ trpm_domains[, 1] == "S1" , c(3, 5)] <- trpm_domains[ trpm_domains[, 1] == "S4" , c(3, 5)]
trpm_domains <- trpm_domains[ ! grepl("S2|S3|S4", trpm_domains[,1] ), ]
trpm_domains[, 1] <- gsub("S1", "S1-S4", trpm_domains[,1])

# Exons
mygene_v7_cd_ex <- mygene_v7_cd
mygene_v7_cd_ex[, 1] <- as.character(mygene_v7_cd_ex[, 1])
for (i in 2:nrow(mygene_v7_cd)) { mygene_v7_cd_ex[i, 2:3] <- c(mygene_v7_cd_ex[i-1, 3] + 1, mygene_v7_cd_ex[i-1, 3] + 1 + diff(unlist(mygene_v7_cd_ex[i, 2:3]))) }

mygene_v7_cd_ex_ls <- vector("list", nrow(mygene_v7_cd_ex))
for (i in 1:nrow(mygene_v7_cd_ex)) {

    myrows <- NULL
    for (j in 1:nrow(trpm_domains)) {
        if (mygene_v7_cd_ex[i, 2] >= trpm_domains[j, 4] & mygene_v7_cd_ex[i, 3] <= trpm_domains[j, 5]) {
            myrows <- rbind(myrows, as.data.frame(c(trpm_domains[j, ], mygene_v7_cd_ex[i,])))
        } else {
            if (mygene_v7_cd_ex[i, 2] <= trpm_domains[j, 5] & mygene_v7_cd_ex[i, 3] >= trpm_domains[j, 5]) {
                mystart <- c(mygene_v7_cd_ex[i, 2], trpm_domains[j, 4]) %>% max()
                myrow <- as.data.frame(c(trpm_domains[j, ], type = mygene_v7_cd_ex[i, 1], start = mystart, end = trpm_domains[j, 5]))
                myrows <- rbind(myrows, myrow)
            }

            if (mygene_v7_cd_ex[i, 3] >= trpm_domains[j, 4] & mygene_v7_cd_ex[i, 3] <= trpm_domains[j, 5]) {
                myrows <- rbind(myrows, as.data.frame(c(trpm_domains[j, ], type = mygene_v7_cd_ex[i, 1], start = trpm_domains[j, 4], end = mygene_v7_cd_ex[i, 3])))
            }

        }
    }

    if (! is.null(myrows)) {
        mygene_v7_cd_ex_ls[[i]] <- myrows
    }
}


#---------------------#
# V10 SmTRPM_PZQ gene #
#---------------------#

# Gene model
mygene_v10 <- mygff_v10[ grepl(mygenes, mygff_v10[, "Parent"]), ] %>% as.data.frame()

mygene_v10_cd <- mygene_v10[, 3:5] %>% as.data.frame()
mygene_v10_cd <- mygene_v10_cd[ mygene_v10_cd[, 1] == "exon", ] %>% unique() %>% .[order(.[,2]), ]
myshift <- mygene_v10_cd[1, 2] - 1
mygene_v10_cd[, 2:3] <- mygene_v10_cd[, 2:3] - myshift

# Isoform model
myiso_v10 <- (mygene_v10[,"Parent"] %>% unlist() %>% unique())[-1]

myiso_v10_cd <- lapply(myiso_v10, function(x) mygene_v10[mygene_v10[,"Parent"] == x, 3:5]) %>% setNames(., myiso_v10)
myiso_v10_cd <- lapply(myiso_v10_cd, function(x) cbind(x[, 1], x[, 2:3] - myshift))
myiso_v10_cd <- lapply(myiso_v10_cd, format_coord)

myiso_v10_bd <- lapply(myiso_v10_cd, function(x) as.character(x[,2]) %>% strsplit(., "-") %>% unlist() %>% as.numeric() %>% range())


# Exon differences between V7 and V10
exon_v7_lg  <- apply(mygene_v7_cd[, 2:3], 1, diff)
exon_v10_lg <- apply(mygene_v10_cd[, 2:3], 1, diff)
exon_diff   <- ! exon_v7_lg %in% exon_v10_lg


#-----------------#
# Expression data #
#-----------------#

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

pdf(paste0(graph_fd, "Fig. 5 v4.pdf"), width = 5.5, heigh = 8, useDingbats = TRUE)

layout(matrix(c(1, 1, 2, 3, 4, 5, 6, 7), ncol = 2, byrow = TRUE), heights = c(0.4, 0.48, 0.48, 0.04))


#-------------------------#
# Gene and isoform models #
#-------------------------#

# Graphic parameters
par(mar=c(2, 4, 2, 0) + 0.1)
myltr <- 1

mymax <- lapply(c(myiso_v7_cd[1:5], myiso_v10_cd), function(x) as.character(x[,2]) %>% strsplit(., "-") %>% unlist() %>% as.numeric() %>% max()) %>% unlist() %>% max()
xlim  <- c(mygene_v7_cd[1, 2], mymax)
ylim  <- c(-length(c(myiso_v7_cd, myiso_v10_cd)) - 1, 0) / 2
plot(1, type = "l", axes = FALSE, ann = FALSE, xlim = xlim, ylim = ylim, frame.plot = FALSE)

# Altered exons between v7 and v10
for(i in 1:nrow(mygene_v7_cd)) {
    if (exon_diff[i]) {
        rect(mygene_v7_cd[i,2], par("usr")[3], mygene_v7_cd[i,3], par("usr")[4], col = "grey", border = "grey")
    }
}

# Gene model
abline(h = 0)

trpm_domains_cd <- trpm_domains[,1] %>% unique() %>% cbind(., ., ., .)
trpm_domains_cd[, -1] <- NA
trpm_domains_cd %<>% as.data.frame()

domain <- ""
for(i in 1:nrow(mygene_v7_cd)) {
    if (is.null(mygene_v7_cd_ex_ls[[i]])) {
        rect(mygene_v7_cd[i,2], -0.1, mygene_v7_cd[i,3], 0.1, col = "black", border = "black")
    } else {
        # Get coordinates for adjustement
        myexon_start  <- mygene_v7_cd[i,2] - 1
        myexon_cd_rel <- mygene_v7_cd_ex[i,2] - 1

        for(j in 1:nrow(mygene_v7_cd_ex_ls[[i]])) {
            # Select data
            myrow <- mygene_v7_cd_ex_ls[[i]][j, ]

            # Ajust coordinates
            mycoord <- myrow[, 8:9] - myexon_cd_rel + myexon_start

            # Adjust color
            myclr <- myrow[,6]
            if (myclr == "white") {
                mybdr <- "black"
            } else {
                mybdr <- myclr
            }

            # Draw rectangle
            rect(mycoord[,1], -0.1, mycoord[,2], 0.1, col = myclr, border = mybdr)

            if (domain != myrow[,1]) {
                domain <- myrow[,1]
                trpm_domains_cd[ trpm_domains_cd[,1] == domain, 2:4] <- c(mycoord, mybdr)
            } else {
                trpm_domains_cd[ trpm_domains_cd[,1] == domain, 3] <- mycoord[,2]
            }

        }
    }
}

# Add annotations
for (i in 1:nrow(trpm_domains_cd)) {
    myshift <- line2user(0, 3)

    lines(trpm_domains_cd[i,2:3], rep(par("usr")[4], 2) + myshift, col = trpm_domains_cd[i,4], xpd = TRUE)
    my_x <- trpm_domains_cd[i,2:3] %>% as.numeric %>% mean
    text(my_x, par("usr")[4] * 2 + myshift, trpm_domains_cd[i,1], col = trpm_domains_cd[i,4], xpd = TRUE)
}

# Isoforms v7
for (i in 1:length(myiso_v7_cd)) {
    genemodel.plot2(model = myiso_v7_cd[[i]], start = 1, bpstop = myiso_v7_bd[[i]][2], orientation = "forward", xaxis = FALSE, ypos = -i / 2, add = TRUE)
}
idx <- i + 1

# Isoforms v10
for (i in 1:length(myiso_v10_cd)) {
    genemodel.plot2(model = myiso_v10_cd[[i]], start = 1, bpstop = myiso_v10_bd[[i]][2], orientation = "forward", xaxis = FALSE, col = "orange", ypos = -(idx + i) / 2, add = TRUE)
}

mylabels <- c("gene", 1:length(myiso_v7_cd), "", 1:length(myiso_v10_cd))
axis(2, at = seq(ylim[2], ylim[1], length.out = length(mylabels)), labels = mylabels, lwd = 0, las = 1)
title(ylab = "Isoform")

mtext(LETTERS[myltr], side = 3, line = 1, at = line2user(par("mar")[2],2), cex = par("cex")*2, adj = 0, padj = 0.5)


#-----------------------------#
# Gene and isoform expression #
#-----------------------------#

# Original plot
for (t in 1:length(mytypes)) {

    par(mar=c(5,4,4,2)+0.1)

    mytype.tmp <- mytypes[t]
    myltr <- myltr + 1

    for (j in 1:nrow(mycomp)) {

        mynm <- paste0(c(mytype.tmp, mycomp[j,]), collapse="-")

        # Loading data
        myDE.tb <- myDE_all[[t]][[mynm]]
        myresults <- myDE.tb[[1]]

        #~~~~~~~~~~~~~~#
        # Volcano plot #
        #~~~~~~~~~~~~~~#

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


    mtext(LETTERS[myltr], side = 3, line = 2, at = line2user(par("mar")[2],2), cex=par("cex")*2, adj = 0, padj = 0.5)

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

    mtext(LETTERS[myltr], side = 3, line = 2, at = line2user(par("mar")[2],2), cex = par("cex")*2, adj = 0, padj = 0.5)
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
