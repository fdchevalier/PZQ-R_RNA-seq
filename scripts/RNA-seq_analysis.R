#!/usr/bin/env Rscript
# Title: RNA-seq_analysis.R
# Version: 0.4
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2020-05-08
# Modified in: 2023-04-24



#==========#
# Comments #
#==========#

# Analyze RNA-seq data from juvenile and adult worms



#==========#
# Versions #
#==========#

# v0.4 - 2023-04-24: adapt the analysis to the S. mansoni v10 genome
# v0.3 - 2021-06-06: file path updated / comparisons added
# v0.2 - 2021-05-13: R object exported / code cleaned
# v0.1 - 2021-01-25: header added / package message handled / info messages added / minor bugs corrected
# v0.0 - 2020-05-08: creation



#==========#
# Packages #
#==========#

cat("Loading packages...\n")

suppressMessages({
    library("EBSeq")
})



#===========#
# Functions #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

source("RNA-seq_analysis_func.R")



#===========#
# Variables #
#===========#

cat("Setting variables...\n")

# Folders
data_fd    <- "../data"
graph_fd   <- "../graphs"
result_fd  <- "../results"
result2_fd <- "../results/3-quantification"


mytypes <- c("gene", "isoform")

# My genes of interest
mygoi3_file <- paste0(result_fd, "/2-QTL/QTL_genes_chr3")
mygenes <- c("Smp_246790", "Smp_317670")

myann_file <- paste0(data_fd, "/genome/Sm_transcript_table_gff-hhpred_2023-04-26.tsv")

# Multi-comparisons
mycomp <- matrix(c("SmLE-PZQ-ER-RNA-adu-m", "SmLE-PZQ-ES-RNA-adu-m",
                   "SmLE-PZQ-ER-RNA-adu-f", "SmLE-PZQ-ES-RNA-adu-f",
                   "SmLE-PZQ-ER-RNA-adu-m", "SmLE-PZQ-ES-RNA-adu-f",
                   "SmLE-PZQ-ES-RNA-adu-m", "SmLE-PZQ-ES-RNA-adu-f",
                   "SmLE-PZQ-ER-RNA-adu-m", "SmLE-PZQ-ES-RNA-juv-m",
                   "SmLE-PZQ-ER-RNA-juv-m", "SmLE-PZQ-ES-RNA-juv-m",
                   "SmLE-PZQ-ER-RNA-juv-f", "SmLE-PZQ-ES-RNA-juv-f",
                   "SmLE-PZQ-ER-RNA-juv-m", "SmLE-PZQ-ER-RNA-adu-m",
                   "SmLE-PZQ-ES-RNA-juv-m", "SmLE-PZQ-ES-RNA-adu-m",
                   "SmLE-PZQ-ER-RNA-juv-f", "SmLE-PZQ-ER-RNA-adu-f",
                   "SmLE-PZQ-ES-RNA-juv-f", "SmLE-PZQ-ES-RNA-adu-f"),
                  ncol=2, byrow=TRUE)

# Graphic variables

## Gene color vector
clr.vec <- c("grey85", "black", "blue", "red", "green")

## Line color for volcano plots
clr.ln <- "brown"



#=================#
# Data processing #
#=================#

cat("Processing data...\n")

# Loading accessory data
mygoi3 <- readLines(mygoi3_file)

mygenes <- c(list(mygoi3), as.list(mygenes))

myann <- read.delim(myann_file, header=TRUE, as.is=TRUE)

myres <- vector("list", length(mytypes))

for (t in 1:length(mytypes)) {

    mytype.tmp <- mytypes[t]

    myres[[t]] <- vector("list", nrow(mycomp))
    names(myres)[t] <- mytype.tmp

    # Loading data
    myfile <- paste0(result2_fd, "/PZQ-ER-ES.", mytype.tmp, ".counts.matrix")
    mydata <- data.matrix(read.table(myfile))

    for (j in 1:nrow(mycomp)) {

        myconditions <- mycomp[j,]

        # Define file name
        myfn <- paste0(c(mytype.tmp, myconditions), collapse="-")

        # Name slot
        names(myres[[t]])[j] <- myfn

        # Subset data
        mycln.names <- make.names(myconditions) # To match column names
        mycln.mat <- sapply(mycln.names, function(x) grep(x, colnames(mydata)))
        mydata.cln <- mydata[, mycln.mat]

        myconditions3 <- rep(NA, length(mycln.mat))

        for (c in 1:ncol(mycln.mat)) {
            mycd.tmp <- paste0("C",c)
            myconditions3[ as.vector(mycln.mat) %in% mycln.mat[,c] ] <- mycd.tmp
        }


        if (mytype.tmp == "isoform") {
            IsoNames     <- rownames(mydata.cln)
            IsoGeneNames <- IsoNames %>% strsplit(., ".", fixed=T) %>% sapply(., function(x) x[1])

            NgList <- GetNg(IsoNames, IsoGeneNames, TrunThre=3)

            IsoNgTrun <- NgList$IsoformNgTrun
        } else {
            IsoNgTrun <- NULL
        }


        myDE.tb <- DE.analysis(mydata.cln, NgVector=IsoNgTrun, Conditions=myconditions3, maxround=10, annotation=myann)

        ## Multi-comparisons
        #myDE.tb <- DE.analysis(mydata.cln, NgVector=IsoNgTrun, Conditions=myconditions2, maxround=10, comparisons=mycomp)

        results <- myDE.tb[[1]]
        myDE    <- myDE.tb[[2]]
        # EBOut   <- myDE.tb[[3]]

        if (! dir.exists(result2_fd)) { dir.create(result2_fd) }
        write.table(as.matrix(results), paste0(result2_fd,"/",myfn,"-all.tsv"), sep="\t", quote=FALSE, col.names=NA)
        write.table(as.matrix(myDE), paste0(result2_fd,"/",myfn,"-DE.tsv"), sep="\t", quote=FALSE, col.names=NA)

        myres[[t]][[j]] <- myDE.tb
    }
}

# Export object for use by other scripts
save(myres, file=paste0(result2_fd, "/myDE.tb.RData"))



#=========#
# Figures #
#=========#

cat("Generating graphs...\n")

if (! dir.exists(graph_fd)) { dir.create(graph_fd) }


for (t in 1:length(mytypes)) {

    mytype.tmp <- mytypes[t]

    for (j in 1:nrow(mycomp)) {
        
        myconditions <- mycomp[j,]

        myfn <- names(myres[[t]])[j]
        
        # Loading data
        myDE.tb <- myres[[t]][[j]]

        myresults <- myDE.tb[[1]]
        EBOut     <- myDE.tb[[3]]

        FCout <- PostFC(EBOut)


        #--------------#
        # Sanity check #
        #--------------#

        pdf(paste0(graph_fd,"/",myfn,"-FC.pdf"))
        PlotPostVsRawFC(EBOut, FCout)
        dev.off()

        pdf(paste0(graph_fd,"/",myfn,"-qqplot.pdf"), width=7*length(myconditions))
        par(mfrow=c(1,length(myconditions)))
        QQP(EBOut)
        DenNHist(EBOut)
        dev.off()

        #--------------#
        # Volcano plot #
        #--------------#

        myclr <- vector("list", length(mygenes)+1)

        myclr <- lapply(myclr, function(x) rep(FALSE, nrow(results)))

        myclr[[1]][ myresults$PPDE > 0.95 & abs(log2(myresults$PostFC)) > 1 ] <- TRUE  # !! make variables !

        for (i in 2:(length(mygenes)+1)) {
            myclr[[i]][ grepl(paste0(mygenes[[i-1]], collapse="|"), rownames(myresults)) ] <- TRUE 
        }

        for (i in 1:length(myclr)) {


            png(paste0(graph_fd,"/",myfn,"-",i,".png"))
            
            plot(log2(myresults$PostFC), myresults$PPDE, xlab="Fold change (log2)", ylab="Posterior probability", pch=19, col=clr.vec[1])
            j <- 1 
            while (j <= i) {
               points(log2(myresults$PostFC)[myclr[[j]]], myresults$PPDE[myclr[[j]]], pch=19, col=clr.vec[j+1])
                j <- j+1
            }

            abline(h=0.95, col=clr.ln, lty=2)
            abline(v=1,    col=clr.ln, lty=2)
            abline(v=-1,   col=clr.ln, lty=2)

            legend("bottomleft", c("QTL chr. 3", "TRP", "TCP1"), col=clr.vec[3:length(clr.vec)], pch=19, bty="n")
            
            dev.off()
        }
    }
}


