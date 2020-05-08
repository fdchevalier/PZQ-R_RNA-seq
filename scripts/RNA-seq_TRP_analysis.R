
#==========#
# Packages #
#==========#

library("EBSeq")
library("matrixStats")
library("gplots")
library("colorspace")



#===========#
# Variables #
#===========#

# Working directory
#setwd(file.path(getwd(), "scripts"))

# Folders
data_fd    <- "../data"
graph_fd   <- "../graphs"
result_fd  <- "../results"
result2_fd <- "../results/3-quantification"


mytypes <- c("gene", "isoform")
mygenes <- c("Smp_246790")

# myconditions <- c("SmLE-PZQ-ER-RNA-f", "SmLE-PZQ-ER-RNA-m", "SmLE-PZQ-ES-RNA-f", "SmLE-PZQ-ES-RNA-m")
# myconditions <- c("SmLE-PZQ-ER-RNA-juv-f", "SmLE-PZQ-ER-RNA-juv-m", "SmLE-PZQ-ES-RNA-juv-f", "SmLE-PZQ-ES-RNA-juv-m", "SmLE-PZQ-ER-RNA-adu-f", "SmLE-PZQ-ER-RNA-adu-m", "SmLE-PZQ-ES-RNA-adu-f", "SmLE-PZQ-ES-RNA-adu-m")
myconditions <- matrix(c("SmLE-PZQ-ER-RNA-juv-f", "ER juvenile female",
                         "SmLE-PZQ-ER-RNA-juv-m", "ER juvenile male",
                         "SmLE-PZQ-ES-RNA-juv-f", "ES juvenile female",
                         "SmLE-PZQ-ES-RNA-juv-m", "ES juvenile male",
                         "SmLE-PZQ-ER-RNA-adu-f", "ER adult female",
                         "SmLE-PZQ-ER-RNA-adu-m", "ER adult male",
                         "SmLE-PZQ-ES-RNA-adu-f", "ES adult female",
                         "SmLE-PZQ-ES-RNA-adu-m", "ES adult male"),
                       ncol=2, byrow=TRUE)
myconditions2 <- make.names(myconditions[,1]) # To match column names

# myclr <- c("aquamarine", "aquamarine4", "coral", "coral4")
myclr <- rainbow_hcl(nrow(myconditions))


#=================#
# Data processing #
#=================#

myres <- vector("list", length(mytypes))

for (t in 1:length(mytypes)) {

    mytype.tmp <- mytypes[t]


    # Load data
    myfile <- paste0(result2_fd, "/PZQ-ER-ES.", mytype.tmp, ".counts.matrix")
    mydata <- read.delim(myfile, header=TRUE, row.names=1)


    # Normalization
    size <- MedianNorm(mydata)
    mydata <- GetNormalizedMat(mydata, size)

    mydata.i <- mydata[ grepl(mygenes, rownames(mydata)), , drop=FALSE]

    # sapply(seq(1,ncol(mydata.i), by=3), function(x) mean(mydata.i[,x:(x+2)]))

    mydata.i.mean <- sapply(myconditions2, function(x) rowMeans(mydata.i[,grep(paste0(x,".*"), colnames(mydata.i)), drop=FALSE]))
    mydata.i.sd   <- sapply(myconditions2, function(x) rowSds(mydata.i[,grep(paste0(x,".*"), colnames(mydata.i)), drop=FALSE]) / sqrt(length(grep(paste0(x,".*"), colnames(mydata.i)))))

    myres[[t]] <- list(mydata.i.mean, mydata.i.sd)
    names(myres)[t] <- mytype.tmp

}



#=========#
# Figures #
#=========#

if (! dir.exists(graph_fd)) { dir.create(graph_fd) }

# Graphs
for (t in 1:length(mytypes)) {

    mytype.tmp <- mytypes[t]

    mydata.i.mean <- myres[[t]][[1]]
    mydata.i.sd <- myres[[t]][[2]]
    
    pdf(paste0(graph_fd,"/",mygenes,"_",mytype.tmp,"_TPM.pdf"), width=3+length(myconditions))
    par(mar=c(8,4,0,1)+0.1)

    if (mytype.tmp == "gene") {
        barplot2(t(mydata.i.mean), ylab="Expression level (TPM)", beside=TRUE, plot.ci=TRUE, ci.u=t(mydata.i.mean+mydata.i.sd), ci.l=t(mydata.i.mean-mydata.i.sd), las=2, col="grey", names.arg=myconditions[,2])
    }

    if (mytype.tmp == "isoform") {
        mynames <- sapply(rownames(mydata.i.mean), function(x) strsplit(x, ":")[[1]][[2]])
        barplot2(t(mydata.i.mean), ylab="Expression level (TPM)", beside=TRUE, plot.ci=TRUE, ci.u=t(mydata.i.mean+mydata.i.sd), ci.l=t(mydata.i.mean-mydata.i.sd), las=2, col=myclr, names.arg=mynames)
        legend("topleft", myconditions[,2], fill=myclr, bty="n")
    }
    dev.off()
}
