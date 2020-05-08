#==========#
# Packages #
#==========#

library("glmpca")
library("colorspace")
library("magrittr")



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


# Graphic variables
# myclr <- c("aquamarine", "aquamarine4", "coral", "coral4")
myclr <- rainbow_hcl(nrow(myconditions))

myseed <- 1598



#=================#
# Data processing #
#=================#

mytype.tmp <- mytypes[1]

# Loading data
myfile <- paste0(result2_fd, "/PZQ-ER-ES.", mytype.tmp, ".counts.matrix")
mydata <- data.matrix(read.table(myfile))

mydata.cln <- mydata[ rowSums(mydata) > 0, ]

set.seed(myseed)
gpca <- glmpca(mydata.cln, L=2)

eig.val <- (gpca$dev^2 / sum(gpca$dev^2) * 100) %>% round(., digit=2)



#=========#
# Figures #
#=========#

if (! dir.exists(graph_fd)) { dir.create(graph_fd) }

myfactors <- gpca$factors

myclr.gpca <- rep(NA, nrow(myfactors))
mypch.gpca <- rep(NA, nrow(myfactors))

for (i in myconditions2) { myclr.gpca[ grep(i, rownames(myfactors)) ] <- myclr[ myconditions2 %in% i ] }
for (i in 1:3) { mypch.gpca[ grep(i, rownames(myfactors)) ] <- i }

pdf(paste0(graph_fd,"/",mytype.tmp,"_PCA.pdf"), width=9)
layout(matrix(1:2, ncol=2), width=c(6.5,3.5))
par(mar=c(4,4,0,0)+0.1)
plot(myfactors, col=myclr.gpca, pch=mypch.gpca, xlab=paste0("Dim 1 (", eig.val[1], "%)"), ylab=paste0("Dim 2 (", eig.val[2], "%)"))
abline(h=0, col="grey")
abline(v=0, col="grey")
# legend("topright", myconditions2, bty="n", col=myclr, pch=1)
plot.new()
par(mar=c(rep(0,4)+0.1))
legend("center", myconditions2, bty="n", col=myclr, pch=1)
dev.off()
