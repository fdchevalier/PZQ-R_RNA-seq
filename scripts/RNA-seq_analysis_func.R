library("EBSeq")


DE.analysis <- function(x, NgVector=NULL, Conditions=NULL, maxround=10, FDR=0.05, fold.change=2, comparisons=NULL, annotation=NULL) {
    
    # Arguments:
    # x             data frame with counts
    # NgVector      NgVector (see doc of EBseq)
    # Conditions    Conditions (see doc of EBseq)
    # maxround      Number of iteration (see doc of EBseq). The best iteration number will be determined
    # FDR           False discovery rate level (see doc of EBseq)
    # fold.change   Threshold of fold change in expression
    # comparisons   Vector or matrix for comparisons in case of multi condition comparisons
    # cond.names    Include conditions names instead of C1, C2, etc.
    # annotation    Table of annotation

    # TODO
    # * Sanity check of object given
    # * Improvment of the second section about multi comparisons (especially the object returned)


    # Normalisation factor (size of the libraries)
    Sizes <- MedianNorm(x)

    if (length(unique(Conditions)) < 3) {

        # Test iteration number to reach convergence
        myit <- 1
        while(myit < 3) {
            # Run model with the given maxround
            EBOut <- EBTest(x, NgVector=NgVector, Conditions=Conditions, sizeFactors=Sizes, maxround=maxround, Print=FALSE)

            # Compute difference in beta values between iteration
            mybeta <- apply(EBOut$Beta, 2, function(x) log10(abs(diff(x))))

            # Optimal number of iteration for convergence
            maxround.tmp <- which(rowSums(mybeta < -3) == ncol(mybeta)) + 1

            if (length(maxround.tmp) == 0) {
                if (any(rowSums(mybeta < -3) == ncol(mybeta))) { stop("Model did not converged after ", nrow(mybeta)+1, " iterations. Please increase the maxround value.") }
                maxround <- maxround+maxround
                myit <- myit+1
            } else {
                message("Convergence reached after ", maxround.tmp[1], " iterations.")
                maxround <- maxround.tmp[1]
                myit <- Inf
            }

        }

        EBOut <- EBTest(x, NgVector=NgVector, Conditions=Conditions, sizeFactors=Sizes, maxround=maxround, Print=FALSE)
        
        PP <- as.data.frame(GetPPMat(EBOut))
        fc_res <- PostFC(EBOut)

        results <- cbind(PP, fc_res$PostFC, fc_res$RealFC,unlist(EBOut$C1Mean)[rownames(PP)], unlist(EBOut$C2Mean)[rownames(PP)])

        cond.names <- colnames(x)[sapply(unique(Conditions), function(x) which(Conditions == x))[1,]]
        colnames(results) <- c("PPEE", "PPDE", "PostFC", "RealFC", paste0(cond.names[1],"_Mean"),paste0(cond.names[2],"_Mean"))

        myDE <- results[ results$PPDE > 1-FDR & abs(log2(results$PostFC)) > log2(fold.change), ]
        if (! is.null(annotation)) {
            mynames <- sapply(rownames(myDE), function(x) strsplit(x, ":")[[1]][2])
            myDE    <- cbind(myDE, t(sapply(mynames, function(x) annotation[ grep(x, annotation[,1]), 2:3][1,])))
        }

        return(list(results, myDE, EBOut))

    } else {

        if (is.null(comparisons)) { stop("comparisons argument missing. Exiting.") }

        mypatterns <- GetPatterns(Conditions)

        mycomp2 <- apply(comparisons, 2, make.names)
        
        # Check if conditions correspond to column names of x

        mypatterns.nb <- matrix(NA, ncol=length(Conditions), nrow=nrow(mycomp2))
        for (i in 1:nrow(mycomp2)) {
            for (j in 1:ncol(mycomp2)) {
                 mypatterns.nb[i, grep(mycomp2[i,j], Conditions) ] <- j
            }
        }

        Parti <- mypatterns[ apply(mypatterns, 1, function(x) paste0(x, collapse="")) %in%  apply(mypatterns.nb, 1, function(x) paste0(x, collapse="")) , ]

        MultiOut <- EBMultiTest(x, NgVector=NgVector, Conditions=Conditions, AllParti=Parti, sizeFactors=Sizes, maxround=maxround)

        MultiPP <- GetMultiPP(MultiOut)
        MultiFC <- GetMultiFC(MultiOut)

        return(MultiOut)
    }

}



