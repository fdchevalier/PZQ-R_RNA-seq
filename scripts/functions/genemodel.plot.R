#'genemodel.plot
#'
#' This function plots a gene model
#'
#' @param model data.frame containing model information. Required columns are "type", "coordinates"
#' @param start start position
#' @param bpstop stop position
#' @param orientation either "foward" or "reverse" indicates the direction of transcription
#' @param xaxis default is TRUE and adds axis above gene model showing position
#' @import graphics
#' @export
#' @examples
#' data(AT5G62640)
#' genemodel.plot(AT5G62640, 25149433, 25152541, "reverse", xaxis=TRUE)

# Source: https://raw.githubusercontent.com/greymonroe/genemodel/c046a58a27402b374471fd6216202c36a5d9e1c0/R/genemodel.plot.R
# Modified in: 2023-08-27


genemodel.plot2 <- function(model, start, bpstop, orientation, col = "steelblue3", xlim, xaxis = TRUE, ypos = NULL, add = FALSE) {

    if (! add) {
        if (xaxis) {
            par(mar=c(1,1,3,1))
        } else {
            par(mar=c(1,1,1,1))
        }
    }

    if (is.null(ypos)) {
        ypos <- 0
    }

    col_hsv <- rgb2hsv(col2rgb(col))
    col1 <- col
    col2 <- hsv(col_hsv[1,], col_hsv[2,] * 0.5, 1)
    col_bd <- hsv(col_hsv[1,], col_hsv[2,], col_hsv[3,] * 0.6)
    
    model <- cbind(model[,1], as.data.frame(stringr::str_split_fixed(model$coordinates, "-", 2)))
    colnames(model) <- c("feature", "start", "bpstop")
    model$start  <- as.numeric(as.character(model$start))
    model$bpstop <- as.numeric(as.character(model$bpstop))

    length <- bpstop-start

    if (orientation == "reverse") {
        model$newstart <- bpstop-model$bpstop+1
        model$newstop  <- bpstop-model$start+1
        model <- model[which(model$feature!="exon"),]
        model <- model[which(model$feature!="ORF"),]
        model <- model[order(model$newstart),]
        
        x <- c(model$newstop[1], model$newstop[1], model$newstart[1], start-.02*length, model$newstart[1])
        y <- c(0,.2,.2,.1,0) + ypos

        if (missing("xlim")) {
            xlim <- c(start-.03*length, bpstop)
        }

    }

    if (orientation == "forward") {
        model$newstart <- start+model$start-1
        model$newstop <- start+model$bpstop-1
        model <- model[which(model$feature!="exon"),]
        model <- model[which(model$feature!="ORF"),]
        model <- model[order(model$newstop, decreasing = T),]
        
        x <- c(model$newstart[1], model$newstart[1], model$newstop[1], bpstop+.02*length, model$newstop[1])
        # y <- c(0,.2,.2,.1,0) + ypos
        y <- c(-0.1, 0.1, 0.1, 0, -0.1) + ypos
        
        if (missing("xlim")) {
            xlim <- c(start, bpstop+.03*length)
        }

    }

    if (! add) {
        plot(1, type = "l", axes = FALSE, ann = FALSE, xlim = xlim, ylim = c(-0.5, 0.5), frame.plot = FALSE)
    }

    for (i in 2:nrow(model)) {
        type <- model$feature[i]

        if (type == "coding_region") {
          rect(model$newstart[i], y[1], model$newstop[i], y[2], col = col1, border = col_bd, lwd = 1)
        }

        if (type == "intron") {
          mid <- mean(c(model$newstart[i], model$newstop[i]))
          segments(x0 = model$newstart[i], y0 = y[4], x1 = mid, y1 = y[2], lwd = 1, col = col_bd)
          segments(x0 = model$newstop[i],  y0 = y[4], x1 = mid, y1 = y[2], lwd = 1, col = col_bd)
        }

        if (type == "3' utr") {
          rect(model$newstart[i], y[1], model$newstop[i], y[2], col = col2, border = col_bd, lwd = 1)
        }

        if (type == "5' utr") {
          rect(model$newstart[i], y[1], model$newstop[i], y[2], col = col2, border = col_bd, lwd = 1)
        }

    }
      
    endtype <- model$feature[1]
    if (endtype == "coding_region") {
        polygon(x, y, border = col_bd , col = col1, lwd = 1) 
    } else {
        polygon(x, y, border = col_bd , col = col2, lwd = 1)
    }

    if (xaxis == T) {
        Axis(side = 3, labels = T, cex.axis = 0.7)
    }
}
