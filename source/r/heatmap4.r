#############################################################################################################
# RPREP: Ribosomal Profiling Reports, an open-source cloud-enabled framework for reproducible
# Ribosomal Profiling data processing, analysis, and result reporting
# 
# https://github.com/emmesgit/RPREP
#
# Copyright (C) 2020 The Emmes Company L.L.C. 
# 
# This program is free software that contains third party software subject to various licenses, 
# namely, the GNU General Public License version 3 (or later), the GNU Affero General Public License 
# version 3 (or later), and the LaTeX Project Public License v.1.3(c). A list of the software contained 
# in this program, including the applicable licenses, can be accessed here: 
# 
# https://github.com/emmesgit/RPREP/blob/master/SOFTWARE.xlsx  
# 
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# To cite this software, please reference <DOI>
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# Program:  heatmap4.r 
# Version:  RPREP 1.0.1
# Author:   Travis L. Jensen, Johannes B. Goll, Sami R. Cherikh, and William F. Hooper
# Purpose:  getting tired of heatmaps taking up 4 plotting spaces -- making one that plots in a single plotting area.
#           Will allow for multiple simple heatmaps on a page
# Input:    N/A
# Output:   N/A
#############################################################################################################

heatmap4 = function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
        distfun = dist, hclustfun = hclust, reorderfun = function(d, 
                w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, 
                "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
        margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 
                1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
        labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
        verbose = getOption("verbose"), ...) 
{
    scale <- if (symm && missing(scale)) 
                "none"
            else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("'x' must be a numeric matrix")
    nr <- di[1L]
    nc <- di[2L]
    if (nr <= 1 || nc <= 1) 
        stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2L) 
        stop("'margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (!doRdend && identical(Colv, "Rowv")) 
        doCdend <- FALSE
    if (is.null(Rowv)) 
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv)) 
        Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
        if (inherits(Rowv, "dendrogram")) 
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x))
            ddr <- as.dendrogram(hcr)
            if (!is.logical(Rowv) || Rowv) 
                ddr <- reorderfun(ddr, Rowv)
        }
        if (nr != length(rowInd <- order.dendrogram(ddr))) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1L:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram")) 
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr != nc) 
                stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm) 
                                        x
                                    else t(x)))
            ddc <- as.dendrogram(hcc)
            if (!is.logical(Colv) || Colv) 
                ddc <- reorderfun(ddc, Colv)
        }
        if (nc != length(colInd <- order.dendrogram(ddc))) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1L:nc
    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow)) 
                if (is.null(rownames(x))) 
                    (1L:nr)[rowInd]
                else rownames(x)
            else labRow[rowInd]
    labCol <- if (is.null(labCol)) 
                if (is.null(colnames(x))) 
                    (1L:nc)[colInd]
                else colnames(x)
            else labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1L, rowMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 1L, sd, na.rm = na.rm)
        x <- sweep(x, 1L, sx, "/", check.margin = FALSE)
    }
    else if (scale == "column") {
        x <- sweep(x, 2L, colMeans(x, na.rm = na.rm), check.margin = FALSE)
        sx <- apply(x, 2L, sd, na.rm = na.rm)
        x <- sweep(x, 2L, sx, "/", check.margin = FALSE)
    }
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, 4)
    lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
            4)
    if (!missing(ColSideColors)) {
        if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
            stop("'ColSideColors' must be a character vector of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        lhei <- c(lhei[1L], 0.2, lhei[2L])
    }
    if (!missing(RowSideColors)) {
        if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
                        1), lmat[, 2] + 1)
        lwid <- c(lwid[1L], 0.2, lwid[2L])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei, 
                "; lmat=\n")
        print(lmat)
    }
    if (!symm || scale != "none") 
        x <- t(x)
    if (revC) {
        iy <- nr:1
        if (doRdend) 
            ddr <- rev(ddr)
        x <- x[, iy]
    }
    else iy <- 1L:nr
    image(1L:nc, 1L:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
                    c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
    axis(1, 1L:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
            cex.axis = cexCol)
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = 1)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
            cex.axis = cexRow)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = 1)
    if (!missing(add.expr)) 
        eval.parent(substitute(add.expr))
}