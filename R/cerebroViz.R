#' A tool to visualize biological data mapped to SVG brain diagrams.
#'
#' \code{cerebroViz} is a tool for visualizing spatiotemporal data in the brain within an anatomical context. The user inputs a matrix and the tool creates SVG diagrams with color mapping reflective of the input data. \code{cerebroViz} supports 30 brain regions used by BrainSpan, GTEx, Roadmap Epigenomics, and more.
#' @param x matrix containing input data. Row names should reflect the appropriate brain regions. Column names may represent different time points or samples and do not require naming
#' @param filePrefix desired prefix for the output SVG files
#' @param palette character vector of color values for visualizing brain regions. Accepts color names, hex values, and RGB values
#' @param timePoint a numeric vector of columns in 'x' to visualize
#' @param divData logical indicating if input data is to be visualized as it diverges from the median. Defaults to false
#' @param secPalette character vector of length three specifying colors for the brain background, brain outline/legend, and svg background in that order
#' @param clamp coefficient to the Median Absolute Deviation. Added and subtracted from the median to identify a range of non-outliers. Values external to this range will 'clamped' to extremes of the non-outlier range
#' @param naHatch logical indicating if regions of missing data should be filled with a cross-hatch pattern to differentiate them from the brain's background
#' @param legend logical indicating if the legend bar is displayed
#' @param customNames dataframe or matrix with 2 columns. The first column for cerebroViz convention names, the second column for custom user names. Cells in input cannot be factors.
#' @param figLabel logical indicating if figure label should be added to the output.
#' @param regLabel logical indicating if region labels should be added to the output.
#' @keywords cerebroViz
#' @import XML
#' @export
#' @examples
#' data(cerebroEx)
#' cerebroViz(cerebroEx)
cerebroViz <- function(x, filePrefix = "cerebroViz_output", palette = NULL,
    timePoint = 1, divData = FALSE, secPalette = c("white", "black", "white"),
    clamp = NULL, naHatch = FALSE, legend = TRUE, customNames = NULL, figLabel = FALSE, regLabel = FALSE) {
    require(XML)

    if (is.null(palette) & divData == FALSE) {
        palette <- c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E", "#7A0177", "#49006A")
    }

    if (is.null(palette) & divData == TRUE) {
        palette <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
    }

    regions <- c("A1C", "AMY", "ANG", "BS", "CAU", "CB", "CNG", "DFC",
        "FL", "HIP", "HTH", "IPC", "ITC", "M1C", "MED", "MFC", "OL", "OFC",
        "PL", "PIT", "PUT", "PON", "S1C", "SN", "STC", "STR", "TL", "THA",
        "V1C", "VFC")

    if (class(x) != "matrix")
        stop("'x' must be of class 'matrix'.")
    if (sum(is.na(rownames(x))) > 0)
        stop("Row names of 'x' must be valid.")
    if (length(rownames(x)) != nrow(x))
        stop("Row names must be supplied for each row in 'x'.")
    if (max(timePoint) > ncol(x))
        stop("'timePoint' invalid")
    if (sum(timePoint%%1 != 0))
        stop("'timePoint' invalid")
    if (length(secPalette) != 3)
        stop("'secPalette' must have length 3.")
    if (!is.null(customNames)) {
        if (ncol(customNames) != 2)
            stop("Unexpected input for customNames.")
    }
    if (is.null(customNames) & sum(rownames(x) %in% regions == FALSE) >
        0)
        warning(paste("Unknown row names in input data: ", paste(rownames(x)[rownames(x) %in%
            regions == FALSE], collapse = ", "), ". Unknown regions will be excluded from visualization. See the help manual for 'customNames' argument.",
            sep = ""))


    if (!is.null(customNames)) {
        if (sum(customNames[, 2] %in% regions) > 0)
            stop(paste("customNames contains region names already used by cerebroViz convention: ",
                paste(customNames[customNames[, 2] %in% regions, 2], collapse = ", "),
                sep = ""))
        for (indA in 1:nrow(customNames)) {
            rownames(x)[rownames(x) == customNames[indA, 2]] <- customNames[indA,
                1]
        }
    }

    inpReg <- rownames(x)

    # creating the vector for 'parent' regions (regions that encompass
    # others) and a warning of overshadowing.
    supReg <- c("BS", "FL", "OL", "PL", "TL", "STR")
    supLog <- matrix(!is.na(x[rownames(x) %in% supReg, timePoint]), ncol = length(timePoint))
    rownames(supLog) <- rownames(x)[rownames(x) %in% supReg]
    userReg <- rowSums(supLog)
    userReg <- names(userReg[userReg > 0])
    if (length(userReg) > 0) {
        warning(paste("The following regions encompass other regions of the brain: ",
            paste(userReg, collapse = ", "), ". Subregions are masked in output.",
            sep = ""))
    }

    xmed <- median(x, na.rm = TRUE)
    xmad <- mad(x, constant = 1, na.rm = TRUE)
    xmin <- min(x, na.rm = TRUE)
    xmax <- max(x, na.rm = TRUE)

    # set the default clamp value (no clamping)
    avoidClamp <- max(abs(xmed - xmin), abs(xmed - xmax))/xmad
    if (is.null(clamp)) {
        clamp <- avoidClamp + 1
    }
    if (clamp <= 0)
        stop("'clamp' must be >0")

    naMatrix <- matrix(data = NA, nrow = length(regions[which(regions %in%
        rownames(x) == FALSE)]), ncol = ncol(x))
    rownames(naMatrix) <- regions[which(regions %in% rownames(x) == FALSE)]
    dataMatrix <- rbind(x, naMatrix)
    dataMatrix <- as.matrix(dataMatrix[order(rownames(dataMatrix)), ])
    xScaled <- cerebroScale(dataMatrix, clamp, divData)
    hexInd <- round(xScaled * 200 + 1)
    f <- colorRampPalette(palette)
    hexVec <- f(201)

    # loop for each timePoint
    svgOuter <- system.file("extdata/svg/brain_outer.svg", package = "cerebroViz")
    svgSlice <- system.file("extdata/svg/brain_slice.svg", package = "cerebroViz")
    for (indA in 1:length(timePoint)) {
        if (ncol(x) > 1) {
            lupiter <- hexInd[, timePoint[indA]]
        } else {
            lupiter <- hexInd[,1]
        }
        oXml <- xmlTreeParse(svgOuter, useInternalNodes = TRUE)
        sXml <- xmlTreeParse(svgSlice, useInternalNodes = TRUE)
        cXml <- c(oXml, sXml)
        cXml <- editSecPal(cXml, secPalette)
        if (naHatch == TRUE) {
            cXml <- editHatch(cXml, lupiter)
        }
        cXml <- editPal(lupiter, cXml, hexVec, naHatch)
        cXml <- maskRegions(cXml, supReg, lupiter)
        cXml <- editLegend(cXml, palette, divData, clamp, legend, xmin,
            xmed, xmad, xmax)
        cXml <- editLabel(cXml, timePoint, figLabel, x, indA, regLabel, inpReg)
        oXml <- cXml[1][[1]]
        sXml <- cXml[2][[1]]
        saveXML(oXml, paste(filePrefix, "_outer_", timePoint[indA], ".svg",
            sep = ""))
        saveXML(sXml, paste(filePrefix, "_slice_", timePoint[indA], ".svg",
            sep = ""))
    }
    message("Success! Your diagrams have been saved.")
}

#' A function to scale sequential and divergent data to a 0:1 range
#'
#' This function scales input data to a 0:1 range. If divData=FALSE it will scale linearly. If divData=TRUE, it will utilize a polylinear scale centered at the median of the input data.
#' @param x input data matrix
#' @param clamp coefficient to the Median Absolute Deviation. Added and subtracted from the median to identify a range of non-outliers. Values external to this range will 'clamped' to extremes of the non-outlier range.
#' @param divData logical indicating if input data is divergent in nature. Default assumes data is sequential.
#' @keywords cerebroScale
#' @export
#' @examples
#' cerebroScale(x, clamp = NULL, divData=FALSE)
# cerebroScale
cerebroScale <- function(x, clamp, divData) {
    xmed <- median(x, na.rm = TRUE)
    xmad <- mad(x, constant = 1, na.rm = TRUE)
    xmin <- min(x, na.rm = TRUE)
    xmax <- max(x, na.rm = TRUE)
    avoidClamp <- max(abs(xmed - xmin), abs(xmed - xmax))/xmad
    fillMatrix <- x

    if (is.null(clamp)) {
        clamp <- avoidClamp + 1
    }
    outlrs <- clamp * xmad
    if (clamp <= 0)
        stop("clamp must be >0")
    pctOL <- round(length(which(x[!is.na(x)] <= (xmed - (outlrs)) | x[!is.na(x)] >=
        (xmed + (outlrs))))/length(x[!is.na(x)]) * 100, 2)

    if (pctOL > 0) {
        warning(paste("The clamp value of ", clamp, " will clamp ", pctOL,
            "% of input values (outliers) to the min or max of the scaled range.",
            sep = ""))
    }

    if (divData == TRUE) {
        abvMed <- x[x >= xmed & x <= (xmed + outlrs) & !is.na(x)]
        belMed <- x[x <= xmed & x >= (xmed - outlrs) & !is.na(x)]
        if (length(which(!is.na(x)))%%2 == 0) {
            # imputing median if even number of data points
            rightsc <- rescale(c(xmed, abvMed), c(0.5, 1))[-1]
            fillMatrix[x >= xmed & x <= (xmed + outlrs) & !is.na(x)] <- rightsc
            leftsc <- rescale(c(xmed, belMed), c(0, 0.5))[-1]
            fillMatrix[x <= xmed & x >= (xmed - outlrs) & !is.na(x)] <- leftsc
            fillMatrix[x < (xmed - outlrs) & !is.na(x)] <- 0
            fillMatrix[x > (xmed + outlrs) & !is.na(x)] <- 1
            xScaled <- fillMatrix
        }
        if ((length(which(!is.na(x))))%%2 == 1) {
            rightsc <- rescale(abvMed, c(0.5, 1))
            fillMatrix[x >= xmed & x <= (xmed + outlrs) & !is.na(x)] <- rightsc
            leftsc <- rescale(belMed, c(0, 0.5))
            fillMatrix[x <= xmed & x >= (xmed - outlrs) & !is.na(x)] <- leftsc
            fillMatrix[x < (xmed - outlrs) & !is.na(x)] <- 0
            fillMatrix[x > (xmed + outlrs) & !is.na(x)] <- 1
            xScaled <- fillMatrix
        }
    }
    if (divData == FALSE) {
        nonoutlrs <- x[x >= (xmed - outlrs) & x <= (xmed + outlrs) & !is.na(x)]
        xsc <- rescale(nonoutlrs, c(0, 1))
        fillMatrix[x >= (xmed - outlrs) & x <= (xmed + outlrs) & !is.na(x)] <- xsc
        fillMatrix[x < (xmed - outlrs) & !is.na(x)] <- 0
        fillMatrix[x > (xmed + outlrs) & !is.na(x)] <- 1
        xScaled <- fillMatrix
    }
    return(xScaled)
}

#' a function used by cerebroViz() to edit the brain outline and brain background.
#'
#' for each xml, remove fill attributes from 'brainBackground' & 'brainOutline' and replace them with the designated colors.
#' @param cXml list containing the xml object for each SVG.
#' @param secPalette character vector of length three specifying colors for the brain outline, brain background, and svg background in that order.
#' @keywords internal
#' @examples
#' editSecPal(cXml, secPalette)
# editSecPal
editSecPal <- function(cXml, secPalette) {
    for (indB in 1:length(cXml)) {
        node <- getNodeSet(cXml[indB][[1]], "//*[@id='brainBackground']")[[1]]
        removeAttributes(node, "fill")
        addAttributes(node, fill = col2hex(secPalette[1]))
        node <- getNodeSet(cXml[indB][[1]], "//*[@id='brainOutline']")[[1]]
        removeAttributes(node, "fill")
        addAttributes(node, fill = col2hex(secPalette[2]))
        node <- getNodeSet(cXml[indB][[1]], "//*[@id='legendRect']")[[1]]
        removeAttributes(node, "stroke")
        addAttributes(node, stroke = col2hex(secPalette[2]))
        node <- getNodeSet(cXml[indB][[1]], "//*[@id='leglableft']")[[1]]
        removeAttributes(node, "fill")
        addAttributes(node, fill = col2hex(secPalette[2]))
        node <- getNodeSet(cXml[indB][[1]], "//*[@id='leglabmid']")[[1]]
        removeAttributes(node, "fill")
        addAttributes(node, fill = col2hex(secPalette[2]))
        node <- getNodeSet(cXml[indB][[1]], "//*[@id='leglabright']")[[1]]
        removeAttributes(node, "fill")
        addAttributes(node, fill = col2hex(secPalette[2]))
        node <- getNodeSet(cXml[indB][[1]], "//*[@id='svgBackground']")[[1]]
        removeAttributes(node, "fill")
        addAttributes(node, fill = col2hex(secPalette[3]))
    }
    return(cXml)
}

#' a function used by cerebroViz() to add cross-hatching to regions with missing data.
#'
#' for each xml, get style and append cross-hatching pattern.
#' @param cXml list containing the xml object for each SVG.
#' @param lupiter hex gradient indices for current iteration (timePoint) in the loop, specified within cerebroViz().
#' @keywords internal
#' @examples
#' editHatch(cXml, lupiter)
# editHatch
editHatch <- function(cXml, lupiter) {
    missNames <- names(lupiter[is.na(lupiter)])
    for (indB in 1:length(cXml)) {
        for (indC in 1:length(missNames)) {
            node <- getNodeSet(cXml[indB][[1]], paste("//*[@id='", missNames[indC],
                "']", sep = ""))[1]
            if (!is.null(node[[1]])) {
                style <- xmlGetAttr(node[[1]], "style")
                style <- paste("fill:url(#hatch00);", style, sep = "")
                addAttributes(node[[1]], style = style)
            }
        }
    }
    return(cXml)
}

#' a function used by cerebroViz() to map input data to brain regions.
#'
#' for each xml, add the appropriate fill attributes for each region
#' @param lupiter hex gradient indices for current iteration (timePoint) in the loop, specified within cerebroViz().
#' @param cXml list containing the xml object for each SVG.
#' @param hexVec character vector of length 201 containing the hex value gradient for visualization.
#' @param naHatch logical indicating if regions of missing data should be filled with a cross-hatch pattern to differentiate them from the brain's background.
#' @keywords internal
#' @examples
#' editPal(lupiter, cXml, hexVec, naHatch)
# editPal
editPal <- function(lupiter, cXml, hexVec, naHatch) {
    fillData <- lupiter[which(!is.na(lupiter))]
    for (indB in 1:length(cXml)) {
        for (indC in 1:length(fillData)) {
            node <- getNodeSet(cXml[indB][[1]], paste("//*[@id='", names(fillData)[indC],
                "']", sep = ""))[1]
            if (!is.null(node[[1]])) {
                removeAttributes(node[[1]], "fill")
                addAttributes(node[[1]], fill = hexVec[fillData[indC]])
            }
        }
    }
    if (naHatch == FALSE) {
        missNames <- lupiter[is.na(lupiter)]
        for (indB in 1:length(cXml)) {
            for (indC in 1:length(missNames)) {
                node <- getNodeSet(cXml[indB][[1]], paste("//*[@id='",
                  names(missNames)[indC], "']", sep = ""))[1]
                if (!is.null(node[[1]])) {
                  removeAttributes(node[[1]], "fill-opacity")
                  addAttributes(node[[1]], "fill-opacity" = 0)
                }
            }
        }
    }
    return(cXml)
}

#' A function used by cerebroViz() to edit the legend.
#'
#' For each xml, map the appropriate colors to the legend.
#' @param cXml list containing the xml object for each SVG.
#' @param palette character vector of color values to use in the visualization. Accepts color names, hex values, and RGB values. For sequential data, it is recommended to use two colors, or a sequence of colors in a gradient. For divergent data, it is recommended to use three colors with a neutral color in the middle.
#' @param divData logical indicating if input data is divergent in nature. Default assumes data is sequential.
#' @param clamp coefficient to the Median Absolute Deviation. Added and subtracted from the median to identify a range of non-outliers. Values external to this range will 'clamped' to extremes of the non-outlier range.
#' @param legend logical indicating if the legend bar should be visible.
#' @param xmin minimum of input data.
#' @param xmed median of input data.
#' @param xmad Median Absolute Deviation of input data (constant=1).
#' @param xmax maximum of input data.
#' @keywords internal
#' @examples
#' editLegend(cXml, palette, divData, clamp, legend, xmin, xmed, xmad, xmax)
editLegend <- function(cXml, palette, divData, clamp, legend, xmin, xmed,
    xmad, xmax) {
    labmin <- round(max(xmin, (xmed - (clamp * xmad))), 3)
    labmax <- round(min(xmax, (xmed + (clamp * xmad))), 3)
    labmed <- round(xmed, 3)
    labels <- c(labmin, labmed, labmax)
    stopoffset <- paste(as.character(seq(0, 100, (100/(length(palette) -
        1)))), "%", sep = "")
    stopcolor <- col2hex(palette)
    for (indB in 1:length(cXml)) {
        gradnode <- getNodeSet(cXml[indB][[1]], "//*[@id='gradient']")
        for (indC in 1:length(palette)) {
            newstop <- newXMLNode("stop", attrs = c(offset = stopoffset[indC],
                "stop-color" = stopcolor[indC]))
            gradnode[[1]] <- addChildren(gradnode[[1]], newstop)
        }
        for (indC in 1:length(labels)) {
            node <- getNodeSet(cXml[indB][[1]], "//*[@class='legendLabel']")[[indC]]
            nv <- paste("\n", labels[indC], "\n", sep = "")
            xmlValue(node) <- nv
            if (divData == FALSE & indC == 2) {
                nv <- paste("\n", "\n", sep = "")
                xmlValue(node) <- nv
            }
        }
        if (legend == FALSE) {
            node <- getNodeSet(cXml[indB][[1]], "//*[@class='legendBar']")[[1]]
            removeAttributes(node, "opacity")
            addAttributes(node, opacity = "0")
        }
    }
    return(cXml)
}

#' A function used by cerebroViz() to lower opacity for superior regions when data is not provided.
#'
#' For each missing superior region, set opacity to 0.
#' @param cXml list containing the xml object for each SVG.
#' @param supReg regions that encompass other brain regions, specified within cerebroViz().
#' @param lupiter hex gradient indices for current iteration (timePoint) in the loop, specified within cerebroViz().
#' @keywords internal
#' @examples
#' maskRegions(cXml, supReg, lupiter)
# maskRegions
maskRegions <- function(cXml, supReg, lupiter) {
    missNames <- names(lupiter[is.na(lupiter)])
    opacdown <- supReg[supReg %in% missNames]
    if (length(opacdown) > 0) {
        for (indB in 1:length(opacdown)) {
            lobename <- opacdown[indB]
            lobenode <- getNodeSet(cXml[1][[1]], paste("//*[@id='", lobename,
                "']", sep = ""))
            if (length(lobenode) > 0) {
                node <- lobenode[[1]]
                removeAttributes(node, "fill-opacity")
                addAttributes(node, "fill-opacity" = "0")
            }
        }
    }
    if (("STR" %in% missNames) & (sum(c("CAU", "PUT") %in% missNames) <
        2)) {
        node <- getNodeSet(cXml[2][[1]], "//*[@id='STR']")[[1]]
        removeAttributes(node, "fill-opacity")
        addAttributes(node, "fill-opacity" = "0")
    }
    return(cXml)
}

#' A function used by cerebroViz() to add a figure label.
#'
#' For each xml, get label text node a fill with column names or timePoint.
#' @param cXml list containing the xml object for each SVG.
#' @param timePoint a numeric vector of columns in 'x' to visualize.
#' @param figLabel logical indicating if figure label should be added to the output.
#' @param x matrix containing input data. Row names should reflect the appropriate brain regions. Column names may represent different time points or samples and do not require naming.
#' @param indA level 'A' for loop index.
#' @param regLabel logical indicating if region labels should be added to the output.
#' @param inpReg names of input regions.
#' @keywords internal
#' @examples
#' editLabel(cXml, timePoint, figLabel, x, indA, regLabel, inpReg)
editLabel <- function(cXml, timePoint, figLabel, x, indA, regLabel, inpReg) {
    if (figLabel == TRUE) {
        if (!is.null(colnames(x))) {
            for (indB in 1:length(cXml)) {
                node <- getNodeSet(cXml[indB][[1]], "//*[@id='labelText']")[1]
                nv <- paste("\n", colnames(x)[timePoint[indA]], "\n", sep = "")
                xmlValue(node[[1]]) <- nv
            }
        } else {
            for (indB in 1:length(cXml)) {
                node <- getNodeSet(cXml[indB][[1]], "//*[@id='labelText']")[1]
                nv <- paste("\n", timePoint[indA], "\n", sep = "")
                xmlValue(node[[1]]) <- nv
            }
        }
    }
    if (regLabel == TRUE) {
        labClass <- paste(inpReg, "_lab", sep="")
        for (indB in 1:length(cXml)) {
            for (indC in 1:length(labClass)) {
                node <- getNodeSet(cXml[indB][[1]], paste("//*[@class='", labClass[indC],"']", sep=""))
                if (length(node) == 2) {
                    for (indD in 1:length(node)) {
                        removeAttributes(node[[indD]], "opacity")
                        addAttributes(node[[indD]], opacity = "1")
                    }
                }
            }
            if ("STR_lab" %in% labClass & indB == 2) {
                node <- getNodeSet(cXml[indB][[1]], "//*[@id='STR_children']")[[1]]
                removeAttributes(node, "opacity")
                addAttributes(node, opacity = "0")
            }
            if ("FL_lab" %in% labClass & indB == 1) {
                node <- getNodeSet(cXml[indB][[1]], "//*[@id='FL_children']")[[1]]
                removeAttributes(node, "opacity")
                addAttributes(node, opacity = "0")
            }
            if ("PL_lab" %in% labClass & indB == 1) {
                node <- getNodeSet(cXml[indB][[1]], "//*[@id='PL_children']")[[1]]
                removeAttributes(node, "opacity")
                addAttributes(node, opacity = "0")
            }
            if ("OL_lab" %in% labClass & indB == 1) {
                node <- getNodeSet(cXml[indB][[1]], "//*[@id='OL_children']")[[1]]
                removeAttributes(node, "opacity")
                addAttributes(node, opacity = "0")
            }
            if ("TL_lab" %in% labClass & indB == 1) {
                node <- getNodeSet(cXml[indB][[1]], "//*[@id='TL_children']")[[1]]
                removeAttributes(node, "opacity")
                addAttributes(node, opacity = "0")
            }
            if ("BS_lab" %in% labClass & indB == 1) {
                node <- getNodeSet(cXml[indB][[1]], "//*[@id='BS_children']")[[1]]
                removeAttributes(node, "opacity")
                addAttributes(node, opacity = "0")
            }
        }
    }
    return(cXml)
}

#' A function from the 'gplots' package used by cerebroViz() to convert color values to hex values.
#'
#' For each input color, convert to a hex value.
#' @param cname color name or value.
#' @keywords internal
#' @examples
#' col2hex(cname = "red")
col2hex <- function(cname)
  {
    colMat <- col2rgb(cname)
    rgb(
        red=colMat[1,]/255,
        green=colMat[2,]/255,
        blue=colMat[3,]/255
        )
  }

  #' A function from the 'gplots' package used by cerebroViz() to convert color values to hex values.
  #'
  #' For each input color, convert to a hex value.
  #' @param resc_inp numeric vector of values to manipulate.
  #' @param to output range (numeric vector of length two).
  #' @param from input range (numeric vector of length two).  If not given, is calculated from the range of ‘resc_inp’.
  #' @keywords internal
  #' @examples
  #' rescale(1:100)
rescale <- function (resc_inp, to = c(0, 1), from = range(resc_inp, na.rm = TRUE, finite = TRUE))
  {
    if (zero_range(from) || zero_range(to)) {
        return(ifelse(is.na(resc_inp), NA, mean(to)))
    }
    (resc_inp - from[1])/diff(from) * diff(to) + to[1]
  }
