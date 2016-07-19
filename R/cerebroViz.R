#' A tool to visualize biological data mapped to SVG brain diagrams.
#'
#' 'cerebroViz' is a tool for visualizing spatiotemporal data in the brain within an anatomical context. The user inputs a matrix and the tool creates SVG diagrams with color mapping reflective of the input data. 'cerebroViz' supports 30 brain regions used by BrainSpan, GTEx, Roadmap Epigenomics, and more.
#' @param x matrix containing input data. Row names should reflect the appropriate brain regions. Column names may represent different time points or samples and do not require naming
#' @param filePrefix desired prefix for the output SVG files
#' @param palette character vector of color values for visualizing brain regions. Accepts color names, hex values, and RGB values
#' @param timePoint a numeric vector of columns in 'x' to visualize
#' @param divData logical indicating if input data is to be visualized as it diverges from the median. Defaults to false
#' @param secPalette character vector of length three specifying colors for the brain background, brain outline/legend, and svg background in that order
#' @param clamp coefficient to the Median Absolute Deviation. Added and subtracted from the median to identify a range of non-outliers. Values external to this range will 'clamped' to extremes of the non-outlier range
#' @param naHatch logical indicating if regions of missing data should be filled with a cross-hatch pattern to differentiate them from the brain's background
#' @param legend logical indicating if the legend bar is displayed
#' @param customNames dataframe or matrix with 2 columns. The first column for cerebroViz convention names, the second column for custom user names. Cells in input cannot be factors
#' @keywords cerebroViz
#' @import XML
#' @import gplots
#' @import scales
#' @import grDevices
#' @import RColorBrewer
#' @export
#' @examples
#' x = t(apply(apply(rbind(matrix((sample(c(-400:600),260)/100),nrow=26,ncol=10),matrix(NA,nrow=4,ncol=10)),2,sample),1,sample))
#' rownames(x) = c("A1C", "CNG", "AMY", "ANG", "BS", "CAU", "CB", "DFC", "FL", "HIP", "HTH", "IPC", "ITC", "M1C", "MED", "MFC", "OL", "OFC", "PL", "PIT", "PUT", "PON", "S1C", "SN", "STC", "STR", "TL", "THA", "V1C", "VFC")
#' cerebroViz(x)
cerebroViz <- function(x, filePrefix = "cerebroViz_output", palette = NULL, timePoint = 1, divData = FALSE, secPalette = c("white","black","white"), clamp = NULL, naHatch = FALSE, legend = TRUE, customNames = NULL){
  require(XML)
  require(gplots)
  require(scales)
  require(grDevices)
  require(RColorBrewer)

  if(is.null(palette) & divData==FALSE){
    palette = brewer.pal(n=9, name="YlOrRd")
  }

  if(is.null(palette) & divData==TRUE){
    palette = rev(brewer.pal(n=11, name="RdYlBu"))
  }

  #creating the master regions vector
  regions = c("A1C", "AMY", "ANG", "BS", "CAU", "CB", "CNG", "DFC", "FL", "HIP", "HTH", "IPC", "ITC", "M1C", "MED", "MFC", "OL", "OFC", "PL", "PIT", "PUT", "PON", "S1C", "SN", "STC", "STR", "TL", "THA", "V1C", "VFC")
################################################ E R R O R   H A N D L I N G ###
  if(class(x)!="matrix") stop("'x' must be of class 'matrix'.")
  if(sum(is.na(rownames(x)))>0) stop("Row names of 'x' must be valid.")
  if(length(rownames(x))!=nrow(x)) stop("Row names must be supplied for each row in 'x'.")
  if(max(timePoint)>ncol(x)) stop("'timePoint' invalid")
  if(sum(timePoint%%1!=0)) stop("'timePoint' invalid")
  if(length(secPalette)!=3) stop("'secPalette' must have length 3.")
  if(!is.null(customNames)){
    if(ncol(customNames)!=2) stop("Unexpected input for customNames.")
  }
  if(is.null(customNames) & sum(rownames(x)%in%regions==FALSE)>0) warning(paste("Unknown row names in input data: ",paste(rownames(x)[rownames(x)%in%regions==FALSE],collapse=", "),". Unknown regions will be excluded from visualization. See the help manual for 'customNames' argument.",sep=""))

  #customNames
  if(!is.null(customNames)){
    if(sum(customNames[,2]%in%regions)>0) stop(paste("customNames contains region names already used by cerebroViz convention: ",paste(customNames[customNames[,2]%in%regions,2],collapse=", "),sep=""))
    for(indA in 1:nrow(customNames)){
      rownames(x)[rownames(x)==customNames[indA,2]]=customNames[indA,1]
    }
  }

#################################################### R E G I O N   S E T U P ###
  #creating the vector for 'parent' regions (regions that encompass others) and a warning of overshadowing.
  supReg = c("BS", "FL", "OL", "PL", "TL", "STR")
  supLog = matrix(!is.na(x[rownames(x) %in% supReg,timePoint]),ncol=length(timePoint))
  rownames(supLog) = rownames(x)[rownames(x)%in%supReg]
  userReg = rowSums(supLog)
  userReg = names(userReg[userReg>0])
  if(length(userReg)>0){
     warning(paste("The following regions encompass other regions of the brain: ", paste(userReg, collapse=", "),". Subregions are masked in output.", sep=""))
  }

  xmed = median(x, na.rm=TRUE)
  xmad = mad(x, constant = 1, na.rm=TRUE)
  xmin = min(x, na.rm=TRUE)
  xmax = max(x, na.rm=TRUE)

  #set the default clamp value (no clamping)
  avoidClamp = max(abs(xmed-xmin),abs(xmed-xmax))/xmad
  if(is.null(clamp)){
    clamp = avoidClamp+1
  }
  if(clamp<=0) stop("'clamp' must be >0")

  naMatrix = matrix(data=NA, nrow=length(regions[which(regions%in%rownames(x)==FALSE)]),ncol=ncol(x))
  rownames(naMatrix) = regions[which(regions%in%rownames(x)==FALSE)]
  dataMatrix = rbind(x, naMatrix)
  dataMatrix = as.matrix(dataMatrix[order(rownames(dataMatrix)),])
  xScaled = cerebroScale(dataMatrix, clamp, divData)
  hexInd = round(xScaled*200+1)
  f = colorRampPalette(palette)
  hexVec = f(201)

  #loop for each timePoint
  svgOuter = system.file("extdata/svg/brain-outer.svg",package="cerebroViz")
  svgSlice = system.file("extdata/svg/brain-slice.svg",package="cerebroViz")
  for(indA in 1:length(timePoint)){
    if(ncol(x)>1) {
      lupiter = hexInd[,timePoint[indA]]
    }
    else {
      lupiter = hexInd
    }
    oXml = xmlTreeParse(svgOuter, useInternalNodes=TRUE)
    sXml = xmlTreeParse(svgSlice, useInternalNodes=TRUE)
    cXml =  c(oXml, sXml)
    cXml = editSecPal(cXml, secPalette)
    if(naHatch==TRUE){
        cXml = editHatch(cXml, lupiter)
    }
    cXml = editPal(lupiter, cXml, hexVec, naHatch)
    cXml = maskRegions(cXml, supReg, lupiter)
    cXml = editLegend(cXml, palette, divData, clamp, legend, xmin, xmed, xmad, xmax)
    oXml = cXml[1][[1]]
    sXml = cXml[2][[1]]
    saveXML(oXml, paste(filePrefix,"_outer_",timePoint[indA],".svg",sep=""))
    saveXML(sXml, paste(filePrefix,"_slice_",timePoint[indA],".svg",sep=""))
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
#' cerebroScale(x, clamp = 100, divData=FALSE)
#cerebroScale
cerebroScale <- function(x, clamp, divData){
  xmed = median(x, na.rm=TRUE)
  xmad = mad(x, constant = 1, na.rm=TRUE)
  xmin = min(x, na.rm=TRUE)
  xmax = max(x, na.rm=TRUE)
  avoidClamp = max(abs(xmed-xmin),abs(xmed-xmax))/xmad
  fillMatrix = x

  if(is.null(clamp)){
    clamp = avoidClamp+1
  }
  outlrs = clamp*xmad
  if(clamp<=0) stop("clamp must be >0")
  pctOL = round(length(which(x[!is.na(x)]<=(xmed-(outlrs)) | x[!is.na(x)]>=(xmed+(outlrs))))/length(x[!is.na(x)])*100,2)

  if(pctOL>0){
    warning(paste("The clamp value of ", clamp," will clamp ",pctOL,"% of input values (outliers) to the min or max of the scaled range.", sep=""))
  }

  if(divData==TRUE){
    abvMed = x[x>=xmed & x<=(xmed+outlrs) & !is.na(x)]
    belMed = x[x<=xmed & x>=(xmed-outlrs) & !is.na(x)]
    if(length(which(!is.na(x))) %% 2 == 0){ #imputing median if even number of data points
      rightsc = rescale(c(xmed,abvMed),c(0.5,1))[-1]
      fillMatrix[x>=xmed & x<=(xmed+outlrs) & !is.na(x)] = rightsc
      leftsc = rescale(c(xmed,belMed),c(0,0.5))[-1]
      fillMatrix[x<=xmed & x>=(xmed-outlrs) & !is.na(x)] = leftsc
      fillMatrix[x<(xmed-outlrs) & !is.na(x)] = 0
      fillMatrix[x>(xmed+outlrs) & !is.na(x)] = 1
      xScaled = fillMatrix
    }
    if((length(which(!is.na(x)))) %% 2 == 1){
      rightsc = rescale(abvMed,c(0.5,1))
      fillMatrix[x>=xmed & x<=(xmed+outlrs) & !is.na(x)] = rightsc
      leftsc = rescale(belMed,c(0,0.5))
      fillMatrix[x<=xmed & x>=(xmed-outlrs) & !is.na(x)] = leftsc
      fillMatrix[x<(xmed-outlrs) & !is.na(x)] = 0
      fillMatrix[x>(xmed+outlrs) & !is.na(x)] = 1
      xScaled = fillMatrix
    }
  }
  if(divData==FALSE){
    nonoutlrs = x[x>=(xmed-outlrs) & x<=(xmed+outlrs) & !is.na(x)]
    xsc = rescale(nonoutlrs,c(0,1))
    fillMatrix[x>=(xmed-outlrs) & x<=(xmed+outlrs) & !is.na(x)] = xsc
    fillMatrix[x<(xmed-outlrs) & !is.na(x)] = 0
    fillMatrix[x>(xmed+outlrs) & !is.na(x)] = 1
    xScaled = fillMatrix
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
#editSecPal
editSecPal <- function(cXml, secPalette){
  for(indB in 1:length(cXml)){
    node = getNodeSet(cXml[indB][[1]], "//*[@id='brainBackground']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(secPalette[1]))
    node = getNodeSet(cXml[indB][[1]], "//*[@id='brainOutline']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(secPalette[2]))
    node = getNodeSet(cXml[indB][[1]], "//*[@id='legendRect']")[[1]]
    removeAttributes(node, "stroke")
    addAttributes(node, stroke=col2hex(secPalette[2]))
    node = getNodeSet(cXml[indB][[1]], "//*[@id='leglableft']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(secPalette[2]))
    node = getNodeSet(cXml[indB][[1]], "//*[@id='leglabmid']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(secPalette[2]))
    node = getNodeSet(cXml[indB][[1]], "//*[@id='leglabright']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(secPalette[2]))
    node = getNodeSet(cXml[indB][[1]], "//*[@id='svgBackground']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(secPalette[3]))
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
#editHatch
editHatch <- function(cXml, lupiter){
  missNames = names(lupiter[is.na(lupiter)])
  for(indB in 1:length(cXml)){
    for(indC in 1:length(missNames)){
      node = getNodeSet(cXml[indB][[1]], paste("//*[@id='",missNames[indC],"']",sep=""))[1]
      if(!is.null(node[[1]])){
        style = xmlGetAttr(node[[1]], "style")
        style = paste("fill:url(#hatch00);",style, sep="")
        addAttributes(node[[1]], style=style)
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
#editPal
editPal <- function(lupiter, cXml, hexVec, naHatch){
  fillData = lupiter[which(!is.na(lupiter))]
  for(indB in 1:length(cXml)){
    for(indC in 1:length(fillData)){
      node = getNodeSet(cXml[indB][[1]], paste("//*[@id='",names(fillData)[indC],"']",sep=""))[1]
      if(!is.null(node[[1]])){
        removeAttributes(node[[1]], "fill")
        addAttributes(node[[1]], fill=hexVec[fillData[indC]])
      }
    }
  }
  if(naHatch==FALSE){
    missNames = lupiter[is.na(lupiter)]
    for(indB in 1:length(cXml)){
      for(indC in 1:length(missNames)){
        node = getNodeSet(cXml[indB][[1]], paste("//*[@id='",names(missNames)[indC],"']",sep=""))[1]
        if(!is.null(node[[1]])){
          removeAttributes(node[[1]], "fill-opacity")
          addAttributes(node[[1]], "fill-opacity"=0)
        }
      }
    }
  }
  return(cXml)
}

#' a function used by cerebroViz() to edit the legend.
#'
#' for each xml, map the appropriate colors to the legend
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
editLegend <- function(cXml, palette, divData, clamp, legend, xmin, xmed, xmad, xmax){
  labmin = round(max(xmin, (xmed-(clamp*xmad))),3)
  labmax = round(min(xmax, (xmed+(clamp*xmad))),3)
  labmed = round(xmed, 3)
  labels = c(labmin, labmed, labmax)
  stopoffset =  paste(as.character(seq(0,100,(100/(length(palette)-1)))),"%", sep="")
  stopcolor = col2hex(palette)
  for(indB in 1:length(cXml)){
    gradnode = getNodeSet(cXml[indB][[1]], "//*[@id='gradient']")
    for(indC in 1:length(palette)){
      newstop = newXMLNode("stop", attrs=c("offset"=stopoffset[indC],"stop-color"=stopcolor[indC]))
      gradnode[[1]] = addChildren(gradnode[[1]],newstop)
    }
    for(indC in 1:length(labels)){
      node = getNodeSet(cXml[indB][[1]], "//*[@class='legendLabel']")[[indC]]
      nv = paste("\n",labels[indC],"\n", sep="")
      xmlValue(node) = nv
      if(divData==FALSE & indC==2){
        nv = paste("\n","\n",sep="")
        xmlValue(node) = nv
      }
    }
    if(legend==FALSE){
      node = getNodeSet(cXml[indB][[1]], "//*[@class='legendBar']")[[1]]
      removeAttributes(node, "opacity")
      addAttributes(node, opacity="0")
    }
  }
  return(cXml)
}

#' a function used by cerebroViz() to lower opacity for superior regions when data is not provided
#'
#' for each missing superior region, set opacity to 0.
#' @param cXml list containing the xml object for each SVG.
#' @param supReg regions that encompass other brain regions, specified within cerebroViz().
#' @param lupiter hex gradient indices for current iteration (timePoint) in the loop, specified within cerebroViz().
#' @keywords internal
#' @examples
#' maskRegions(cXml, supReg, lupiter)
#maskRegions
maskRegions <- function(cXml, supReg, lupiter){
  missNames = names(lupiter[is.na(lupiter)])
  opacdown = supReg[supReg%in%missNames]
  if(length(opacdown)>0){
    for(indB in 1:length(opacdown)){
      lobename = opacdown[indB]
      lobenode = getNodeSet(cXml[1][[1]], paste("//*[@id='",lobename,"']",sep=""))
        if(length(lobenode)>0){
          node = lobenode[[1]]
          removeAttributes(node,"fill-opacity")
          addAttributes(node, "fill-opacity"="0")
        }
    }
  }
  if(("STR"%in%missNames) & (sum(c("CAU","PUT") %in% missNames)<2)){
    node = getNodeSet(cXml[2][[1]], "//*[@id='STR']")[[1]]
    removeAttributes(node, "fill-opacity")
    addAttributes(node, "fill-opacity"="0")
  }
  return(cXml)
}
