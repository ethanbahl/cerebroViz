#' A tool to visualize biological data mapped to SVG brain diagrams.
#'
#' 'cerebroViz' is a tool for visualizing spatiotemporal data in the brain within an anatomical context. The user inputs a matrix and the tool creates SVG diagrams with color mapping reflective of the input data. 'cerebroViz' supports 30 brain regions used by BrainSpan, GTEx, Roadmap Epigenomics, and more.
#' @param x matrix containing input data. Row names should reflect the appropriate brain regions. Column names may represent different time points or samples and do not require naming
#' @param outfile desired prefix for the output SVG files
#' @param palette character vector of color values for visualizing brain regions. Accepts color names, hex values, and RGB values
#' @param timePoint a numeric vector of columns in 'x' to visualize
#' @param divData logical indicating if input data is to be visualized as it diverges from the median. Defaults to false
#' @param secondaryPalette character vector of length three specifying colors for the brain background, brain outline/legend, and svg background in that order
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
cerebroViz <- function(x, outfile = "cerebroViz_output", palette = NULL, timePoint = 1, divData = FALSE, secondaryPalette = c("white","black","white"), clamp = NULL, naHatch = FALSE, legend = TRUE, customNames = NULL){
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
  if(length(secondaryPalette)!=3) stop("'secondaryPalette' must have length 3.")
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
  srg = c("BS", "FL", "OL", "PL", "TL", "STR")
  suplog = matrix(!is.na(x[rownames(x) %in% srg,timePoint]),ncol=length(timePoint))
  rownames(suplog) = rownames(x)[rownames(x)%in%srg]
  usrg = rowSums(suplog)
  usrg = names(usrg[usrg>0])
  if(length(usrg)>0){
     warning(paste("The following regions encompass other regions of the brain: ", paste(usrg, collapse=", "),". Subregions are masked in output.", sep=""))
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

  #set regions w/ no data to NA
  #name the rows, join the users data with the NA data, alphabetize
  NAmatrix = matrix(data=NA, nrow=length(regions[which(regions%in%rownames(x)==FALSE)]),ncol=ncol(x))
  rownames(NAmatrix) = regions[which(regions%in%rownames(x)==FALSE)]
  datMat = rbind(x, NAmatrix)
  datMat = as.matrix(datMat[order(rownames(datMat)),])
  x_scaled = cerebroScale(datMat, clamp, divData)
  hexInd = round(x_scaled*200+1)
  f = colorRampPalette(palette)
  hexVec = f(201)

  #loop for each timePoint
  svg_outer = system.file("extdata/svg/brain-outer.svg",package="cerebroViz")
  svg_slice = system.file("extdata/svg/brain-slice.svg",package="cerebroViz")
  for(indA in 1:length(timePoint)){
    if(ncol(x)>1) {
      lupiter = hexInd[,timePoint[indA]]
    }
    else {
      lupiter = hexInd
    }
    xmll = xmlTreeParse(svg_outer, useInternalNodes=TRUE)
    xmls = xmlTreeParse(svg_slice, useInternalNodes=TRUE)
    xmlc =  c(xmll, xmls)
    xmlc = editsecondaryPalette(xmlc, secondaryPalette)
    if(naHatch==TRUE){
        xmlc = editnaHatch(xmlc, lupiter)
    }
    xmlc = editpalette(lupiter, xmlc, hexVec, naHatch)
    xmlc = maskRegions(xmlc, srg, lupiter)
    xmlc = editLegend(xmin, xmed, clamp, xmad, xmax, xmlc, palette, legend, divData)
    xmll = xmlc[1][[1]]
    xmls = xmlc[2][[1]]
    saveXML(xmll, paste(outfile,"_outer_",timePoint[indA],".svg",sep=""))
    saveXML(xmls, paste(outfile,"_slice_",timePoint[indA],".svg",sep=""))
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
  fill_matrix = x

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
    abvmed = x[x>=xmed & x<=(xmed+outlrs) & !is.na(x)]
    belmed = x[x<=xmed & x>=(xmed-outlrs) & !is.na(x)]
    if(length(which(!is.na(x))) %% 2 == 0){ #imputing median if even number of data points
      rsc = rescale(c(xmed,abvmed),c(0.5,1))[-1]
      fill_matrix[x>=xmed & x<=(xmed+outlrs) & !is.na(x)] = rsc
      lsc = rescale(c(xmed,belmed),c(0,0.5))[-1]
      fill_matrix[x<=xmed & x>=(xmed-outlrs) & !is.na(x)] = lsc
      fill_matrix[x<(xmed-outlrs) & !is.na(x)] = 0
      fill_matrix[x>(xmed+outlrs) & !is.na(x)] = 1
      x_scaled = fill_matrix
    }
    if((length(which(!is.na(x)))) %% 2 == 1){
      rsc = rescale(abvmed,c(0.5,1))
      fill_matrix[x>=xmed & x<=(xmed+outlrs) & !is.na(x)] = rsc
      lsc = rescale(belmed,c(0,0.5))
      fill_matrix[x<=xmed & x>=(xmed-outlrs) & !is.na(x)] = lsc
      fill_matrix[x<(xmed-outlrs) & !is.na(x)] = 0
      fill_matrix[x>(xmed+outlrs) & !is.na(x)] = 1
      x_scaled = fill_matrix
    }
  }
  if(divData==FALSE){
    nonoutlrs = x[x>=(xmed-outlrs) & x<=(xmed+outlrs) & !is.na(x)]
    xsc = rescale(nonoutlrs,c(0,1))
    fill_matrix[x>=(xmed-outlrs) & x<=(xmed+outlrs) & !is.na(x)] = xsc
    fill_matrix[x<(xmed-outlrs) & !is.na(x)] = 0
    fill_matrix[x>(xmed+outlrs) & !is.na(x)] = 1
    x_scaled = fill_matrix
  }
  return(x_scaled)
}

#' a function used by cerebroViz() to edit the brain outline and brain background.
#'
#' for each xml, remove fill attributes from 'brainBackground' & 'brainOutline' and replace them with the designated colors.
#' @param xmlc list containing the xml object for each SVG.
#' @param secondaryPalette character vector of length three specifying colors for the brain outline, brain background, and svg background in that order.
#' @keywords internal
#' @examples
#' editsecondaryPalette(xmlc, secondaryPalette)
#editsecondaryPalette
editsecondaryPalette <- function(xmlc, secondaryPalette){
  for(indB in 1:length(xmlc)){
    node = getNodeSet(xmlc[indB][[1]], "//*[@id='brainBackground']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(secondaryPalette[1]))
    node = getNodeSet(xmlc[indB][[1]], "//*[@id='brainOutline']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(secondaryPalette[2]))
    node = getNodeSet(xmlc[indB][[1]], "//*[@id='legendRect']")[[1]]
    removeAttributes(node, "stroke")
    addAttributes(node, stroke=col2hex(secondaryPalette[2]))
    node = getNodeSet(xmlc[indB][[1]], "//*[@id='leglableft']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(secondaryPalette[2]))
    node = getNodeSet(xmlc[indB][[1]], "//*[@id='leglabmid']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(secondaryPalette[2]))
    node = getNodeSet(xmlc[indB][[1]], "//*[@id='leglabright']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(secondaryPalette[2]))
    node = getNodeSet(xmlc[indB][[1]], "//*[@id='svgBackground']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(secondaryPalette[3]))
  }
    return(xmlc)
}

#' a function used by cerebroViz() to add cross-hatching to regions with missing data.
#'
#' for each xml, get style and append cross-hatching pattern.
#' @param xmlc list containing the xml object for each SVG.
#' @param lupiter hex gradient indices for current iteration (timePoint) in the loop, specified within cerebroViz().
#' @keywords internal
#' @examples
#' editnaHatch(xmlc, lupiter)
#editnaHatch
editnaHatch <- function(xmlc, lupiter){
  missNames = names(lupiter[is.na(lupiter)])
  for(indB in 1:length(xmlc)){
    for(indC in 1:length(missNames)){
      node = getNodeSet(xmlc[indB][[1]], paste("//*[@id='",missNames[indC],"']",sep=""))[1]
      if(!is.null(node[[1]])){
        style = xmlGetAttr(node[[1]], "style")
        style = paste("fill:url(#hatch00);",style, sep="")
        addAttributes(node[[1]], style=style)
      }
    }
  }
  return(xmlc)
}

#' a function used by cerebroViz() to map input data to brain regions.
#'
#' for each xml, add the appropriate fill attributes for each region
#' @param lupiter hex gradient indices for current iteration (timePoint) in the loop, specified within cerebroViz().
#' @param xmlc list containing the xml object for each SVG.
#' @param hexVec character vector of length 201 containing the hex value gradient for visualization.
#' @param naHatch logical indicating if regions of missing data should be filled with a cross-hatch pattern to differentiate them from the brain's background.
#' @keywords internal
#' @examples
#' editpalette(lupiter, xmlc, hexVec, naHatch)
#editpalette
editpalette <- function(lupiter, xmlc, hexVec, naHatch){
  fillData = lupiter[which(!is.na(lupiter))]
  for(indB in 1:length(xmlc)){
    for(indC in 1:length(fillData)){
      node = getNodeSet(xmlc[indB][[1]], paste("//*[@id='",names(fillData)[indC],"']",sep=""))[1]
      if(!is.null(node[[1]])){
        removeAttributes(node[[1]], "fill")
        addAttributes(node[[1]], fill=hexVec[fillData[indC]])
      }
    }
  }
  if(naHatch==FALSE){
    missNames = lupiter[is.na(lupiter)]
    for(indB in 1:length(xmlc)){
      for(indC in 1:length(missNames)){
        node = getNodeSet(xmlc[indB][[1]], paste("//*[@id='",names(missNames)[indC],"']",sep=""))[1]
        if(!is.null(node[[1]])){
          removeAttributes(node[[1]], "fill-opacity")
          addAttributes(node[[1]], "fill-opacity"=0)
        }
      }
    }
  }
  return(xmlc)
}

#' a function used by cerebroViz() to edit the legend.
#'
#' for each xml, map the appropriate colors to the legend
#' @param xmin minimum of input data.
#' @param xmed median of input data.
#' @param clamp coefficient to the Median Absolute Deviation. Added and subtracted from the median to identify a range of non-outliers. Values external to this range will 'clamped' to extremes of the non-outlier range.
#' @param xmad Median Absolute Deviation of input data (constant=1).
#' @param xmax maximum of input data.
#' @param xmlc list containing the xml object for each SVG.
#' @param palette character vector of color values to use in the visualization. Accepts color names, hex values, and RGB values. For sequential data, it is recommended to use two colors, or a sequence of colors in a gradient. For divergent data, it is recommended to use three colors with a neutral color in the middle.
#' @param legend logical indicating if the legend bar should be visible.
#' @param divData logical indicating if input data is divergent in nature. Default assumes data is sequential.
#' @keywords internal
#' @examples
#' editLegend(xmin, xmed, clamp, xmad, xmax, xmlc, palette, legend, divData)
#editLegend()
editLegend <- function(xmin, xmed, clamp, xmad, xmax, xmlc, palette, legend, divData){
  labmin = round(max(xmin, (xmed-(clamp*xmad))),3)
  labmax = round(min(xmax, (xmed+(clamp*xmad))),3)
  labmed = round(xmed, 3)
  labels = c(labmin, labmed, labmax)
  stopoffset =  paste(as.character(seq(0,100,(100/(length(palette)-1)))),"%", sep="")
  stopcolor = col2hex(palette)
  for(indB in 1:length(xmlc)){
    gradnode = getNodeSet(xmlc[indB][[1]], "//*[@id='gradient']")
    for(indC in 1:length(palette)){
      newstop = newXMLNode("stop", attrs=c("offset"=stopoffset[indC],"stop-color"=stopcolor[indC]))
      gradnode[[1]] = addChildren(gradnode[[1]],newstop)
    }
    for(indC in 1:length(labels)){
      node = getNodeSet(xmlc[indB][[1]], "//*[@class='legendLabel']")[[indC]]
      nv = paste("\n",labels[indC],"\n", sep="")
      xmlValue(node) = nv
      if(divData==FALSE & indC==2){
        nv = paste("\n","\n",sep="")
        xmlValue(node) = nv
      }
    }
    if(legend==FALSE){
      node = getNodeSet(xmlc[indB][[1]], "//*[@class='legendBar']")[[1]]
      removeAttributes(node, "opacity")
      addAttributes(node, opacity="0")
    }
  }
  return(xmlc)
}

#' a function used by cerebroViz() to lower opacity for superior regions when data is not provided
#'
#' for each missing superior region, set opacity to 0.
#' @param xmlc list containing the xml object for each SVG.
#' @param srg regions that encompass other brain regions, specified within cerebroViz().
#' @param lupiter hex gradient indices for current iteration (timePoint) in the loop, specified within cerebroViz().
#' @keywords internal
#' @examples
#' maskRegions(xmlc, srg, lupiter)
#maskRegions
maskRegions <- function(xmlc, srg, lupiter){
  missNames = names(lupiter[is.na(lupiter)])
  opacdown = srg[srg%in%missNames]
  if(length(opacdown)>0){
    for(indB in 1:length(opacdown)){
      lobename = opacdown[indB]
      lobenode = getNodeSet(xmlc[1][[1]], paste("//*[@id='",lobename,"']",sep=""))
        if(length(lobenode)>0){
          node = lobenode[[1]]
          removeAttributes(node,"fill-opacity")
          addAttributes(node, "fill-opacity"="0")
        }
    }
  }
  if(("STR"%in%missNames) & (sum(c("CAU","PUT") %in% missNames)<2)){
    node = getNodeSet(xmlc[2][[1]], "//*[@id='STR']")[[1]]
    removeAttributes(node, "fill-opacity")
    addAttributes(node, "fill-opacity"="0")
  }
  return(xmlc)
}
