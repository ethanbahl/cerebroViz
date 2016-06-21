#' A tool to visualize biological data mapped to SVG brain diagrams.
#'
#' 'cerebroViz' is a data mapping tool for visualizing spatiotemporal data in the brain. The user inputs a matrix of data and the tool outputs publication quality SVG diagrams with color mapping reflective of the input data. 'cerebroViz' supports 30 brain regions used by BrainSpan, GTex, Roadmap Epigneomics, and more.
#' @param x a matrix object containing the data to map. Rownames should reflect the appropriate brain region. Columns may represent different time points or replicates.
#' @param timepoint a numeric vector of columns in 'x' to visualize.
#' @param outfile the desired prefix for the output SVG files.
#' @param regCol character vector of color values to use in the visualization. Accepts color names, hex values, and RGB values. For sequential data, it is recommended to use two colors, or a sequence of colors in a gradient. For divergent data, it is recommended to use three colors with a neutral color in the middle.
#' @param svgCol character vector of length three specifying colors for the brain outline, brain background, and svg background in that order.
#' @param divergent.data logical indicating if input data is divergent in nature. Default assumes data is sequential.
#' @param clamp coefficient to the Median Absolute Deviation. Added and subtracted from the median to identify a range of non-outliers. Values external to this range will 'clamped' to extremes of the non-outlier range.
#' @param cross.hatch logical indicating if regions of missing data should be filled with a cross-hatch pattern to differentiate them from the brain's background.
#' @param legend.toggle logical indicating if the legend bar should be visible.
#' @param customNames dataframe with 2 columns. The first column for cerebroViz convention names, the second column for custom user names.
#' @keywords cerebroViz
#' @import XML
#' @import gplots
#' @import scales
#' @import grDevices
#' @export
#' @examples
#' x = t(apply(apply(rbind(matrix((sample(c(-400:600),260)/100),nrow=26,ncol=10),matrix(NA,nrow=4,ncol=10)),2,sample),1,sample))
#' rownames(x) = c("A1C", "CNG", "AMY", "ANG", "BS", "CAU", "CB", "DFC", "FCX", "HIP", "HTH", "IPC", "ITC", "M1C", "MED", "MFC", "OCX", "OFC", "PCX", "PIT", "PUT", "PON", "S1C", "SN", "STC", "STR", "TCX", "THA", "V1C", "VFC")
#' cerebroViz(x)
cerebroViz = function(x, timepoint=1, outfile = "cerebroViz_output", regCol = c("#FEECE4","#ee2d26"), svgCol = c("white","black","white"), divergent.data=FALSE, clamp=NULL, cross.hatch=FALSE, legend.toggle=TRUE, customNames=NULL){
  require(XML)
  require(gplots)
  require(scales)
  require(grDevices)

  #creating the master regions vector
  regions = c("A1C", "CNG", "AMY", "ANG", "BS", "CAU", "CB", "DFC", "FCX", "HIP", "HTH", "IPC", "ITC", "M1C", "MED", "MFC", "OCX", "OFC", "PCX", "PIT", "PUT", "PON", "S1C", "SN", "STC", "STR", "TCX", "THA", "V1C", "VFC")
################################################ E R R O R   H A N D L I N G ###
  if(class(x)!="matrix") stop("'x' must be of class 'matrix'")
  if(sum(is.na(rownames(x)))>0) stop("rownames of 'x' must be valid")
  if(length(rownames(x))!=nrow(x)) stop("rownames must be supplied for each row in 'x'")
  if(max(timepoint)>ncol(x)) stop("'timepoint' invalid")
  if(sum(timepoint%%1!=0)) stop("'timepoint' invalid")
  if(length(svgCol)!=3) stop("'svgCol' must have length 3")
  if(length(regCol)==3 & divergent.data==FALSE | length(regCol)==2 & divergent.data==TRUE){
    warning("recommended usage: 2 colors (regCol) for sequential data and 3 colors for divergent data.")
  }
  if(!is.null(customNames)){
    if(ncol(customNames)!=2) stop("unexpected input for customNames")
  }
  if(is.null(customNames) & sum(rownames(x)%in%regions==FALSE)>0) warning(paste("Unknown rownames in input data: ",paste(rownames(x)[rownames(x)%in%regions==FALSE],collapse=", "),". Unknown regions will be excluded from visualization. See the help page for 'customNames' argument.",sep=""))

#################################################### R E G I O N   S E T U P ###
  #creating the vector for 'parent' regions (regions that encompass others) and a warning of overshadowing.
  srg = c("BS", "FCX", "OCX", "PCX", "TCX", "STR")
  suplog = !is.na(x[rownames(x) %in% srg,timepoint])
  if(length(timepoint)==1){
    usrg = names(suplog[suplog==TRUE])
  }
  if(length(timepoint)>1){
    usrg = rowSums(suplog)
    usrg = names(usrg[usrg>0])
  }
  if(length(usrg)>0){
     warning(paste("The following regions encompass other regions of the brain: ", paste(usrg, collapse=", "),". Subregions are masked in output.", sep=""))
  }

  xmed = median(x, na.rm=TRUE)
  xmad = mad(x, constant = 1, na.rm=TRUE)
  xmin = min(x, na.rm=TRUE)
  xmax = max(x, na.rm=TRUE)

###################################################### C U S T O M N A M E S ###
  #customNames
  if(!is.null(customNames)){
    if(sum(customNames[,2]%in%regions)>0) stop(paste("customNames contains region names already used in cerebroViz convention: ",paste(customNames[customNames[,2]%in%regions,2],collapse=", "),sep=""))
    for(i in 1:nrow(customNames)){
      rownames(x)[rownames(x)==customNames[i,2]]=customNames[i,1]
    }
  }

################################################################## C L A M P ###
  #set the default clamp value (no clamping)
  avoidClamp = max(abs(xmed-xmin),abs(xmed-xmax))/xmad
  if(is.null(clamp)){
    clamp = avoidClamp+.01
  }
  if(clamp<=0) stop("clamp must be >0")
  pctOL = round(length(which(x<=(xmed-(clamp*xmad)) | x>=(xmed+(clamp*xmad))))/length(x)*100,2)
  if(pctOL>0){
    warning(paste("The clamp value of ", clamp," will clamp ",pctOL,"% of input values (outliers) to the minimum/maximum colors. Minimum and maximum values displayed on figure legends represent the values outliers are clamped to.", sep=""))
  }

################################################################ D A T M A T ###
  #set regions w/ no data to NA
  #name the rows, join the users data with the NA data, alphabetize
  NAmatrix = matrix(data=NA, nrow=length(regions[which(regions%in%rownames(x)==FALSE)]),ncol=ncol(x))
  rownames(NAmatrix) = regions[which(regions%in%rownames(x)==FALSE)]
  datMat = rbind(x, NAmatrix)
  datMat = as.matrix(datMat[order(rownames(datMat)),])
  x_scaled = cerebroScale(datMat, clamp, divergent.data)
  hexInd = round(x_scaled*200+1)
  f = colorRampPalette(regCol)
  hexVec = f(201)

############################################################ B I G   L O O P ###
  lobesvg = system.file("extdata/svg/brainlobe.svg",package="cerebroViz")
  sagsvg = system.file("extdata/svg/brainsagittal.svg",package="cerebroViz")
  for(j in 1:length(timepoint)){
    if(ncol(x)>1) {
      tmp = hexInd[,timepoint[j]]
    }
    else {
      tmp = hexInd
    }
    xmll = xmlTreeParse(lobesvg, useInternalNodes=TRUE)
    xmls = xmlTreeParse(sagsvg, useInternalNodes=TRUE)
    xmlc =  c(xmll, xmls)
    xmlc = editSvgCol(xmlc, svgCol)
    if(cross.hatch==TRUE){
        xmlc = editCrossHatch(xmlc, tmp)
    }
    xmlc = editRegCol(tmp, xmlc, hexVec, cross.hatch)
    xmlc = unmaskRegions(xmlc, srg, tmp)
    xmlc = editLegend(xmin, xmed, clamp, xmad, xmax, xmlc, regCol, legend.toggle, divergent.data)
    xmll = xmlc[1][[1]]
    xmls = xmlc[2][[1]]
    saveXML(xmll, paste(outfile,"_outer_",timepoint[j],".svg",sep=""))
    saveXML(xmls, paste(outfile,"_slice_",timepoint[j],".svg",sep=""))
    message("Success! Your diagrams have been saved.")
    }
}

#' A function to scale sequential and divergent data to a 0:1 range
#'
#' This function scales input data to a 0:1 range. If divergent.data=FALSE it will scale linearly. If divergent.data=TRUE, it will utilize a polylinear scale centered at the median of the input data.
#' @param x input data matrix
#' @param clamp coefficient to the Median Absolute Deviation. Added and subtracted from the median to identify a range of non-outliers. Values external to this range will 'clamped' to extremes of the non-outlier range.
#' @param divergent.data logical indicating if input data is divergent in nature. Default assumes data is sequential.
#' @keywords cerebroScale
#' @export
#' @examples
#' cerebroScale(x, clamp = 100, divergent.data=FALSE)
#cerebroScale
cerebroScale = function(x, clamp, divergent.data){
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
  pctOL = round(length(which(x<=(xmed-(outlrs)) | x>=(xmed+(outlrs))))/length(x)*100,2)
  if(pctOL>0){
    warning(paste("The clamp value of ", clamp," will clamp ",pctOL,"% of input values (outliers) to 0 or 1.", sep=""))
  }

  if(divergent.data==TRUE){
    abvmed = x[x>=xmed & x<=(xmed+outlrs) & !is.na(x)]
    belmed = x[x<=xmed & x>=(xmed-outlrs) & !is.na(x)]
    if(length(which(!is.na(x))) %% 2 == 0){ #imputing median if even number of data points
      rsc = rescale(c(xmed,abvmed),c(0.5,1))[-1]
      fill_matrix[x>=xmed & x<=(xmed+outlrs) & !is.na(x)] = rsc
      lsc = rescale(c(xmed,belmed),c(0,0.5))[-1]
      fill_matrix[x<=xmed & x>=(xmed-outlrs) & !is.na(x)] = lsc
      x_scaled = fill_matrix
    }
    if((length(which(!is.na(x)))) %% 2 == 1){
      rsc = rescale(abvmed,c(0.5,1))
      fill_matrix[x>=xmed & x<=(xmed+outlrs) & !is.na(x)] = rsc
      lsc = rescale(belmed,c(0,0.5))
      fill_matrix[x<=xmed & x>=(xmed-outlrs) & !is.na(x)] = lsc
      x_scaled = fill_matrix
    }
  }
  if(divergent.data==FALSE){
      x_scaled = rescale(x, to=c(0,1), from=range(x, na.rm=TRUE, finite=TRUE))
  }
  return(x_scaled)
}

#' a function used by cerebroViz() to edit the brain outline and brain background.
#'
#' for each xml, remove fill attributes from 'brainBackground' & 'brainOutline' and replace them with the designated colors.
#' @param xmlc list containing the xml object for each SVG.
#' @param svgCol character vector of length three specifying colors for the brain outline, brain background, and svg background in that order.
#' @keywords internal
#' @examples
#' editSvgCol(xmlc, svgCol)
#editSvgCol
editSvgCol = function(xmlc, svgCol){
  for(k in 1:length(xmlc)){
    node = getNodeSet(xmlc[k][[1]], "//*[@id='brainBackground']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(svgCol[1]))
    node = getNodeSet(xmlc[k][[1]], "//*[@id='brainOutline']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(svgCol[2]))
    node = getNodeSet(xmlc[k][[1]], "//*[@id='svgBackground']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(svgCol[3]))
  }
    return(xmlc)
}

#' a function used by cerebroViz() to add cross-hatching to regions with missing data.
#'
#' for each xml, get style and append cross-hatching pattern.
#' @param xmlc list containing the xml object for each SVG.
#' @param tmp hex gradient indices for current iteration (timepoint) in the loop, specified within cerebroViz().
#' @keywords internal
#' @examples
#' editCrossHatch(xmlc, tmp)
#editCrossHatch
editCrossHatch = function(xmlc, tmp){
  nhatch = names(tmp[is.na(tmp)])
  for(k in 1:length(xmlc)){
    for(m in 1:length(nhatch)){
      node = getNodeSet(xmlc[k][[1]], paste("//*[@id='",nhatch[m],"']",sep=""))[1]
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
#' @param tmp hex gradient indices for current iteration (timepoint) in the loop, specified within cerebroViz().
#' @param xmlc list containing the xml object for each SVG.
#' @param hexVec character vector of length 201 containing the hex value gradient for visualization.
#' @param cross.hatch logical indicating if regions of missing data should be filled with a cross-hatch pattern to differentiate them from the brain's background.
#' @keywords internal
#' @examples
#' editRegCol(tmp, xmlc, hexVec, cross.hatch)
#editRegCol
editRegCol = function(tmp, xmlc, hexVec, cross.hatch){
  nfill = tmp[which(!is.na(tmp))]
  for(k in 1:length(xmlc)){
    for(m in 1:length(nfill)){
      node = getNodeSet(xmlc[k][[1]], paste("//*[@id='",names(nfill)[m],"']",sep=""))[1]
      if(!is.null(node[[1]])){
        removeAttributes(node[[1]], "fill")
        addAttributes(node[[1]], fill=hexVec[nfill[m]])
      }
    }
  }
  if(cross.hatch==FALSE){
    nhatch = tmp[is.na(tmp)]
    for(k in 1:length(xmlc)){
      for(m in 1:length(nhatch)){
        node = getNodeSet(xmlc[k][[1]], paste("//*[@id='",names(nhatch)[m],"']",sep=""))[1]
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
#' @param regCol character vector of color values to use in the visualization. Accepts color names, hex values, and RGB values. For sequential data, it is recommended to use two colors, or a sequence of colors in a gradient. For divergent data, it is recommended to use three colors with a neutral color in the middle.
#' @param legend.toggle logical indicating if the legend bar should be visible.
#' @param divergent.data logical indicating if input data is divergent in nature. Default assumes data is sequential.
#' @keywords internal
#' @examples
#' editLegend(xmin, xmed, clamp, xmad, xmax, xmlc, regCol, legend.toggle, divergent.data)
#editLegend()
editLegend = function(xmin, xmed, clamp, xmad, xmax, xmlc, regCol, legend.toggle, divergent.data){
  labmin = round(max(xmin, (xmed-(clamp*xmad))),3)
  labmax = round(min(xmax, (xmed+(clamp*xmad))),3)
  labmed = round(xmed, 3)
  labels = c(labmin, labmed, labmax)
  stopoffset =  paste(as.character(seq(0,100,(100/(length(regCol)-1)))),"%", sep="")
  stopcolor = col2hex(regCol)
  for(k in 1:length(xmlc)){
    gradnode = getNodeSet(xmlc[k][[1]], "//*[@id='gradient']")
    for(i in 1:length(regCol)){
      newstop = newXMLNode("stop", attrs=c("offset"=stopoffset[i],"stop-color"=stopcolor[i]))
      gradnode[[1]] = addChildren(gradnode[[1]],newstop)
    }
    for(m in 1:length(labels)){
      node = getNodeSet(xmlc[k][[1]], "//*[@class='legendLabel']")[[m]]
      nv = paste("\n",labels[m],"\n", sep="")
      xmlValue(node) = nv
      if(divergent.data==FALSE & m==2){
        nv = paste("\n","\n",sep="")
        xmlValue(node) = nv
      }
    }
    if(legend.toggle==FALSE){
      node = getNodeSet(xmlc[k][[1]], "//*[@class='legendBar']")[[1]]
      removeAttributes(node, "opacity")
      addAttributes(node, opacity="0")
    }
  }
  return(xmlc)
}

#' a function used by cerebroViz() to unmask subregions when superior regions are not supplied in input data.
#'
#' for each missing superior region, set opacity to 0.
#' @param xmlc list containing the xml object for each SVG.
#' @param srg regions that encompass other brain regions, specified within cerebroViz().
#' @param tmp hex gradient indices for current iteration (timepoint) in the loop, specified within cerebroViz().
#' @keywords internal
#' @examples
#' unmaskRegions(xmlc, srg, tmp)
#unmaskRegions
unmaskRegions = function(xmlc, srg, tmp){
  nhatch = names(tmp[is.na(tmp)])
  opacdown = srg[srg%in%nhatch]
  if(length(opacdown)>0){
    for(m in 1:length(opacdown)){
      lobename = opacdown[m]
      lobenode = getNodeSet(xmlc[1][[1]], paste("//*[@id='",lobename,"']",sep=""))
        if(length(lobenode)>0){
          node = lobenode[[1]]
          removeAttributes(node,"fill-opacity")
          addAttributes(node, "fill-opacity"="0")
        }
    }
  }
  if(("STR"%in%nhatch) & (sum(c("CAU","PUT") %in% nhatch)>0)){
    node = getNodeSet(xmlc[2][[1]], "//*[@id='STR']")[[1]]
    removeAttributes(node, "fill-opacity")
    addAttributes(node, "fill-opacity"="0")
  }
  return(xmlc)
}
