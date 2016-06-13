#' A tool to visualize biological data mapped to SVG brain diagrams.
#'
#' 'cerebroViz' is a data mapping tool for visualizing spatiotemporal data in the brain. The user inputs a matrix of data and the tool outputs publication quality SVG diagrams with color mapping reflective of the input data. 'cerebroViz' supports 30 brain regions used by BrainSpan, GTex, Roadmap Epigneomics, and more.
#' @param x a matrix object containing the data to map. Rownames should reflect the appropriate brain region. Columns may represent different time points or replicates.
#' @param timepoint a numeric vector of columns in 'x' to visualize.
#' @param outfile the desired prefix for the output SVG files.
#' @param regCol a character vector of colors to use in the visualization. Accepts color names, hex values, and RGB values. For sequential data, it is recommended to use two colors. Three colors, with a neutral color in the middle, are suggested for divergent data.
#' @param brainCol a character vector of length two specifying colors for the outline and background of the brain.
#' @param backgroundCol a character vector of length one specifying the color of the SVG background.
#' @param divergent.data logical indicating if input data is divergent in nature. Default assumes data is sequential.
#' @param clamp a numeric vector of length one specifying the value to multiply by MAD to determine outliers that will be 'clamped' down to prevent skewed visualizations.
#' @param cross.hatch logical indicating if regions lacking data should be cross-hatched to differentiate them from the brain's background.
#' @param legend.toggle logical indicating if the legend bar should be visible.
#' @param custom.names logical indicating if custom naming conventions are used for the input data. If TRUE, user should complete the mapping spreadsheet.
#' @keywords cerebroViz
#' @import XML
#' @import gplots
#' @import scales
#' @import grDevices
#' @export
#' @examples
#' x = t(apply(apply(rbind(matrix((sample(c(-400:600),260)/100),nrow=26,ncol=10),matrix(NA,nrow=4,ncol=10)),2,sample),1,sample))
#' rownames(x) = c("A1C", "ACC", "AMY", "ANG", "BS", "CAU", "CB", "DFC", "FCX", "HIP", "HTH", "IPC", "ITC", "M1C", "MED", "MFC", "OCX", "OFC", "PCX", "PIT", "PUT", "PON", "S1C", "SN", "STC", "STR", "TCX", "THA", "V1C", "VFC")
#' cerebroViz(x, regCol=c("blue","grey","red"))
cerebroViz = function(x, timepoint=1, outfile = "cerebroViz_output", regCol = c("blue","grey","red"), brainCol = c("white","black"), backgroundCol="white", divergent.data=FALSE, clamp=NULL, cross.hatch=FALSE, legend.toggle=TRUE, custom.names=FALSE){
  require(XML)
  require(gplots)
  require(scales)
  require(grDevices)

################################################ E R R O R   H A N D L I N G ###
  if(class(x)!="matrix") stop("'x' must be of class 'matrix'")
  if(sum(is.na(rownames(x)))>0) stop("rownames of 'x' must be valid")
  if(length(rownames(x))!=nrow(x)) stop("rownames must be supplied for each row in 'x'")
  if(max(timepoint)>ncol(x)) stop("'timepoint' invalid")
  if(sum(timepoint%%1!=0)) stop("'timepoint' invalid")
  if(length(brainCol)!=2) stop("'brainCol' must have length 2")
  if(length(backgroundCol)!=1) stop("'backgroundCol' must have length 1")

#################################################### R E G I O N   S E T U P ###
  #creating the master regions vector
  regions = c("A1C", "ACC", "AMY", "ANG", "BS", "CAU", "CB", "DFC", "FCX", "HIP", "HTH", "IPC", "ITC", "M1C", "MED", "MFC", "OCX", "OFC", "PCX", "PIT", "PUT", "PON", "S1C", "SN", "STC", "STR", "TCX", "THA", "V1C", "VFC")

  #creating the vector for 'parent' regions (regions that encompass others)
  #warn the user that their data for these regions will overwrite visualization for the subregions
  srg = c("BS", "FCX", "OCX", "PCX", "TCX", "STR")
  usrg = srg[srg%in%rownames(x)]
  if(length(usrg)>0){
     warning(paste("The following brain regions encompass other regions of the brain: ", paste(usrg, collapse=", "),". Subregions will be masked in the output.", sep=""))
  }

  xmed = median(x, na.rm=TRUE)
  xmad = mad(x, constant = 1, na.rm=TRUE)
  xmin = min(x, na.rm=TRUE)
  xmax = max(x, na.rm=TRUE)

###################################################### S P R E A D S H E E T ###
  #spreadsheet - naming conventions
  if(custom.names==TRUE){
    dat = read.table(system.file("extdata/map/suggestedmapping.txt", package="cerebroViz"), fill=TRUE, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings="")
    unknames = rownames(x)[which(rownames(x)%in%regions==FALSE)]
    for(i in 1:length(unknames)){
      unkindex = grep(unknames[i], dat[,6])
      inpindex = which(rownames(x)==unknames[i])
      rownames(x)[inpindex]=dat[unkindex,1]
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

################################################## N O R M A L I Z A T I O N ###
  hexInd = cerebroScale(datMat, clamp, xmed, xmad, divergent.data)

################################################################ H E X V E C ###
  f = colorRampPalette(regCol)
  hexVec = f(201)

############################################################ B I G   L O O P ###
 lobesvg = system.file("extdata/svg/brainlobe.svg",package="cerebroViz")
 sagsvg = system.file("extdata/svg/brainsagittal.svg",package="cerebroViz")

 for(j in 1:length(timepoint)){
   #timepoint selection and filling in cross hatching for missing values
   if(ncol(x)>1) {
     tmp = hexInd[,timepoint[j]]
   }
   else {
     tmp = hexInd
   }

   xmll = xmlTreeParse(lobesvg, useInternalNodes=TRUE)
   xmls = xmlTreeParse(sagsvg, useInternalNodes=TRUE)
   xmlc =  c(xmll, xmls)

############################################################ B R A I N C O L ###
  xmlc = edit.brainCol(xmlc, brainCol)

######################################################## C R O S S H A T C H ###
  if(cross.hatch==TRUE){
    xmlc = edit.crossHatch(xmlc)
  }

################################################################ R E G C O L ###
  xmlc = edit.regCol(tmp, xmlc, hexVec)

############################################################## M A S K R E G ###
  xmlc = edit.maskReg(xmlc, x, usrg)

################################################################ L E G E N D ###
  xmlc = edit.legend(xmin, xmed, clamp, xmad, xmax, xmlc, hexVec, legend.toggle, divergent.data)

############################################################## S A V E X M L ###
  xmll = xmlc[1][[1]]
  xmls = xmlc[2][[1]]
  saveXML(xmll, paste(outfile,"_outer_",timepoint[j],".svg",sep=""))
  saveXML(xmls, paste(outfile,"_slice_",timepoint[j],".svg",sep=""))
  message("Success! Your diagrams have been saved.")
  }
}

#' A data scaling function used by cerebroViz()
#'
#' This function scales the data passed to cerebroViz and translates data points to color values.
#' @param datMat
#' @param clamp
#' @param xmed
#' @param xmad
#' @param divergent.data
#' @keywords scale
#' @export
#' @examples
#' cerebroScale(datMat,clamp,xmed,xmad,divergent.data)
cerebroScale = function(datMat, clamp, xmed, xmad, divergent.data){
  tmp = datMat
  ol = clamp*xmad
  if(divergent.data==TRUE){
    abvmed = tmp[datMat>=xmed & datMat<=(xmed+ol) & !is.na(datMat)]
    belmed = tmp[datMat<=xmed & datMat>=(xmed-ol) & !is.na(datMat)]
    if(length(which(!is.na(tmp))) %% 2 == 0){ #imputing median if even number of data points
      rsc = round(rescale(c(xmed,abvmed),c(101,201)))[-1]
      tmp[datMat>=xmed & datMat<=(xmed+ol) & !is.na(datMat)] = rsc
      lsc = round(rescale(c(xmed,belmed),c(1,101)))[-1]
      tmp[datMat<=xmed & datMat>=(xmed-ol) & !is.na(datMat)] = lsc
      hexInd = tmp
    }
    if((length(which(!is.na(tmp)))) %% 2 == 1){
      rsc = round(rescale(abvmed,c(101,201)))
      tmp[datMat>=xmed & datMat<=(xmed+ol) & !is.na(datMat)] = rsc
      lsc = round(rescale(belmed,c(1,101)))
      tmp[datMat<=xmed & datMat>=(xmed-ol) & !is.na(datMat)] = lsc
      hexInd = tmp
    }
  }
  if(divergent.data==FALSE){
      hexInd = round(rescale(datMat, to=c(1, 201), from=range(datMat, na.rm=TRUE, finite=TRUE)))
  }
  return(hexInd)
}

#' a function used by cerebroViz() to edit the brain outline and brain background.
#'
#' for each xml, remove fill attributes from 'brainBackground' & 'brainOutline' and replace them with the designated colors.
#' @param xmlc
#' @param brainCol
#' @keywords brainCol
#' @examples
#' edit.brainCol(xmlc, brainCol)
#edit.brainCol
edit.brainCol = function(xmlc, brainCol){
  for(k in 1:length(xmlc)){
    node = getNodeSet(xmlc[k][[1]], "//*[@id='brainBackground']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(brainCol[1]))
    node = getNodeSet(xmlc[k][[1]], "//*[@id='brainOutline']")[[1]]
    removeAttributes(node, "fill")
    addAttributes(node, fill=col2hex(brainCol[2]))
  }
    return(xmlc)
}

#' a function used by cerebroViz() to add cross-hatching to regions with missing data.
#'
#' for each xml, get style and append cross-hatching pattern.
#' @param xmlc
#' @keywords cross.hatch
#' @examples
#' edit.crossHatch(xmlc)
#edit.crossHatch
edit.crossHatch = function(xmlc){
  nhatch = names(tmp[which(is.na(tmp))])
  if("STR"%in%nhatch & "CAU"%in%nhatch==FALSE | "PUT"%in%nhatch==FALSE){
    nhatch = nhatch[-(which(nhatch=="STR"))]
  }
  for(k in 1:length(xmlc)){
    for(m in 1:length(nhatch)){
      node = getNodeSet(xmlc[k][[1]], paste("//*[@id='",nhatch[m],"']",sep=""))[1]
      if(is.null(node[[1]])==FALSE){
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
#' @param tmp
#' @param xmlc
#' @param hexVec
#' @keywords regCol
#' @examples
#' edit.regCol(tmp,xmlc,hexVec)
#edit.regCol
edit.regCol = function(tmp, xmlc, hexVec){
  nfill = tmp[which(!is.na(tmp))]
  for(k in 1:length(xmlc)){
    for(m in 1:length(nfill)){
      node = getNodeSet(xmlc[k][[1]], paste("//*[@id='",names(nfill)[m],"']",sep=""))[1]
      if(is.null(node[[1]])==FALSE){
        removeAttributes(node[[1]], "fill")
        addAttributes(node[[1]], fill=hexVec[nfill[m]])
      }
    }
  }
  nhatch = tmp[which(is.na(tmp))]
  for(k in 1:length(xmlc)){
    for(m in 1:length(nhatch)){
      node = getNodeSet(xmlc[k][[1]], paste("//*[@id='",names(nhatch)[m],"']",sep=""))[1]
      if(is.null(node[[1]])==FALSE){
        removeAttributes(node[[1]], "fill-opacity")
        addAttributes(node[[1]], "fill-opacity"=0)
      }
    }
  }
  return(xmlc)
}

#' a function used by cerebroViz() to edit the legend.
#'
#' for each xml, map the appropriate colors to the legend
#' @param xmin
#' @param xmed
#' @param clamp
#' @param xmad
#' @param xmax
#' @param xmlc
#' @param hexVec
#' @param legend.toggle
#' @param divergent.data
#' @keywords legend
#' @examples
#' edit.legend(xmin, xmed, clamp, xmad, xmax, xmlc, hexVec, legend.toggle)
#edit.legend()
edit.legend = function(xmin, xmed, clamp, xmad, xmax, xmlc, hexVec, legend.toggle, divergent.data){
  labmin = round(max(xmin, (xmed-(clamp*xmad))),3)
  labmax = round(min(xmax, (xmed+(clamp*xmad))),3)
  labmed = round(xmed, 3)
  labels = c(labmin, labmed, labmax)
  for(k in 1:length(xmlc)){
    node = getNodeSet(xmlc[k][[1]], "//*[@class='legend']")
    for(m in 1:length(node)){
      removeAttributes(node[[m]],"fill")
      addAttributes(node[[m]],fill=hexVec[m])
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

#' a function used by cerebroViz() to mask subregions when superior regions are supplied in input data.
#'
#' for each superior region, get children region nodes and set opacity to 0.
#' @param xmlc
#' @param x
#' @keywords usrg
#' @examples
#' edit.maskReg(xmlc,x,usrg)
#edit.maskReg
edit.maskReg = function(xmlc, x, usrg){
  if(length(usrg)>0){
    for(m in 1:length(usrg)){
      lobename = usrg[m]
      lobenode = getNodeSet(xmlc[1][[1]], paste("//*[@class='",lobename,"']",sep=""))
        for(lobeind in 1:length(lobenode)){
          if(length(lobenode)>0){
          node = lobenode[[lobeind]]
          removeAttributes(node,"fill-opacity")
          addAttributes(node, "fill-opacity"=0)
        }
      }
    }
  }
  if("STR"%in%rownames(x)){
    node = getNodeSet(xmlc[2][[1]], "//*[@id='STR']")[[1]]
    removeAttributes(node, "fill-opacity")
    addAttributes(node, "fill-opacity"=1)
  }
  return(xmlc)
}
