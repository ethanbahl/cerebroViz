## ---- echo=TRUE, fig.show='hold', tidy = TRUE----------------------------
library('cerebroViz')
data(ex1)
ex1

## ---- echo=TRUE----------------------------------------------------------
cerebroViz(ex1)


## ---- echo=TRUE----------------------------------------------------------
library('cerebroViz')

cerebroViz(ex1, timepoint = c(1:3))

## ------------------------------------------------------------------------
library('cerebroViz')

dir.create('custom_directory/')

cerebroViz(ex1, outfile = 'custom_directory/a_custom_filename', timepoint = c(1:3))

list.files('custom_directory/')


## ---- echo=TRUE----------------------------------------------------------
library('cerebroViz')

summary(as.vector(ex1))

cerebroViz(ex1, divergent.data = TRUE, outfile = 'divergent')

## ---- echo=TRUE----------------------------------------------------------
ex1 -> ex1_out

ex1_out["MED",1] <- 55

## ---- echo=TRUE----------------------------------------------------------
cerebroViz(ex1_out, outfile = 'outlier')

## ---- echo=TRUE----------------------------------------------------------
cerebroViz(ex1_out, clamp = 3, outfile = 'outlier_clamped')

## ---- echo=TRUE----------------------------------------------------------
ex1 -> ex1_NA
ex1_NA["PON",1] <- NA

cerebroViz(ex1_NA, cross.hatch = TRUE, outfile = 'cross_hatch')

## ---- echo=TRUE----------------------------------------------------------
library(cerebroViz)

cerebroViz(ex1, regCol = c("coral", "antiquewhite", "cornflowerblue"), outfile = "custom_scale")

## ---- echo=TRUE----------------------------------------------------------
library(cerebroViz)

cerebroViz(ex1, svgCol = c("darkgrey", "white", "lightgrey"), outfile = "custom_background")

## ---- echo=TRUE----------------------------------------------------------
library(cerebroViz)
library(RColorBrewer)

cerebroViz(ex1, regCol = brewer.pal(5, "BrBG"), divergent.data = TRUE, outfile = "brewer")

## ---- echo=TRUE----------------------------------------------------------
head(ex1)

head(cerebroScale(ex1, clamp = NULL, divergent.data = FALSE))

## ---- echo=TRUE----------------------------------------------------------
summary(as.vector(ex1))

summary(as.vector(cerebroScale(ex1, clamp = NULL, divergent.data = FALSE)))

summary(as.vector(cerebroScale(ex1, clamp = 3, divergent.data = TRUE)))

## ---- echo=TRUE, fig.width=7, fig.height=7, fig.show='hold'--------------
cerebroScale(ex1, clamp = 3, divergent.data = TRUE) -> ex1_scaled

heatmap(ex1_scaled, Colv = NA, col = c("#D73027", "#FC8D59", "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4"))

cerebroViz(ex1, clamp = 3, divergent.data = TRUE, outfile = "scaled")

## ---- echo = TRUE, results='hide'----------------------------------------
  data(regionMap)
  regionMap

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(regionMap)

## ---- echo=TRUE----------------------------------------------------------
rownames(ex1)[c(13, 18)]
rownames(ex1)[c(13, 18)] <- c('medulla', 'pons')
cerebroViz(ex1, outfile = "missing_name")

## ---- echo = TRUE--------------------------------------------------------
matrix(c("PON", "MED","pons", "medulla"), ncol=2 ) -> cnm
cnm
cerebroViz(ex1, customNames = cnm, outfile = "custom_name")

