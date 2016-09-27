## ---- echo=TRUE, fig.show='hold', warning=FALSE, tidy=TRUE---------------
library('cerebroViz')
data(cerebroEx)
head(cerebroEx)[,c(1:7)]

## ---- echo = FALSE, message=FALSE, warning=FALSE-------------------------
cerebroViz(cerebroEx)

## ---- echo=TRUE, warning=FALSE-------------------------------------------
cerebroViz(cerebroEx)


## ---- echo=TRUE, warning=FALSE-------------------------------------------
library('cerebroViz')

cerebroViz(cerebroEx, timePoint = c(1, 5, 9))

## ---- warning=FALSE------------------------------------------------------
library('cerebroViz')

dir.create('custom_directory/')

cerebroViz(cerebroEx, filePrefix = 'custom_directory/a_custom_filename', timePoint = c(1, 5, 9))

list.files('custom_directory/')


## ---- echo=TRUE, warning=FALSE-------------------------------------------
library('cerebroViz')

summary(as.vector(cerebroEx))

cerebroViz(cerebroEx, divData = TRUE, timePoint = 11, filePrefix = 'divergent')

## ---- echo=TRUE, warning=FALSE-------------------------------------------
ex1_out <- cerebroEx

ex1_out["CB",11] <- 55

## ---- echo=TRUE, warning=FALSE-------------------------------------------
cerebroViz(ex1_out, timePoint = 11, filePrefix = 'outlier')

## ---- echo=TRUE, warning=FALSE-------------------------------------------
cerebroViz(ex1_out, clamp = 3, timePoint = 11, filePrefix = 'outlier_clamped')

## ---- echo=TRUE, warning=FALSE-------------------------------------------
ex_dif <- cerebroEx[,1:3] - cerebroEx2[,1:3]

cerebroViz(ex_dif, filePrefix = 'ex_dif')

## ---- echo=TRUE, warning=FALSE-------------------------------------------
cerebroEx -> ex1

cerebroViz(ex1, regLabel = TRUE, filePrefix = 'regLabel')

## ---- echo=TRUE, warning=FALSE-------------------------------------------
cerebroEx -> ex1

cerebroViz(ex1, figLabel = TRUE, timePoint = c(5, 9), filePrefix = 'figLabel')

## ---- echo=TRUE, warning=FALSE-------------------------------------------
cerebroEx -> ex1

cerebroViz(ex1, naHatch = TRUE, filePrefix = 'naHatch')

## ---- echo=TRUE, warning=FALSE-------------------------------------------
cerebroEx -> ex1

cerebroViz(ex1, legend = FALSE, filePrefix = 'legend')

## ---- echo=TRUE, warning=FALSE-------------------------------------------
library(cerebroViz)
library(RColorBrewer)

cerebroViz(cerebroEx, palette = rev(brewer.pal(11, "PiYG")), divData = TRUE, timePoint = 11, filePrefix = "palette")

## ---- echo=TRUE, warning=FALSE-------------------------------------------
library(cerebroViz)

cerebroViz(cerebroEx, palette = c("cornflowerblue", "antiquewhite", "coral"), timePoint = 11, filePrefix = "custom_palette")

## ---- echo=TRUE, warning=FALSE-------------------------------------------
library(cerebroViz)

cerebroViz(cerebroEx, secPalette = c("darkgrey", "white", "lightgrey"), timePoint = 11, filePrefix = "secPalette")

## ---- echo=TRUE,warning=FALSE--------------------------------------------
head(cerebroEx)[,c(1:7)]

head(cerebroScale(cerebroEx, clamp = NULL, divData = FALSE))[,c(1:7)]

## ---- echo=TRUE, warning=FALSE-------------------------------------------
summary(as.vector(cerebroEx))

summary(as.vector(cerebroScale(cerebroEx, clamp = NULL, divData = FALSE)))

summary(as.vector(cerebroScale(cerebroEx, clamp = 3, divData = TRUE)))

## ---- echo=TRUE, fig.height=7, fig.show='hold', fig.width=7, warning=FALSE----
cerebroScale(cerebroEx, clamp = NULL, divData = TRUE) -> ex1_scaled

heatmap(ex1_scaled, Colv = NA, scale = "none", col = rev(brewer.pal(11, "RdYlBu")))

cerebroViz(cerebroEx, clamp = NULL, divData = TRUE, filePrefix = "scaled")

## ---- echo = TRUE, results='hide', warning=FALSE-------------------------
  data(regionMap)
  regionMap

## ---- echo=FALSE, results='asis', warning=FALSE--------------------------
knitr::kable(regionMap)

## ---- echo=TRUE, warning=TRUE--------------------------------------------
rownames(cerebroEx)[c(3, 14)]
rownames(cerebroEx)[c(3, 14)] <- c('CBC', 'MD')
cerebroViz(cerebroEx, filePrefix = "missing_name", timePoint = 9)

## ---- echo=TRUE, warning=FALSE-------------------------------------------
matrix(c("CB", "THA", "CBC", "MD"), ncol=2 ) -> cnm
cnm
cerebroViz(cerebroEx, customNames = cnm, filePrefix = "custom_name", timePoint = 9)

## ---- include=FALSE------------------------------------------------------
data("cerebroEx")
cerebroEx

## ---- warning=FALSE------------------------------------------------------
library(cerebroViz)
dir.create("gif_directory")
cerebroViz(cerebroEx2, timePoint = c(1:50), divData = TRUE,  filePrefix = "gif_directory/gif", figLabel = TRUE)


