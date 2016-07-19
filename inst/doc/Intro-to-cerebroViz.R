## ---- echo=TRUE, fig.show='hold', tidy = TRUE----------------------------
library('cerebroViz')
data(cerebro_example)
head(cerebro_example)

## ---- echo=TRUE----------------------------------------------------------
cerebroViz(cerebro_example)


## ---- echo=TRUE----------------------------------------------------------
library('cerebroViz')

cerebroViz(cerebro_example, timePoint = c(1:3))

## ------------------------------------------------------------------------
library('cerebroViz')

dir.create('custom_directory/')

cerebroViz(cerebro_example, outfile = 'custom_directory/a_custom_filename', timePoint = c(1:3))

list.files('custom_directory/')


## ---- echo=TRUE----------------------------------------------------------
library('cerebroViz')

summary(as.vector(cerebro_example))

cerebroViz(cerebro_example, divData = TRUE, outfile = 'divergent')

## ---- echo=TRUE----------------------------------------------------------
cerebro_example -> ex1_out

ex1_out["MED",1] <- 55

## ---- echo=TRUE----------------------------------------------------------
cerebroViz(ex1_out, outfile = 'outlier')

## ---- echo=TRUE----------------------------------------------------------
cerebroViz(ex1_out, clamp = 3, outfile = 'outlier_clamped')

## ---- echo=TRUE----------------------------------------------------------
cerebro_example -> ex1_NA
ex1_NA["PON",1] <- NA

cerebroViz(ex1_NA, naHatch = TRUE, outfile = 'cross_hatch')

## ---- echo=TRUE----------------------------------------------------------
library(cerebroViz)
library(RColorBrewer)

cerebroViz(cerebro_example, palette = brewer.pal(5, "BrBG"), divData = TRUE, outfile = "brewer")

## ---- echo=TRUE----------------------------------------------------------
library(cerebroViz)

cerebroViz(cerebro_example, palette = c("coral", "antiquewhite", "cornflowerblue"), outfile = "custom_scale")

## ---- echo=TRUE----------------------------------------------------------
library(cerebroViz)

cerebroViz(cerebro_example, secondaryPalette = c("darkgrey", "white", "lightgrey"), outfile = "custom_background")

## ---- echo=TRUE----------------------------------------------------------
head(cerebro_example)

head(cerebroScale(cerebro_example, clamp = NULL, divData = FALSE))

## ---- echo=TRUE----------------------------------------------------------
summary(as.vector(cerebro_example))

summary(as.vector(cerebroScale(cerebro_example, clamp = NULL, divData = FALSE)))

summary(as.vector(cerebroScale(cerebro_example, clamp = 3, divData = TRUE)))

## ---- echo=TRUE, fig.width=7, fig.height=7, fig.show='hold'--------------
cerebroScale(cerebro_example, clamp = NULL, divData = TRUE) -> ex1_scaled

heatmap(ex1_scaled, Colv = NA, scale = "none", col = brewer.pal(7, "RdYlBu"))

cerebroViz(cerebro_example, clamp = NULL, divData = TRUE, outfile = "scaled")

## ---- echo = TRUE, results='hide'----------------------------------------
  data(regionMap)
  regionMap

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(regionMap)

## ---- echo=TRUE----------------------------------------------------------
rownames(cerebro_example)[c(13, 18)]
rownames(cerebro_example)[c(13, 18)] <- c('medulla', 'pons')
cerebroViz(cerebro_example, outfile = "missing_name")

## ---- echo = TRUE--------------------------------------------------------
matrix(c("PON", "MED","pons", "medulla"), ncol=2 ) -> cnm
cnm
cerebroViz(cerebro_example, customNames = cnm, outfile = "custom_name")

