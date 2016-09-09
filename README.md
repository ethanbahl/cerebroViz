# cerebroViz
cerebroViz is a data mapping tool for visualizing spatiotemporal
    data in the brain. The user inputs a matrix of data and the tool outputs
    publication quality SVG diagrams with color mapping reflective of the input
    data. cerebroViz supports 30 brain regions used by BrainSpan, GTEx, Roadmap
    Epigenomics, and more.

# cerebroViz Installation
cerebroViz can be downloaded directly through the repository on GitHub or within R with the following commands.
```
library(devtools)
install_github("ethanbahl/cerebroViz")
```
The vignette can be accessed in R with the following command.
```
vignette(topic="intro_cerebroViz", package="cerebroViz")
```
Additionally, the vignette is available at http://ethanbahl.github.io/cerebroViz/
