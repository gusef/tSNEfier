# tSNEfier
Simple tool that reads a gene expression dataset, performs single sample gene
set variance analysis and then displays the data in a tSNE plot. Groups of samples
in the tSNE can be selected and used for differential expression analysis.

## Install from Github
```{r install, eval=FALSE}
library(devtools)
install_github("gusef/tSNEfier")

library(tSNEfier)
library(d3Toolbox)
tSNEfier()

```
## Running tSNEfier
First we load a R/Bioconductor ExpressionSet or SummarizedExperiment containing the gene expression set we want to analyze and a .gmt file that contains all the gene sets to consider. Then we run GSVA to project the dataset into pathway space:
![example](https://raw.github.com/gusef/d3Toolbox/master/resources/loading.png)

Once that is done we can run tSNE on the second page. Notice that you can not only select the parameters for tSNE, but also the pathway space tSNE should be run in. In this case we use a custom immune geneset that contains roughly 2000 genes. 
![example](https://raw.github.com/gusef/d3Toolbox/master/resources/tSNE.png)


The tSNE plot can then be overlayed by covariates provided in the ExpressionSet, specified pathways and single genes that are present in the selected gene set below. Each point can be hovered over to identify potential outliers.
![example](https://raw.github.com/gusef/d3Toolbox/master/resources/overlay.gif)

The pathway space for tSNE can be changed which leads to a different set of clustering:
![example](https://raw.github.com/gusef/d3Toolbox/master/resources/second_tSNE.png)
![example](https://raw.github.com/gusef/d3Toolbox/master/resources/second_tSNE2.png)


And finally d3Toolbox allows us to select a group of dots. We can use this group to run differential expression against the rest of the samples and do the same type of analysis in the pathway space.
![example](https://raw.github.com/gusef/d3Toolbox/master/resources/diff_genes)







