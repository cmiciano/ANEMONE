# About ANEMONE
Differentially expressed genes are commonly clustered by gene expression to learn more about gene function. However, there is a lack of tools that have the ability to explore genes by their regulators. We developed ANEMONE, an intuitive interactive web tool that comes in two parts. The first tool creates clusters of genes based on their cis-regulatory elements within promoters visualized with heatmaps or principal component analysis plots. The second tool provides an interactive regulatory network of genes based on known transcription factors that regulate them and can be used to explore known interactions.
<br>

ANEMONE is hosted at anemone.salk.edu for online use. It was made using RShiny under the R computing environment. It relies on the following packages:
* shiny
* ggplot2
* gplots
* heatmaply
* shinyHeatmaply
* bsplus
* htmltools
* shinythemes
* GeneOverlap
* visNetwork
* igraph
* RColorBrewer

<br>

Questions and comments can be directed to micianoc@gmail.com or mshokhirev@salk.edu

# Running ANEMONE locally
ANEMONE can be run locally on your own Shiny Server. Instructions to install Shiny Server can be found here.
https://shiny.rstudio.com/articles/shiny-server.html
https://docs.rstudio.com/shiny-server

* Create a folder (for example ANEMONElocal) in the location you will host ANEMONE and transfer the ANEMONE files to said folder
* Your directory should look like /srv/shiny-server/ANEMONElocal and contain the files
  * data (a folder with necessary files for ANEMONE)
  * app.R
  * convertSymbols.R
  * genHeatmap.R
  * subsetGenes.R
* Your local version of ANEMONE should be accessible at http://[local_server_url]:3838/ANEMONElocal/

Alternatively, ANEMONE can be run through GitHub using R Studio.

Once you have the above packages installed on your computer, run the following command in your R Studio console

```
shiny::runGithub('ANEMONE', 'cmiciano')

```

This will deploy a local version of ANEMONE in your web browser