# About ANEMONE
Differentially expressed genes are commonly clustered by gene expression to learn more about gene function. However, there is a lack of tools that have the ability to explore genes by their regulators We developed ANEMONE, an intuitive interactive web tool that comes in two parts. The first tool creates clusters of genes based on their cis-regulatory elements within promoters visualized with heatmaps or principal component analysis plots. The second tool provides an interactive regulatory network of genes based on known transcription factors that regulate them can be used to explore known interactions.

<br>
<br>
ANEMONE is hosted at anemone.salk.edu for online use. It was made using RShiny under the R computing environment (Chang et al., 2019). It is dependent on packages include gplots (Gregory et al., 2019), ggplot2 (Wickham, 2016), heatmaply (Galili et al, 2018), and visNetwork (Almende et al., 2019). Questions and comments can be directed to micianoc@gmail.com or mshokhirev@salk.edu

<br>
<br>
<br>

Almende, B et al. (2019) visNetwork: Network Visualization using 'vis.js' Library. R package version 2.0.8. https://CRAN.R-project.org/package=visNetwork

Chang, W. et al. (2019) shiny: Web Application Framework for R. R package version 1.4.0. https://CRAN.R-project.org/package=shiny

Galili, T. et al. (2018) heatmaply: an R package for creating interactive cluster heatmaps for online publishing. Bioinformatics, 34, 1600–1602. 

Warnes, G. et al. (2019) gplots: Various R Programming Tools for Plotting Data. R package version 3.0.1.1. https://CRAN.R-project.org/package=gplots
