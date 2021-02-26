# **Goheatmap**
Heatmap with GO Terms

### Installation

The **development** version can be installed from GitHub using:

``` r
devtools::install_github("ilwookkim/GOheatmap")
```

### Usage

``` r
library(GOheatmap)
```

**Load example data**
selected STAD-TCGA RNAseq data (gene set for Transcriptional regulation by TP53; 50 patients from TP53 wildtype and 50 from mutation)
``` r 
mat.file <- system.file("extdata", "mat.Rdata", package="GOheatmap")
load(mat.file)
```

**Pre-treatment**

Remove NA, if necessary.

``` r

mat <- data.frame(na.omit(countdata))
knitr::kable(head(mat[, 1:4], 3), "simple")

# column: samples, row: genes (HGNC symbol)

         TCGA.CD.8536.01   TCGA.HU.A4G3.01   TCGA.BR.8058.01   TCGA.CG.5722.01
------  ----------------  ----------------  ----------------  ----------------
AKT1               11071             11013              7099              3754
AKT2                5643              2964              6347              5408
APAF1               2102              2265              2689               626
```

**Run goheatmap**

Parameters k (number of clustering), n_go (number of terms to display), sources [details here](https://biit.cs.ut.ee/gprofiler/page/apis), cor (TRUE for spearman's correlation coefficient, FALSE for z-score normalization from variance stabilizing transformed matrix (DESeq2)), title (Title of heatmap)

``` r
goheatmap(mat, k = 3, n_go = 3, sources = "GO:BP", cor.s = TRUE, title = "GOheatmap")
goheatmap(mat, k = 3, n_go = 3, sources = "KEGG", cor.s = TRUE, title = "GOheatmap")
goheatmap(mat, k = 3, n_go = 3, sources = "GO:BP", cor.s = FALSE, title = "GOheatmap")
```

<img src="inst/extdata/example_go.bp.png"/>
<img src="inst/extdata/example_kegg.png"/>
<img src="inst/extdata/example_go.bp2.png"/>
