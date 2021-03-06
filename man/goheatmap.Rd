% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/goheatmap.R
\name{goheatmap}
\alias{goheatmap}
\title{goheatmap}
\usage{
goheatmap(
  mat,
  anno = NA,
  k = 3,
  n_go = 3,
  sources = "GO:BP",
  cor.s = TRUE,
  title = "Goheatmap"
)
}
\arguments{
\item{mat}{numeric matrix of the values.}

\item{anno}{When}

\item{k}{number of groups (cutree).}

\item{n_go}{number of GO Terms to display.}

\item{sources}{Term sources from g:profiler - GO Terms, KEGG, Reactome, WikiPathways, Transfac, miRTarBase, Human Protein Atlas, CORUM protein complexes, Human Phenotype Ontology ("GO:MF","GO:CC","GO:BP","KEGG","REAC","WP","TF","MIRNA","HPA","CORUM","HP")}

\item{cor.s}{TRUE for correlation coefficient matrix. FALSE for z-score normalization of matrix. Default to TRUE.}

\item{title}{Title of heatmap}
}
\description{
ComplexHeatmap with GO Terms using gprofiler2 package.
}
\examples{
mat.file <- system.file("extdata", "mat.Rdata", package="GOheatmap")
load(mat.file)
goheatmap(mat, k = 3, n_go = 3, sources = "GO:BP", cor.s = TRUE)
}
