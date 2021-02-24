#' goheatmap
#'
#' ComplexHeatmap with GO Terms using gprofiler2 package.
#' @param mat numeric matrix of the values.
#' @param k number of groups (cutree).
#' @param n_go number of GO Terms to display.
#' @param sources Term sources from g:profiler - GO Terms, KEGG, Reactome, WikiPathways, Transfac, miRTarBase, Human Protein Atlas, CORUM protein complexes, Human Phenotype Ontology ("GO:MF","GO:CC","GO:BP","KEGG","REAC","WP","TF","MIRNA","HPA","CORUM","HP")
#' @param cor  Use correlation coefficient matrix. Default to TRUE.
#' @param title Title of heatmap
#' @examples
#' goheatmap(mat, k = 3, n_go = 3, sources = "GO:BP")
#' @export
#' @import ComplexHeatmap gprofiler2 dendextend magrittr ellipse
#' @importFrom dplyr filter top_n
#' @importFrom grid grid.text grid.rect gpar

goheatmap <- function(mat, k = 3, n_go = 3, sources = "GO:BP", cor = TRUE, title = "Goheatmap"){
  if(cor){
    mat <- mat[, !sapply(mat, function(x) { stats::sd(x) == 0} )]
    mat <- cor(t(mat), method = "spearman")
  } else {
    mat <- data.frame(t(mat))
  }

  ht <- as.dendrogram(hclust(dist(mat)), method = "average")
  ht <- color_labels(ht, k = k)
  ht <- color_branches(ht, k = k)
  labels_ht <- labels_colors(ht)
  for (i in 1:k) {
    clusters <- names(labels_ht[which(factor(labels_ht) == unique(labels_ht)[i])])
    go <- gost(clusters)
    if(is.null(go)){
      print("No information of genes. Check the matirx")
    }
    gp_mod = go$result[,c("source", "term_id", "term_name", "p_value")]
    row.names(gp_mod) = gp_mod$term_id
    gp_term <- filter(gp_mod, grepl(sources, source, ignore.case=TRUE))
    n_go_st <- gp_term %>% top_n(-n_go)
    if(!dim(n_go_st)[1] == 0){
      assign(paste0("lt",i), unlist(lapply(1:n_go, function(x) paste0(n_go_st$term_id[x]))))
    } else {
      assign(paste0("lt",i), unlist(lapply(1:n_go, function(x) paste0("NA"))))
    }
  }
  text_list <- lapply(1:k, function(x) get(paste0("lt",x)))
  ha = rowAnnotation(foo = anno_empty(border = FALSE,
                                      width = max_text_width(unlist(text_list)) + unit(4, "mm")))
  draw(Heatmap(mat, name = title,
               cluster_rows = ht,
               row_split = k,
               column_split = k,
               right_annotation = ha,
               border = TRUE ,
               show_row_names = FALSE,
               show_column_names = FALSE,
               column_title = NULL,
               row_title = NULL))
  for(i in 1:k) {
    decorate_annotation("foo", slice = i, {
      grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(2, "mm"), just = "left")
    })}
}
