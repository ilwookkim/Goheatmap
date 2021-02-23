#' goheatmap
#'
#' ComplexHeatmap with GO Terms using gprofiler2 package.
#' @param mat numeric matrix of the values.
#' @param k number of groups (cutree).
#' @param n_go number of GO Terms to display.
#' @param sources Term sources from g:profiler - GO Terms, KEGG, Reactome, WikiPathways, Transfac, miRTarBase, Human Protein Atlas, CORUM protein complexes, Human Phenotype Ontology ("GO:MF","GO:CC","GO:BP","KEGG","REAC","WP","TF","MIRNA","HPA","CORUM","HP")
#' @examples
#' goheatmap(mat, k = 3, n_go = 3, sources = "GO:BP")
#' @export
#' @import ComplexHeatmap gprofiler2 dendextend

goheatmap <- function(mat, k = 3, n_go = 3, sources = "GO:BP"){
  ht <- as.dendrogram(hclust(dist(mat)), method = "average")
  ht <- color_labels(ht, k = k)
  ht <- color_branches(ht, k = k)
  labels_ht <- labels_colors(ht)

  for (i in 1:k) {
    clusters <- names(labels_ht[which(factor(labels_ht) == unique(labels_ht)[i])])
    gp <- gost(clusters)
    gp_mod = gp$result[,c("source", "term_id", "term_name", "p_value")]
    row.names(gp_mod) = gp_mod$term_id
    gp_term <- filter(gp_mod, grepl(sources, source, ignore.case=TRUE))
    n_go_st <- gp_term %>% top_n(-n_go)
    if(!dim(n_go_st)[1] == 0){
      assign(paste0("lt",i), unlist(lapply(1:n_go, function(x) paste0(n_go_st$term_id[x]))))
    } else {
      assign(paste0("lt",i), unlist(lapply(1:n_go, function(x) paste0("NA"))))
    }
  }
  # text list for GO Terms
  text_list <- lapply(1:k, function(x) get(paste0("lt",x)))

  ha = rowAnnotation(foo = anno_empty(border = FALSE,
                                      width = max_text_width(unlist(text_list)) + unit(1, "mm")))

  Heatmap(cor_df, name = "cor_df", cluster_rows = ht, row_split = k,
          column_split = k,  right_annotation = ha,border = TRUE , show_row_names = FALSE, show_column_names = FALSE)
  for(i in 1:k) {
    decorate_annotation("foo", slice = i, {
      grid.rect(x = 0, width = unit(1, "mm"), gp = gpar(fill = i, col = NA), just = "left")
      grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(1, "mm"), just = "left")
    })
  }

}
