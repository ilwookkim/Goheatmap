#' goheatmap
#'
#' ComplexHeatmap with GO Terms using gprofiler2 package.
#' @param mat numeric matrix of the values.
#' @param anno When
#' @param k number of groups (cutree).
#' @param n_go number of GO Terms to display.
#' @param sources Term sources from g:profiler - GO Terms, KEGG, Reactome, WikiPathways, Transfac, miRTarBase, Human Protein Atlas, CORUM protein complexes, Human Phenotype Ontology ("GO:MF","GO:CC","GO:BP","KEGG","REAC","WP","TF","MIRNA","HPA","CORUM","HP")
#' @param cor.s  TRUE for correlation coefficient matrix. FALSE for z-score normalization of matrix. Default to TRUE.
#' @param title Title of heatmap
#' @examples
#' mat.file <- system.file("extdata", "mat.Rdata", package="GOheatmap")
#' load(mat.file)
#' goheatmap(mat, k = 3, n_go = 3, sources = "GO:BP", cor.s = TRUE)
#' @export
#' @import ComplexHeatmap gprofiler2 dendextend magrittr DESeq2 RColorBrewer
#' @importFrom dplyr filter top_n
#' @importFrom grid grid.text grid.rect gpar unit
#' @importFrom grDevices colorRampPalette
#' @importFrom stats as.dendrogram cor dist hclust

goheatmap <- function(mat, anno= NA, k = 3, n_go = 3, sources = "GO:BP", cor.s = TRUE, title = "Goheatmap"){
  if(cor.s){
    mat <- mat[, !sapply(mat, function(x) { stats::sd(x) == 0} )]
    mat <- cor(t(mat), method = "spearman")
  } else {
    mat <- varianceStabilizingTransformation(as.matrix(round(mat)))
    mat <- t(scale(t(mat)))
  }

  ht <- as.dendrogram(hclust(dist(mat)), method = "average")
  ht <- color_labels(ht, k = k)
  ht <- color_branches(ht, k = k)
  labels_ht <- labels_colors(ht)
  for (i in 1:k) {
    clusters <- names(labels_ht[which(factor(labels_ht) == unique(labels_ht)[i])])
    tryCatch({
      go <- gost(clusters)
      gp_mod = go$result[,c("source", "term_id", "term_name", "p_value")]
      row.names(gp_mod) = gp_mod$term_id
      gp_term <- filter(gp_mod, grepl(sources, source, ignore.case=TRUE))
      n_go_st <- gp_term %>% top_n(-n_go)

    }, error = function(e) print("There is no term"))
    if(!dim(n_go_st)[1] == 0){
      assign(paste0("lt",i), unlist(lapply(1:n_go, function(x) paste0(n_go_st$term_id[x]))))
    } else {
      assign(paste0("lt",i), unlist(lapply(1:n_go, function(x) paste0("NA"))))
    }
  }
  text_list <- lapply(1:k, function(x) get(paste0("lt",x)))
  ha = rowAnnotation(foo = anno_empty(border = FALSE,
                                      width = max_text_width(unlist(text_list)) + unit(4, "mm")))
  if(cor.s){
    draw(Heatmap(mat, name = title,
                 cluster_rows = ht,
                 row_split = k,
                 column_split = k,
                 right_annotation = ha,
                 border = TRUE ,
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 column_title = NULL,
                 row_title = NULL),
         column_title = title,
         column_title_gp = gpar(fontsize = 15, fontface = "bold"))
    for(i in 1:k) {
      decorate_annotation("foo", slice = i, {
        grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(2, "mm"), just = "left")
      })}
  } else {
    ha_row = NULL
    tryCatch({
      if(!is.null(anno) && !cor.s){
        top <- colorRampPalette(brewer.pal(12,"Paired"))(length(levels(factor(anno[,2]))))
        names(top) <- levels(factor(anno[,2]))
        ha_row <- HeatmapAnnotation(df = anno[,2], annotation_label	= colnames(anno)[2],which="column", col = list(condition = top))
      }
    }, error = function(e) print("No information of condition"))


    draw(Heatmap(mat, name = title,
                 cluster_rows = ht,
                 row_split = k,
                 top_annotation = ha_row,
                 right_annotation = ha,
                 border = TRUE ,
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 column_title = NULL,
                 row_title = NULL),
         column_title = title,
         column_title_gp = gpar(fontsize = 15, fontface = "bold"))
    for(i in 1:k) {
      decorate_annotation("foo", slice = i, {
        grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(2, "mm"), just = "left")
      })}

  }

}
