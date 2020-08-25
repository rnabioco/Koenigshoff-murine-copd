library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(Seurat) #v3
library(readxl)
library(here)
library(DropletUtils)
library(scbp)
library(presto)
library(qs)
theme_set(theme_cowplot())


# ----------------------------------------------------

proj_dir <- here()
data_dir <- file.path(proj_dir, "data", "cellranger", "2020-08-21")
doc_dir <- file.path(proj_dir, "docs")

fig_dir <- "figs"
mkrs_dir <- "markers"
tbls_dir <- "tables"
obj_dir <- "objects"
walk(c(fig_dir, mkrs_dir, tbls_dir, obj_dir),
     dir.create, showWarnings = F)

out_dirs <- dplyr::lst(fig_dir, mkrs_dir, tbls_dir, obj_dir)

#' write out top n markers
#'
top_marker_matrix <-  function(mkrs,
                               n = NULL,
                               min_pct = 0.10,
                               pvalue = 0.05,
                               only_pos = TRUE,
                               col_to_keep = "logFC"){
  to_keep <- mkrs %>%
    filter(padj < pvalue,
           pct_in > min_pct) %>%
    group_by(group) %>%
    arrange(padj, logFC, .by_group = TRUE)

  if(!is.null(n)){
    to_keep <- to_keep %>%
      slice(1:n)
  }

  if(only_pos){
    to_keep <- filter(to_keep, logFC > 0)
  }

  features_to_keep <- to_keep %>%
    pull(feature) %>%
    unique()

  res <- filter(mkrs, feature %in% features_to_keep) %>%
    select(feature, group, !!sym(col_to_keep)) %>%
    pivot_wider(names_from = "group",
              values_from = all_of(col_to_keep))

  res
}


get_marker_summaries <- function(so,
                             group_col = "seurat_clusters",
                             outdir = "./",
                             prefix = NULL,
                             min_pct = 0.10,
                             pvalue = 0.05,
                             only_pos = TRUE,
                             tsv_output_dir = NULL,
                             xlsx_output_dir = NULL){

  full_markers <- presto::wilcoxauc(so, group_col)

  mkrs <- filter(full_markers, padj < pvalue,  pct_in > min_pct)
  if(only_pos){
    mkrs <- filter(mkrs, logFC > 0)
  }
  mkrs <- mkrs %>%
    group_by(group) %>%
    arrange(padj,
            desc(logFC),
            .by_group = TRUE)

  if(is.null(prefix)){
    prefix <- ""
  } else {
    prefix <- str_c(prefix, "_")
  }

  if(is.null(tsv_output_dir)){
    tsv_output_dir <- outdir
  }

  if(is.null(xlsx_output_dir)){
    xlsx_output_dir <- outdir
  }

  mkrs %>%
    write_tsv(file.path(tsv_output_dir,
                        str_c(prefix, "cluster_markers.tsv")))

  mkrs %>%
    ungroup() %>%
    split(., .$group) %>%
    write_markers_xlsx(.,
                       file.path(xlsx_output_dir,
                                 str_c(prefix, "cluster_markers.xlsx")))

  dplyr::lst(full_markers, mkrs)
}



get_alra_assay <- function(so, file_name, overwrite = FALSE){

  ## only used to add to cellbrowser
  if(overwrite || !file.exists(file_name)){
    so <- RunALRA(so, setDefaultAssay = FALSE)
    gc()
    alra_assay <- so@assays$alra
    qs::qsave(alra_assay, file_name)
  } else {
    alra_assay <- qs::qread(file_name)
  }
  alra_assay
}


# summarize_clustering(sobj,
#                      cluster_prefix = "RNA_snn.res",
#                      padj_cutoff = 0.05,
#                      logfc_cutoff = 0,
#                      pct_in_cutoff = 10) {
#   k_clustering <- stringr::str_subset(colnames(sobj@meta.data),
#                                       paste0("^", cluster_prefix))
#   names(k_clustering) <- k_clustering
#   de_res <- purrr::imap_dfr(k_clustering,
#                             ~presto::wilcoxauc(sobj, .x) %>%
#                               filter(padj < padj_cutoff,
#                                      logFC > logfc_cutoff,
#                                      pct_in > pct_in_cutoff),
#                             .id = "clustering")
#
#   de_res_summary <- de_res %>%
#     group_by(clustering, group) %>%
#     summarize(n = n()) %>%
#     ungroup() %>%
#     split(., .$clustering) %>%
#     imap(., function(x, y) {
#       colnames(x)[colnames(x) == "group"] <- y
#       colnames(x)[colnames(x) == "n"] <- str_c("n_de_genes_", y)
#       select(x, -clustering)})
#
#   de_res_cell_summary <- imap(de_res_summary, function(x, id) {
#     left_join(get_metadata(sobj), x) %>%
#       select(cell, matches(str_c("n_de_genes_", id))) %>%
#       column_to_rownames("cell")
#   })
#
#   for (i in seq_along(de_res_cell_summary)) {
#     sobj <- AddMetaData(sobj, de_res_cell_summary[[i]])
#   }
#
#   node_df <- de_res %>%
#     group_by(clustering, group) %>%
#     summarize(n = n()) %>%
#     ungroup() %>%
#     mutate(node = str_c(clustering, "C", group)) %>%
#     as.data.frame()
#
#   custom_clustree()
# }
#
#
# custom_clustree <- function(layout, node_df) {
#
#   # Add the statistic column (or whatever you want to show)
#   iord <- match(layout$node, node_df$node)
#   layout$statistic <- node_df$n[iord]
#
#   # Plot the graph (modified from clustree.R)
#   gg <- ggraph(layout)
#
#   # Add the edges
#   gg <- gg + geom_edge_link(arrow = arrow(length = unit(1.5 * 5, "points"),
#                                           ends = "last"),
#                             end_cap = circle(9.5 * 1.5, "points"),
#                             start_cap = circle(9.5 * 1.5, "points"),
#                             aes_(colour = ~count,
#                                  alpha = ~in_prop,
#                                  edge_width = ~is_core)) +
#     scale_edge_width_manual(values = c(1.5, 1.5),
#                             guide = "none") +
#     scale_edge_colour_gradientn(colours = viridis::viridis(256)) +
#     scale_edge_alpha(limits = c(0, 1))
#
#   # Add the node points (replace "statistic" with the column you want to show)
#   gg <- gg + clustree:::add_node_points("statistic", "size", 1, colnames(layout))
#
#   # Add the node text
#   gg <- gg + geom_node_text(aes_(label = ~cluster), size = 3,
#                             colour = "black")
#
#   # Plot theme
#   gg <- gg + scale_size(range = c(4, 15)) +
#     ggraph::theme_graph(base_family = "",
#                         plot_margin = ggplot2::margin(2, 2, 2, 2))
#   gg
# }


get_orthologs <- function(vec,
                          ortholog_table = here("dbases", "mouse_orthologs.tsv")){

  if(!file.exists(ortholog_table)){
    library(biomaRt)
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

    orthologs <- getBM(mart = ensembl,
                       attributes = c("external_gene_name",
                                      "mmusculus_homolog_associated_gene_name"))

    write_tsv(orthologs, ortholog_table)
  }
  orthologs <- read_tsv(ortholog_table,
                        col_types = cols(
                          external_gene_name = col_character(),
                          mmusculus_homolog_associated_gene_name = col_character()
                        ))
  orthologs <- filter(orthologs,
                      !is.na(mmusculus_homolog_associated_gene_name)) %>%
    group_by(external_gene_name) %>%
    mutate(n_orthos_h = n()) %>%
    group_by(mmusculus_homolog_associated_gene_name) %>%
    mutate(n_orthos_m = n()) %>%
    filter(n_orthos_h == 1, n_orthos_m == 1) %>%
    dplyr::select(-starts_with("n_orthos")) %>%
    ungroup()
  orthologs <- mutate(orthologs,
                      ortho_id = str_c(external_gene_name, ":",
                                       mmusculus_homolog_associated_gene_name))

  tibble(external_gene_name = vec) %>%
    left_join(orthologs, by = "external_gene_name") %>%
    na.omit() %>%
    pull(mmusculus_homolog_associated_gene_name)

}
