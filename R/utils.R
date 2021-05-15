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
library(scuttle)
library(scDblFinder)
library(scran)
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

#' Examine doublet
#'
#' @param so seurat object
#' @param clusters name of column indicating cluster ids
#' @param samples name of column indicating sample
#' @param ndim # of dimensions to use for single cell doublet detection
#' @param k_val # of neighbors to use for single cell doublet detection
score_doublets <- function(so,
                           cluster_col = "seurat_clusters",
                           sample_col = "orig.ident",
                           ndim = 20,
                           k_val = 50) {
  sample_df <- split(so@meta.data, so@meta.data[[sample_col]])

  res <- map(sample_df,
             ~{tmp <- subset(so, cells = rownames(.x))

             n_clusters <- length(unique(tmp[[cluster_col]]))
             if(n_clusters < 3) {
               dbl.out <- NULL
             } else {
               dbl.out <- scDblFinder::findDoubletClusters(tmp@assays$RNA@data,
                                                           clusters = tmp@meta.data[[cluster_col]],
                                                           assay.type = "logcounts")
             }

             tmp <- FindVariableFeatures(tmp, nfeatures = 2000)
             dbl_scores <- scDblFinder::computeDoubletDensity(tmp@assays$RNA@data,
                                                              subset.row = VariableFeatures(tmp),
                                                              dims = ndim,
                                                              k = k_val)
             dbl_scores <- data.frame(row.names = colnames(tmp),
                                      doublet_scores = dbl_scores)
             list(cluster = dbl.out,
                  cell = dbl_scores)
             })
  so@misc[["doublets"]] <- res

  doublet_score <- map_dfr(so@misc$doublets, ~.x$cell)
  doublet_score <- doublet_score[colnames(so), ]
  so$doublet_score <- doublet_score
  so
}

#' Identify outliers using MAD
#'
#' @param vec vector of data
#' @param groups vector of groups for each data point, passed to batch arg of
#' scuttle::isOutlier
#' @param min if data is less than min, then do not consider an outlier
#' @param max if data is greater than min, then always consider an outlier
check_outliers <- function(vec, groups = NULL, min = NULL, max = NULL, ...){
  res <- scuttle::isOutlier(vec, batch = groups, ...)
  if(!is.null(max)){
    res[vec > max] <- TRUE
  }
  if(!is.null(min)){
    res[vec < min] <- FALSE
  }
  res
}

add_subcluster_ids <- function(so,
                               parent_values = "coarse_cell_type",
                               child_values = "refined_clusters",
                               str_format = "{parent_value} ({child_value})") {
  p <- so@meta.data[[parent_values]]
  names(p) <- rownames(so@meta.data)
  v <- so@meta.data[[child_values]]

  p_vals <- unique(p)
  res <- lapply(p_vals,
                function(x){
                  p_subset <- p[which(p == x)]
                  child_value <- as.integer(factor(v[which(p == x)]))
                  parent_value = x
                  res <- stringr::str_glue(str_format)
                  names(res) <- names(p_subset)
                  res
                }) %>%
    unlist()
  res[names(p)]
}

plot_montage <- function(sobj,
                         feature_list,
                         ncols = 6,
                         outfile = NULL,
                         plot_fxn = plot_umap,
                         plot_args = list(minimal_theme = TRUE)) {

  all_plts <- list()
  rel_heights <- c()
  for(i in seq_along(feature_list)){

    features <- feature_list[[i]]
    id <- names(feature_list)[i]
    n_plts <- length(features)

    plts <- do.call(plot_fxn, args = c(list(seurat_obj = sobj,
                                            feature = features),
                                       plot_args))

    title_plt <- ggplot() +
      theme_void() +
      labs(tag = str_wrap(id, 20)) +
      theme(plot.tag.position = c(0.75, 0.5),
            plot.tag = element_text(angle = 0,
                                    size = 18,
                                    face = "bold"))

    title_plts <- c(list(title_plt), vector("list", ncols - 1))

    if (is.ggplot(plts)) {
      plts <- list(plts)
    }
    plts <- c(title_plts, plts)

    if ((n_plts %% ncols) != 0) {
      to_add <- ncols - (n_plts %% ncols)
      plt_lst <- c(plts, vector("list", to_add))
    } else {
      plt_lst <- plts
    }


    all_plts <- c(all_plts, plt_lst)
    nrows <- ceiling(length(plt_lst) / ncols)
    sub_plot_rel_heights <- c(0.5, rep(1, nrows - 1))
    rel_heights <- c(rel_heights, sub_plot_rel_heights)
  }


  p <- plot_grid(plotlist = all_plts,
                 nrow = ceiling(length(all_plts) / ncols),
                 ncol = ncols,
                 scale = 0.8,
                 rel_heights = rel_heights)


  if (!is.null(outfile)){
    save_plot(outfile,
              p,
              nrow = ceiling(length(all_plts) / ncols),
              ncol = ncols,
              base_asp = 1,
              base_height = 2,
              limitsize = FALSE)
  }


  p

}

