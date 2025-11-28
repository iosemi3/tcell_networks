#!/usr/bin/env Rscript
## ============================================================================
## 05_rewiring_vs_pseudotime.R
##
## - Derive a rewiring score for backbone TFs:
##     Rewiring_score ~ SS(set_i, set_j) * node_centrality(backbone, subtype_i)
## - Integrate with Monocle3 pseudotime per subtype
## - Plot backbone rewiring along activation/differentiation axis
## ============================================================================

library(igraph)
library(dplyr)
library(reshape2)
library(ggplot2)

## -------------------- Input objects ----------------------------------------

# graph list (same as in previous scripts)
graphs <- readRDS(
  "/Users/sebas/Desktop/tcell network paper submission/rds_files/tbl_graphs.RDS"
)

# backbone TFs (keep identical to 02/04)
backbone_bois <- c("BCL6","CREM","HMGB2","MAF","PBX4","RBPJ","STAT4","ZNF292")

# if you saved all_results_with_backbone_presence from script 04, load it:
# all_results_with_backbone_presence <- readRDS(".../all_results_with_backbone_presence.rds")

## -------------------- SS coefficient + subtype-level gene sets -------------

szymkiewicz_simpson <- function(set1, set2) {
  if (length(set1) == 0 || length(set2) == 0) return(0)
  length(intersect(set1, set2)) / min(length(set1), length(set2))
}

# Collapse communities-with-backbones (from 04) into subtype-level gene sets
# all_results_with_backbone_presence: list per walktrap object → list of groups
extract_and_compare_genes_v2 <- function(all_results) {
  group_genes   <- list()
  group_drivers <- list()
  
  for (result_list in all_results) {
    for (group_name in names(result_list)) {
      simple_name <- gsub("^CLUSTER_WALKTRAP_", "", group_name)
      simple_name <- sub("_Group[0-9]+.*$", "", simple_name)
      simple_name <- gsub("_graph_obj", "", simple_name)
      
      genes_here   <- result_list[[group_name]]$Genes
      drivers_here <- result_list[[group_name]]$Present_Backbone_Drivers
      
      if (!is.null(group_genes[[simple_name]])) {
        group_genes[[simple_name]] <- unique(c(group_genes[[simple_name]], genes_here))
      } else {
        group_genes[[simple_name]] <- genes_here
      }
      
      if (!is.null(group_drivers[[simple_name]])) {
        group_drivers[[simple_name]] <- unique(c(group_drivers[[simple_name]], drivers_here))
      } else {
        group_drivers[[simple_name]] <- drivers_here
      }
    }
  }
  
  subtypes <- names(group_genes)
  
  # build long df of SS values between subtype gene sets
  final_df <- data.frame(
    Subtype    = character(),
    Comparison = character(),
    SS_Value   = numeric(),
    Backbone   = character(),
    stringsAsFactors = FALSE
  )
  
  for (i in subtypes) {
    for (j in subtypes) {
      if (i == j) next
      ss_value <- szymkiewicz_simpson(group_genes[[i]], group_genes[[j]])
      
      if (!is.null(group_drivers[[i]]) && length(group_drivers[[i]]) > 0) {
        for (backbone in group_drivers[[i]]) {
          final_df <- rbind(
            final_df,
            data.frame(
              Subtype    = i,
              Comparison = j,
              SS_Value   = ss_value,
              Backbone   = backbone,
              stringsAsFactors = FALSE
            )
          )
        }
      } else {
        final_df <- rbind(
          final_df,
          data.frame(
            Subtype    = i,
            Comparison = j,
            SS_Value   = ss_value,
            Backbone   = NA_character_,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  final_df
}

overlaps_melteddf <- extract_and_compare_genes_v2(all_results_with_backbone_presence)

## -------------------- Node-level metric per backbone & subtype -------------

# Final version: “actual_edges_to_degree_centrality” → effectively degree-based
actual_edges_to_degree_centrality <- function(graph, node) {
  in_degree  <- degree(graph, node, mode = "in")
  out_degree <- degree(graph, node, mode = "out")
  total_degree <- in_degree + out_degree
  
  neighbors_vec <- unique(c(
    neighbors(graph, node, mode = "out"),
    neighbors(graph, node, mode = "in")
  ))
  
  subgraph <- induced_subgraph(graph, neighbors_vec)
  actual_edges <- gsize(subgraph)  # kept for completeness, not used in ratio
  
  if (total_degree > 0) {
    ratio <- total_degree
  } else {
    ratio <- 0
  }
  ratio
}

compute_edges_to_degree_centrality_matrix <- function(graphs, nodes) {
  graph_names <- gsub("_graph_obj$", "", names(graphs))
  ratio_matrix <- matrix(
    NA_real_,
    nrow = length(graphs),
    ncol = length(nodes),
    dimnames = list(graph_names, nodes)
  )
  
  for (i in seq_along(graphs)) {
    g <- graphs[[i]]
    for (node in nodes) {
      if (node %in% V(g)$name) {
        ratio_matrix[i, node] <- actual_edges_to_degree_centrality(g, node)
      }
    }
  }
  
  ratio_matrix
}

edges_to_degree_centrality_matrix <- compute_edges_to_degree_centrality_matrix(
  graphs,
  backbone_bois
)

## -------------------- Rewiring score = SS * node_metric --------------------

# add rewiring score to overlaps_melteddf
overlaps_melteddf$Rewiring_score <- NA_real_

# ensure rownames of matrix are subtype labels used in overlaps_melteddf
rownames(edges_to_degree_centrality_matrix) <- gsub(
  "_graph_obj$", "",
  rownames(edges_to_degree_centrality_matrix)
)

for (i in seq_len(nrow(overlaps_melteddf))) {
  subtype      <- overlaps_melteddf$Subtype[i]
  backbone_gene <- overlaps_melteddf$Backbone[i]
  ss_value     <- overlaps_melteddf$SS_Value[i]
  
  if (!is.na(backbone_gene) &&
      subtype %in% rownames(edges_to_degree_centrality_matrix) &&
      backbone_gene %in% colnames(edges_to_degree_centrality_matrix)) {
    
    node_metric <- edges_to_degree_centrality_matrix[subtype, backbone_gene]
    overlaps_melteddf$Rewiring_score[i] <- ss_value * node_metric
  }
}

# min–max normalization
valid_idx <- !is.na(overlaps_melteddf$Rewiring_score)
min_rs <- min(overlaps_melteddf$Rewiring_score[valid_idx])
max_rs <- max(overlaps_melteddf$Rewiring_score[valid_idx])

overlaps_melteddf$Norm_Rewiring_score <- NA_real_
overlaps_melteddf$Norm_Rewiring_score[valid_idx] <-
  (overlaps_melteddf$Rewiring_score[valid_idx] - min_rs) /
  (max_rs - min_rs)

## Optional: export
# write.csv(
#   overlaps_melteddf,
#   "/corgi/sebas/tcell_multi/rewiring/overlaps_melteddf_rewiring.csv",
#   row.names = FALSE
# )

## -------------------- Pseudotime per subtype (Monocle3) --------------------

library(SeuratWrappers)
library(monocle3)
library(Seurat)
library(ggridges)

# Load your integrated CD4+ T-cell object for pseudotime (as in your code)
genes_NKT <- readRDS(
  "/Users/sebas/Desktop/tcell network paper submission/rds_files/T_cells_noCC_DICE_ATAC_final_peaks to genes_NKT.RDS"
)

DefaultAssay(genes_NKT) <- "RNA"

integrated <- genes_NKT

integrated <- ScaleData(integrated)
integrated <- FindVariableFeatures(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30, reduction.name = "UMAP")
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated)
Idents(integrated) <- integrated@meta.data$SingleR.labels

cds <- as.cell_data_set(integrated)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# align Monocle clusters with Seurat identities
cds@clusters$UMAP$clusters <- Idents(genes_NKT)
cds@int_colData@listData$reducedDims$UMAP <- genes_NKT@reductions$umap_RNA@cell.embeddings

# order cells along trajectory (root = naive CD4 T cells)
cellInfo <- data.frame(seuratCluster = Idents(genes_NKT))
root_cells <- rownames(cellInfo)[cellInfo$seuratCluster == "T cells, CD4+, naive"]

cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)

# pseudotime per cell
cds$monocle3_pseudotime <- pseudotime(cds)
pseudotime.RNA <- data.frame(colData(cds))

# quick sanity plots (optional)
# plot_cells(cds,
#            color_cells_by = "pseudotime",
#            label_groups_by_cluster = FALSE,
#            label_branch_points = FALSE,
#            label_roots = FALSE,
#            label_leaves = FALSE)

## average pseudotime per subtype (SingleR labels)
pseudotime.cellsRNA <- pseudotime.RNA
pseudotime.cellsRNA$ident <- pseudotime.cellsRNA$SingleR.labels

avg_pseudo_df <- pseudotime.cellsRNA %>%
  dplyr::group_by(ident) %>%
  dplyr::summarise(avg_pseudo = mean(monocle3_pseudotime, na.rm = TRUE), .groups = "drop")

# harmonise ident → Subtype names
avg_pseudo_df$ident <- gsub("^T cells, CD4\\+\\,\\s*", "", avg_pseudo_df$ident)
avg_pseudo_df <- avg_pseudo_df %>%
  mutate(ident = case_when(
    ident == "naive TREG"        ~ "Treg_naive",
    ident == "NKT cells"         ~ "NK",
    ident == "naive, stimulated" ~ "naive_stim",
    ident == "memory TREG"       ~ "Treg_mem",
    TRUE                         ~ ident
  )) %>%
  dplyr::rename(Subtype = ident)

## -------------------- Merge rewiring and pseudotime ------------------------

pseudo_rewscore_df <- merge(avg_pseudo_df, overlaps_melteddf, by = "Subtype")

# optional custom order for plotting
custom_order <- c("naive", "Treg_naive", "naive_stim",
                  "Th2", "Th17", "Treg_mem", "Th1_17", "Th1", "NK")

pseudo_rewscore_df$Subtype <- factor(
  pseudo_rewscore_df$Subtype,
  levels = custom_order
)

# some backbone–subtype combinations might be missing; you can pad zeros if needed
add_missing_zeros <- function(df, backbone, subtype) {
  if (!any(df$Backbone == backbone & df$Subtype == subtype)) {
    zero_row <- data.frame(
      Subtype            = factor(subtype, levels = levels(df$Subtype)),
      avg_pseudo         = NA_real_,
      Comparison         = NA_character_,
      SS_Value           = 0,
      Backbone           = backbone,
      Rewiring_score     = 0,
      Norm_Rewiring_score = 0,
      stringsAsFactors   = FALSE
    )
    df <- rbind(df, zero_row)
  }
  df
}

backbones_missing_naive <- c("STAT4", "MAF", "CREM", "RBPJ", "ZNF292")
for (b in backbones_missing_naive) {
  pseudo_rewscore_df <- add_missing_zeros(pseudo_rewscore_df, b, "naive")
}

## -------------------- Plots: rewiring vs subtype / pseudotime --------------

# per-backbone smooth curves
ggplot(pseudo_rewscore_df,
       aes(x = Subtype, y = Norm_Rewiring_score, group = Backbone, color = Backbone)) +
  stat_smooth(method = "loess", se = FALSE, size = 1.2) +
  theme_minimal() +
  labs(
    title = "Normalized rewiring score by backbone TF",
    x     = "Subtype (ordered)",
    y     = "Normalized rewiring score",
    color = "Backbone"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# facetted version
ggplot(pseudo_rewscore_df,
       aes(x = Subtype, y = Norm_Rewiring_score, group = Backbone, color = Backbone)) +
  stat_smooth(method = "loess", se = FALSE, size = 1.2) +
  facet_wrap(~ Backbone, nrow = 4, ncol = 2) +
  theme_minimal() +
  labs(
    title = "Normalized rewiring score per backbone TF",
    x     = "Subtype (ordered)",
    y     = "Normalized rewiring score"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## -------------------- Average rewiring per subtype/backbone -----------------

avg_rewiring_score <- pseudo_rewscore_df %>%
  group_by(Subtype, Backbone) %>%
  summarise(
    Avg_Norm_Rewiring_score = mean(Norm_Rewiring_score, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(avg_rewiring_score,
       aes(x = Subtype, y = Avg_Norm_Rewiring_score, group = Backbone, color = Backbone)) +
  geom_line(size = 1.2) +
  geom_point() +
  theme_minimal() +
  labs(
    title = "Average normalized rewiring score per backbone TF",
    x     = "Subtype (ordered)",
    y     = "Average normalized rewiring score",
    color = "Backbone"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Optional export
# write.csv(
#   pseudo_rewscore_df,
#   "/corgi/sebas/tcell_multi/rewiring/pseudo_rewscore_df.csv",
#   row.names = FALSE
# )

## ============================================================================
## 03c. Overlap of backbone-bearing (and driver-bearing) communities
##  - For each subtype and backbone TF, collect communities containing that TF
##  - Compute Szymkiewicz–Simpson overlap between these communities
##  - Visualize global overlap distribution and selected alluvial plots
##  - Repeat for driver TFs
## ============================================================================

library(igraph)
library(reshape2)
library(ggplot2)
library(ggalluvial)
library(ggsci)
library(pheatmap)

## ---------------------------------------------------------------------------
## Helper: overlap coefficient (Szymkiewicz–Simpson)
## ---------------------------------------------------------------------------

overlap_coefficient <- function(set1, set2) {
  if (length(set1) == 0 || length(set2) == 0) {
    return(0)
  }
  intersection_size <- length(intersect(set1, set2))
  min_set_size <- min(length(set1), length(set2))
  intersection_size / min_set_size
}

## ---------------------------------------------------------------------------
## A) BACKBONE COMMUNITIES: pick communities containing backbone.bois
##    and build overlap matrix
## ---------------------------------------------------------------------------

paired_elements_modularity1 <- character()  # subtype
paired_elements_modularity2 <- character()  # backbone TF
backbone_communities        <- list()       # named by "<subtype>_<TF>"

for (obj_name in modularity_groups_ls) {
  # e.g. "Th1_graph_obj_wlktrp_out" -> "Th1_graph_obj"
  stp <- gsub("_graph_obj_wlktrp_out.*", "", obj_name)
  subtype <- stp
  
  comm_list <- get(obj_name)  # list of communities (character vectors)
  
  for (grp in comm_list) {
    # only reasonably sized communities, and must contain at least one backbone TF
    if (length(grp) > 10 && any(backbone.bois %in% grp)) {
      backbone.presence <- backbone.bois[backbone.bois %in% grp]
      
      for (boi in backbone.presence) {
        namekey <- paste0(subtype, "_", boi)  # e.g. "Th1_graph_obj_STAT4"
        
        paired_elements_modularity1 <- c(paired_elements_modularity1, subtype)
        paired_elements_modularity2 <- c(paired_elements_modularity2, boi)
        
        # store community genes in a list keyed by "<subtype>_<TF>"
        backbone_communities[[namekey]] <- grp
      }
    }
  }
}

paired_elements_modularity <- data.frame(
  paired_elements_modularity1 = paired_elements_modularity1,
  paired_elements_modularity2 = paired_elements_modularity2,
  stringsAsFactors = FALSE
)

# Build overlap matrix for backbone communities
bb_keys <- names(backbone_communities)
n_bb    <- length(bb_keys)

overlap_backbone <- matrix(
  0,
  nrow = n_bb,
  ncol = n_bb,
  dimnames = list(bb_keys, bb_keys)
)

for (i in seq_len(n_bb)) {
  set1 <- backbone_communities[[bb_keys[i]]]
  for (j in seq_len(n_bb)) {
    if (i == j) next
    set2 <- backbone_communities[[bb_keys[j]]]
    overlap_backbone[i, j] <- overlap_coefficient(set1, set2)
  }
}

overlap_backbone_df <- as.data.frame(overlap_backbone)

## Optional: export full backbone overlap matrix as supplemental CSV
# write.csv(
#   overlap_backbone_df,
#   "/corgi/sebas/tcell_multi/overlap_matrix_backbones/overlap_backbone_communities.csv"
# )

## ---------------------------------------------------------------------------
## Distribution of overlap coefficients (Figure: histogram)
## ---------------------------------------------------------------------------

value_range <- range(overlap_backbone, na.rm = TRUE)
hist_values <- hist(overlap_backbone, breaks = "FD", plot = FALSE)
hist_df     <- data.frame(values = hist_values$mids, counts = hist_values$counts)

npg_palette <- pal_npg()(3)

ggplot(hist_df, aes(x = values, y = counts, fill = values)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors = npg_palette[1:3]) +
  labs(x = "Overlap coefficient", y = "Frequency") +
  theme_minimal() +
  theme(
    axis.text  = element_text(face = "bold", size = 20),
    axis.title = element_text(face = "bold", size = 20),
    panel.grid = element_blank(),
    legend.position = "none"
  )

## ---------------------------------------------------------------------------
## High / low overlap subsets (OPTIONAL)
## Note: high_overlap_backbone / low_overlap_backbone themselves are not
##       used in the alluvial plots; alluvial uses backbone_functions.
## ---------------------------------------------------------------------------

molten_overlap_backbone <- melt(overlap_backbone, varnames = c("Var1", "Var2"),
                                value.name = "value")

high_overlap_backbone <- subset(
  molten_overlap_backbone,
  value > 0.9 & value != 1
)

low_overlap_backbone <- subset(
  molten_overlap_backbone,
  value < 0.1 & value != 0
)

## ---------------------------------------------------------------------------
## Alluvial plots for selected backbones (uses backbone_functions)
## backbone_functions: DF with columns: subtype, termGO, gene, loggrsize
## ---------------------------------------------------------------------------

# High-overlap selection: naive and Treg_naive, excluding HMGB2 / ZNF292
high_overlap_backbone_functions <- backbone_functions[
  backbone_functions$subtype %in% c("naive", "Treg_naive"),
]

high_overlap_backbone_functions <- high_overlap_backbone_functions[
  !high_overlap_backbone_functions$gene %in% c("HMGB2", "ZNF292"),
]

top <- ggplot(
  data = high_overlap_backbone_functions,
  aes(axis1 = subtype, axis2 = termGO, axis3 = gene, y = loggrsize)
) +
  scale_x_discrete(limits = c("Subtype", "Function", "Gene"),
                   expand = c(0.2, 0.02)) +
  geom_alluvium(aes(fill = gene), width = 1/3, linewidth = 1.5) +
  scale_fill_npg() +
  geom_stratum(aes(fill = gene), width = 1/3, linewidth = 1.5) +
  coord_flip() +
  theme_void() +
  theme(legend.position = "none")

top +
  ggrepel::geom_text_repel(
    aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), NA),
        fontface = "bold"),
    stat = "stratum", size = 6, direction = "y", nudge_x = +0.3
  ) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(after_stat(x) == 2, as.character(after_stat(stratum)), NA),
        fontface = "bold"),
    stat = "stratum", size = 6, direction = "y", nudge_x = +0.4
  ) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(after_stat(x) == 3, as.character(after_stat(stratum)), NA),
        fontface = "bold"),
    stat = "stratum", size = 6, direction = "y", nudge_x = +0.3
  )

# Low-overlap functions (complement of above within backbone_functions)
low_overlap_backbone_functions <- backbone_functions[
  !rownames(backbone_functions) %in% rownames(high_overlap_backbone_functions),
]

low <- ggplot(
  data = low_overlap_backbone_functions,
  aes(axis1 = subtype, axis2 = termGO, axis3 = gene, y = loggrsize)
) +
  scale_x_discrete(limits = c("Subtype", "Function", "Gene"),
                   expand = c(0.2, 0.02)) +
  geom_alluvium(aes(fill = gene), width = 1/3, linewidth = 1.5) +
  scale_fill_npg() +
  geom_stratum(aes(fill = gene), width = 1/3, linewidth = 1.5) +
  coord_flip() +
  theme_void() +
  theme(legend.position = "none")

low +
  ggrepel::geom_text_repel(
    aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), NA),
        fontface = "bold"),
    stat = "stratum", size = 6, direction = "y", nudge_x = +0.3
  ) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(after_stat(x) == 2, as.character(after_stat(stratum)), NA),
        fontface = "bold"),
    stat = "stratum", size = 6, direction = "y", nudge_x = +0.4
  ) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(after_stat(x) == 3, as.character(after_stat(stratum)), NA),
        fontface = "bold"),
    stat = "stratum", size = 6, direction = "y", nudge_x = +0.3
  )

## ---------------------------------------------------------------------------
## Export all backbone communities (per subtype) for Cytoscape / Supplemental
## ---------------------------------------------------------------------------

counter <- 1L

for (obj_name in modularity_groups_ls) {
  stp     <- gsub("_graph_obj_wlktrp_out.*", "", obj_name)
  subtype <- stp
  
  comm_list <- get(obj_name)
  for (grp in comm_list) {
    if (length(grp) > 10 && any(backbone.bois %in% grp)) {
      to.write <- grp
      nameout  <- paste0(subtype, "_", counter, ".csv")
      # write.csv(
      #   to.write,
      #   file = paste0(
      #     "/corgi/sebas/tcell_multi/multiome_local_seb_20221207/backbone_communities/",
      #     nameout
      #   ),
      #   row.names = FALSE
      # )
      counter <- counter + 1L
    }
  }
}

## Optional global heatmap of backbone community overlaps
pheatmap::pheatmap(overlap_backbone,
                   main = "All vs All backbone-community overlap")

## Helper: heatmap for single subtype vs all (rows filtered)
create_heatmap <- function(term) {
  subset_rows <- grepl(term, rownames(overlap_backbone_df))
  subset_df   <- overlap_backbone_df[subset_rows, , drop = FALSE]
  pheatmap::pheatmap(
    subset_df,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean"
  )
}

# Example:
# create_heatmap("Th1_17")

## ---------------------------------------------------------------------------
## B) DRIVER COMMUNITIES: same logic, but for driver.bois
## ---------------------------------------------------------------------------

paired_elements_modularity1 <- character()
paired_elements_modularity2 <- character()
driver_communities          <- list()

for (obj_name in modularity_groups_ls) {
  stp     <- gsub("_graph_obj_wlktrp_out.*", "", obj_name)
  subtype <- stp
  
  comm_list <- get(obj_name)
  for (grp in comm_list) {
    if (length(grp) > 10 && any(driver.bois %in% grp)) {
      driver.presence <- driver.bois[driver.bois %in% grp]
      for (boi in driver.presence) {
        namekey <- paste0(subtype, "_", boi)
        paired_elements_modularity1 <- c(paired_elements_modularity1, subtype)
        paired_elements_modularity2 <- c(paired_elements_modularity2, boi)
        driver_communities[[namekey]] <- grp
      }
    }
  }
}

paired_elements_modularity_driver <- data.frame(
  paired_elements_modularity1 = paired_elements_modularity1,
  paired_elements_modularity2 = paired_elements_modularity2,
  stringsAsFactors = FALSE
)

drv_keys <- names(driver_communities)
n_drv    <- length(drv_keys)

overlap_driver <- matrix(
  0,
  nrow = n_drv,
  ncol = n_drv,
  dimnames = list(drv_keys, drv_keys)
)

for (i in seq_len(n_drv)) {
  set1 <- driver_communities[[drv_keys[i]]]
  for (j in seq_len(n_drv)) {
    if (i == j) next
    set2 <- driver_communities[[drv_keys[j]]]
    overlap_driver[i, j] <- overlap_coefficient(set1, set2)
  }
}

overlap_driver_df <- as.data.frame(overlap_driver)

## Optional: export driver overlap
# write.csv(
#   overlap_driver_df,
#   "/corgi/sebas/tcell_multi/overlap_matrix_drivers/overlap_driver_communities.csv"
# )

pheatmap::pheatmap(overlap_driver,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   main = "Driver-community overlap")

## Example: filter out rows/cols containing 'naive' and replot
rows_no_naive <- grep("naive", rownames(overlap_driver), invert = TRUE)
cols_no_naive <- grep("naive", colnames(overlap_driver), invert = TRUE)
filtered_overlap_driver <- overlap_driver[rows_no_naive, cols_no_naive, drop = FALSE]
pheatmap::pheatmap(filtered_overlap_driver)

hist(filtered_overlap_driver)







