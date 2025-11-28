#!/usr/bin/env Rscript
## ============================================================================
## 04_modularity_communities_GO_enrichment.R
##
## Standalone script:
##   - reads graph list (Pando/xgboost networks) from RDS
##   - runs community detection (walktrap + optional infomap / labelprop)
##   - computes GO enrichment per community (GOBP / GOMF)
##   - annotates communities that contain backbone TFs or driver TFs
##   - provides helper to extract info per gene (for sanity checks)
## ============================================================================

library(igraph)
library(gprofiler2)
library(reshape2)
library(ggplot2)
library(dplyr)

## -------------------- Load graphs -------------------------------------------

graphs <- readRDS(
  "/Users/sebas/Desktop/tcell network paper submission/rds_files/tbl_graphs.RDS"
)
# graphs: named list of igraph/tbl_graph networks
# names(graphs) typically like "<subtype>_graph_obj"

## -------------------- Generic community-detection wrapper (optional) --------

detect_communities <- function(graph, algorithm = c(
  "walktrap", "fastgreedy", "louvain",
  "labelpropagation", "infomap", "leadingeigenvector"
)) {
  algorithm <- match.arg(algorithm)
  if (algorithm == "walktrap") {
    cluster_walktrap(graph)
  } else if (algorithm == "fastgreedy") {
    cluster_fast_greedy(as.undirected(graph))
  } else if (algorithm == "louvain") {
    cluster_louvain(as.undirected(graph))
  } else if (algorithm == "labelpropagation") {
    cluster_label_prop(as.undirected(graph))
  } else if (algorithm == "infomap") {
    cluster_infomap(graph)
  } else if (algorithm == "leadingeigenvector") {
    cluster_leading_eigen(as.undirected(graph))
  }
}

## -------------------- INFOMAP: community size distributions (optional) ------

par(mfrow = c(3, 3))  # 3x3 grid; adjust to your number of graphs

for (nm in names(graphs)) {
  g_current <- graphs[[nm]]
  
  com_infomap <- cluster_infomap(g_current)
  community_sizes <- sizes(com_infomap)
  
  hist(
    community_sizes,
    main = paste("Community Sizes (Infomap) -", nm),
    xlab = "Size", ylab = "Frequency",
    col = "lightblue", breaks = 7
  )
  print(table(community_sizes))
}

## -------------------- LABEL PROPAGATION: community sizes (optional) ---------

par(mfrow = c(3, 3))  # reset grid if needed

for (nm in names(graphs)) {
  g_current <- graphs[[nm]]
  
  com_label <- cluster_label_prop(as.undirected(g_current))
  community_sizes <- sizes(com_label)
  
  hist(
    community_sizes,
    main = paste("Community Sizes (LabelProp) -", nm),
    xlab = "Size", ylab = "Frequency",
    col = "lightblue", breaks = 7
  )
  print(table(community_sizes))
}

## -------------------- WALKTRAP communities (main) ---------------------------

# We use walktrap for the main downstream GO analysis
community_list_walktrap <- lapply(graphs, cluster_walktrap)
names(community_list_walktrap) <- paste0("CLUSTER_WALKTRAP_", names(graphs))

## -------------------- GO enrichment on communities --------------------------

# overlap/Jaccard utilities kept in case needed later
overlap_coefficient <- function(set1, set2) {
  intersection_size <- length(intersect(set1, set2))
  intersection_size / sqrt(length(set1) * length(set2))
}

jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union_size   <- length(union(set1, set2))
  intersection / union_size
}

# Generic GO-enrichment over communities
perform_GO_enrichment <- function(community_obj,
                                  label_prefix,
                                  organism = "hsapiens",
                                  sources = c("GO:BP"),
                                  pval_cutoff = 0.05) {
  membership_vec <- membership(community_obj)
  results_list <- list()
  
  for (group_id in unique(membership_vec)) {
    group_members <- names(which(membership_vec == group_id))
    
    if (length(group_members) > 10) {
      gostres <- try(
        gost(
          query            = group_members,
          exclude_iea      = TRUE,
          user_threshold   = pval_cutoff,
          correction_method = "fdr",
          domain_scope     = "annotated",
          sources          = sources,
          organism         = organism
        ),
        silent = TRUE
      )
      
      if (inherits(gostres, "try-error") || is.null(gostres$result)) {
        next
      }
      
      gprofout <- gostres$result[, c("p_value", "term_id", "term_name")]
      gprofout <- gprofout[gprofout$p_value < pval_cutoff, , drop = FALSE]
      
      if (nrow(gprofout) > 0) {
        gprofout <- gprofout[1, , drop = FALSE]
        results_list[[paste0(label_prefix, "_Group", group_id)]] <- gprofout
      }
    }
  }
  
  results_list
}

# GOBP enrichment per walktrap-community object
all_results_GOBP <- lapply(
  names(community_list_walktrap),
  function(lbl) {
    perform_GO_enrichment(
      community_obj = community_list_walktrap[[lbl]],
      label_prefix  = lbl,
      sources       = c("GO:BP")
    )
  }
)
names(all_results_GOBP) <- names(community_list_walktrap)

# GOMF enrichment per walktrap-community object
all_results_GOMF <- lapply(
  names(community_list_walktrap),
  function(lbl) {
    perform_GO_enrichment(
      community_obj = community_list_walktrap[[lbl]],
      label_prefix  = lbl,
      sources       = c("GO:MF")
    )
  }
)
names(all_results_GOMF) <- names(community_list_walktrap)

## -------------------- Enrichment restricted to communities with backbone ----

# Backbone genes (ideally should match script 02; here we define explicitly
# or you can source backbone.bois from script 02 output)
backbone_bois <- c("BCL6","CREM","HMGB2","MAF","PBX4","RBPJ","STAT4","ZNF292")

perform_GO_enrichment_with_backbone_presence <- function(community_obj,
                                                         backbone_genes,
                                                         label_prefix,
                                                         organism = "hsapiens",
                                                         sources = c("GO:MF"),
                                                         pval_cutoff = 0.05) {
  membership_vec <- membership(community_obj)
  results_list <- list()
  
  for (group_id in unique(membership_vec)) {
    group_members <- names(which(membership_vec == group_id))
    present_backbone_genes <- intersect(group_members, backbone_genes)
    
    if (length(present_backbone_genes) > 0 && length(group_members) > 10) {
      gostres <- try(
        gost(
          query            = group_members,
          exclude_iea      = TRUE,
          user_threshold   = pval_cutoff,
          correction_method = "fdr",
          domain_scope     = "annotated",
          sources          = sources,
          organism         = organism
        ),
        silent = TRUE
      )
      
      if (inherits(gostres, "try-error") || is.null(gostres$result)) {
        next
      }
      
      gprofout <- gostres$result[, c("p_value", "term_id", "term_name")]
      gprofout <- gprofout[gprofout$p_value < pval_cutoff, , drop = FALSE]
      
      if (nrow(gprofout) > 0) {
        gprofout <- gprofout[1, , drop = FALSE]
        results_list[[paste0(label_prefix, "_Group", group_id)]] <- list(
          GO_Results              = gprofout,
          Genes                   = group_members,
          Present_Backbone_Drivers = present_backbone_genes
        )
      }
    }
  }
  
  results_list
}

all_results_with_backbone_presence <- lapply(
  names(community_list_walktrap),
  function(lbl) {
    perform_GO_enrichment_with_backbone_presence(
      community_obj  = community_list_walktrap[[lbl]],
      backbone_genes = backbone_bois,
      label_prefix   = lbl,
      sources        = c("GO:MF")
    )
  }
)
names(all_results_with_backbone_presence) <- names(community_list_walktrap)

## -------------------- Enrichment restricted to communities with drivers -----

NKT_drivers        <- c("ZNF732", "RPS8", "KLRB1", "GNPTAB")
Th1_drivers        <- c("EOMES", "CCL5", "CENPE", "SERPINB9", "PAM", "SSBP2", "AHI1", "DTNB")
Th2_drivers        <- c("KLF2", "ABI1", "ARHGEF3")
Th17_drivers       <- c("HLF", "PVT1", "HMGCR", "MAST4", "LRRFIP2", "EDEM3", "RANBP9",
                        "MAP3K4", "GABPB1", "SLC16A7", "SLC9A9", "PHLPP1", "PBX1")
Th1_17_drivers     <- c("HLF", "EOMES", "GMNN", "ZC3HAV1")
mTreg_drivers      <- c("SOX13", "PLIN2", "TXNIP", "BCAS3", "TMEM260")
nTreg_drivers      <- c("TFEC", "DACH1", "ESR1", "DUSP4", "ST6GALNAC3", "CARMIL1")
naive.stim_drivers <- c("MYB", "CFLAR", "MBD5", "STAMBPL1", "HEATR5A", "KCNQ5")

# collect all driver genes into one vector (your original list_drivers logic)
driver_genes_all <- unique(c(
  NKT_drivers, Th1_drivers, Th2_drivers, Th17_drivers, Th1_17_drivers,
  mTreg_drivers, nTreg_drivers, naive.stim_drivers
))

perform_GO_enrichment_with_driver_presence <- function(community_obj,
                                                       driver_genes,
                                                       label_prefix,
                                                       organism = "hsapiens",
                                                       sources = c("GO:MF"),
                                                       pval_cutoff = 0.05) {
  membership_vec <- membership(community_obj)
  results_list <- list()
  
  for (group_id in unique(membership_vec)) {
    group_members <- names(which(membership_vec == group_id))
    present_drivers <- intersect(group_members, driver_genes)
    
    if (length(present_drivers) > 0 && length(group_members) > 10) {
      gostres <- try(
        gost(
          query            = group_members,
          exclude_iea      = TRUE,
          user_threshold   = pval_cutoff,
          correction_method = "fdr",
          domain_scope     = "annotated",
          sources          = sources,
          organism         = organism
        ),
        silent = TRUE
      )
      
      if (inherits(gostres, "try-error") || is.null(gostres$result)) {
        next
      }
      
      gprofout <- gostres$result[, c("p_value", "term_id", "term_name")]
      gprofout <- gprofout[gprofout$p_value < pval_cutoff, , drop = FALSE]
      
      if (nrow(gprofout) > 0) {
        gprofout <- gprofout[1, , drop = FALSE]
        results_list[[paste0(label_prefix, "_Group", group_id)]] <- list(
          GO_Results              = gprofout,
          Genes                   = group_members,
          Present_Backbone_Drivers = present_drivers
        )
      }
    }
  }
  
  results_list
}

all_results_with_driver_presence <- lapply(
  names(community_list_walktrap),
  function(lbl) {
    perform_GO_enrichment_with_driver_presence(
      community_obj = community_list_walktrap[[lbl]],
      driver_genes  = driver_genes_all,
      label_prefix  = lbl,
      sources       = c("GO:MF")
    )
  }
)
names(all_results_with_driver_presence) <- names(community_list_walktrap)

## -------------------- Helper: extract info for a given gene -----------------

extract_info_for_gene <- function(nested_list, target_gene) {
  extracted_info <- data.frame(
    Group                  = character(),
    GO_Results             = I(list()),
    Genes                  = I(list()),
    Present_Backbone_Drivers = I(list()),
    stringsAsFactors       = FALSE
  )
  
  for (outer in nested_list) {
    for (inner_name in names(outer)) {
      inner_list <- outer[[inner_name]]
      if (target_gene %in% inner_list$Genes) {
        extracted_info <- rbind(
          extracted_info,
          data.frame(
            Group                  = inner_name,
            GO_Results             = I(list(inner_list$GO_Results)),
            Genes                  = I(list(inner_list$Genes)),
            Present_Backbone_Drivers = I(list(inner_list$Present_Backbone_Drivers)),
            stringsAsFactors       = FALSE
          )
        )
      }
    }
  }
  
  extracted_info
}

## Example usage:
# a <- extract_info_for_gene(all_results_with_backbone_presence, "HMGB2")
# a

## ============================================================================
## 4) Overlap of communities containing backbones / drivers (SS coefficient)
## ============================================================================

## ---- 4.1 Szymkiewicz–Simpson coefficient -----------------------------------

szymkiewicz_simpson <- function(set1, set2) {
  if (length(set1) == 0 || length(set2) == 0) return(0)
  length(intersect(set1, set2)) / min(length(set1), length(set2))
}

## ---- 4.2 Overlap of GO-enriched backbone communities -----------------------
## Uses: all_results_with_backbone_presence (list of lists produced above)

extract_and_compare_genes_v1 <- function(all_results) {
  group_genes   <- list()
  group_drivers <- list()
  
  # all_results: list (per walktrap-cluster), each an inner list of groups
  for (result_list in all_results) {
    for (group_name in names(result_list)) {
      group_genes[[group_name]]   <- result_list[[group_name]]$Genes
      group_drivers[[group_name]] <- result_list[[group_name]]$Present_Backbone_Drivers
    }
  }
  
  group_names <- names(group_genes)
  
  # build more informative labels: subtype + backbone drivers as suffix
  modified_group_names <- sapply(group_names, function(name) {
    short_name <- gsub("^CLUSTER_WALKTRAP_", "", name)
    short_name <- sub("_Group[0-9]+.*$", "", short_name)
    
    drv_vec <- group_drivers[[name]]
    drv_vec <- drv_vec[!is.na(drv_vec) & nzchar(drv_vec)]
    driver_suffix <- if (length(drv_vec) > 0) {
      paste0("_", paste(unique(drv_vec), collapse = "_"))
    } else {
      ""
    }
    
    paste0(short_name, driver_suffix)
  })
  
  # overlap matrix (Szymkiewicz–Simpson)
  overlap_matrix <- matrix(
    0,
    nrow = length(group_genes),
    ncol = length(group_genes),
    dimnames = list(modified_group_names, modified_group_names)
  )
  
  for (i in seq_along(group_genes)) {
    for (j in seq_along(group_genes)) {
      overlap_matrix[i, j] <- szymkiewicz_simpson(
        group_genes[[group_names[i]]],
        group_genes[[group_names[j]]]
      )
    }
  }
  
  overlap_matrix
}

overlap_matrix_backbone_GO <- extract_and_compare_genes_v1(all_results_with_backbone_presence)

# Heatmap
melted_overlap_matrix <- reshape2::melt(overlap_matrix_backbone_GO)

ggplot(melted_overlap_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0.5, limit = c(0, 1), space = "Lab",
    name = "Overlap\nCoefficient"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(vjust = 1, face = "bold"),
    axis.title  = element_blank(),
    legend.title.align = 0.5
  ) +
  coord_fixed()

## Optional:
# write.csv(
#   melted_overlap_matrix,
#   "/corgi/sebas/tcell_multi/overlap_matrix_backbones_drivers/overlap_mtx_backbone_GO_SS.csv",
#   row.names = FALSE
# )

## ---- 4.3 Overlap on all backbone-containing communities (no GO filter) -----
## Uses: community_list_walktrap + backbone_bois + subtype-specific drivers

# we already have NKT_drivers, Th1_drivers, ..., naive.stim_drivers
list_drivers <- grep("_drivers$", ls(), value = TRUE)

# collect all communities that:
#   - have >10 genes
#   - contain at least one backbone TF
build_all_group_genes <- function(community_list, backbone_genes) {
  out <- list()
  for (lbl in names(community_list)) {
    comm <- community_list[[lbl]]
    membership_vec <- membership(comm)
    
    for (group_id in unique(membership_vec)) {
      group_members <- names(which(membership_vec == group_id))
      if (length(group_members) > 10 &&
          length(intersect(group_members, backbone_genes)) > 0) {
        gname <- paste0(lbl, "_Group", group_id)
        out[[gname]] <- group_members
      }
    }
  }
  out
}

all_group_genes <- build_all_group_genes(
  community_list_walktrap,
  backbone_bois
)

group_names <- names(all_group_genes)
overlap_matrix_groups <- matrix(
  0,
  nrow = length(group_names),
  ncol = length(group_names),
  dimnames = list(group_names, group_names)
)

for (i in seq_along(group_names)) {
  for (j in seq_along(group_names)) {
    if (i != j) {
      overlap_matrix_groups[i, j] <- szymkiewicz_simpson(
        all_group_genes[[group_names[i]]],
        all_group_genes[[group_names[j]]]
      )
    }
  }
}

# helper to build labels including subtype + backbone + drivers + group size
create_new_name_with_drivers <- function(group_name,
                                         all_group_genes,
                                         backbone_genes,
                                         list_drivers) {
  group_backbones <- intersect(all_group_genes[[group_name]], backbone_genes)
  group_size      <- length(all_group_genes[[group_name]])
  
  # subtype from group_name: "CLUSTER_WALKTRAP_<subtype>_graph_obj_GroupX"
  subtype <- sub("^CLUSTER_WALKTRAP_", "", group_name)
  subtype <- sub("_graph_obj_Group[0-9]+$", "", subtype)
  
  driver_obj_name <- paste0(subtype, "_drivers")
  if (driver_obj_name %in% list_drivers) {
    group_drivers <- get(driver_obj_name)
  } else {
    group_drivers <- character(0)
  }
  
  present_drivers <- intersect(all_group_genes[[group_name]], group_drivers)
  combined_genes  <- unique(c(group_backbones, present_drivers))
  
  if (length(combined_genes) > 0) {
    paste0(subtype, "_", paste(combined_genes, collapse = "_"), "_", group_size)
  } else {
    paste0(subtype, "_", group_size)
  }
}

new_names_with_drivers <- sapply(
  names(all_group_genes),
  create_new_name_with_drivers,
  all_group_genes,
  backbone_bois,
  list_drivers
)

rownames(overlap_matrix_groups) <- new_names_with_drivers
colnames(overlap_matrix_groups) <- new_names_with_drivers

melted_overlap_matrix_groups <- reshape2::melt(overlap_matrix_groups)

ggplot(melted_overlap_matrix_groups, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0.5, limit = c(0, 1), space = "Lab",
    name = "Overlap\nCoefficient"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(vjust = 1, face = "bold"),
    axis.title  = element_blank(),
    legend.title.align = 0.5
  ) +
  coord_fixed()

## Optional:
# write.csv(
#   melted_overlap_matrix_groups,
#   "/corgi/sebas/tcell_multi/overlap_matrix_backbones_drivers/overlap_mtx_backbone_allComm_SS.csv",
#   row.names = FALSE
# )

## ---- 4.4 Combined overlap: drivers + backbones + GO context ----------------
## Uses: all_results_with_driver_presence + all_results_with_backbone_presence

extract_and_compare_genes <- function(all_results_with_driver_presence,
                                      all_results_with_backbone_presence) {
  group_genes    <- list()
  group_drivers  <- list()
  group_go_terms <- list()
  
  # from backbone_presence: genes, backbone drivers, GO terms
  for (result_list in all_results_with_backbone_presence) {
    for (group_name in names(result_list)) {
      group_genes[[group_name]]   <- result_list[[group_name]]$Genes
      group_drivers[[group_name]] <- result_list[[group_name]]$Present_Backbone_Drivers
      
      go_results <- result_list[[group_name]]$GO_Results
      if (!is.null(go_results) && nrow(go_results) > 0) {
        group_go_terms[[group_name]] <- go_results$term_id[1]
      } else {
        group_go_terms[[group_name]] <- NA_character_
      }
    }
  }
  
  # genes of interest from driver presence
  genes_of_interest <- unique(unlist(
    lapply(all_results_with_driver_presence, function(x) {
      unlist(lapply(x, function(y) y$Present_Backbone_Drivers))
    })
  ))
  
  group_names <- names(group_genes)
  
  modified_group_names <- sapply(group_names, function(name) {
    short_name <- gsub("^CLUSTER_WALKTRAP_", "", name)
    short_name <- sub("_Group[0-9]+.*$", "", short_name)
    
    drv_vec <- unique(c(
      group_drivers[[name]],
      genes_of_interest[genes_of_interest %in% group_genes[[name]]]
    ))
    drv_vec <- drv_vec[!is.na(drv_vec) & nzchar(drv_vec)]
    
    driver_suffix <- if (length(drv_vec) > 0) {
      paste0("_", paste(drv_vec, collapse = "_"))
    } else {
      ""
    }
    
    go_suffix <- if (!is.na(group_go_terms[[name]])) {
      paste0("_", group_go_terms[[name]])
    } else {
      ""
    }
    
    paste0(short_name, driver_suffix, go_suffix)
  })
  
  overlap_matrix <- matrix(
    0,
    nrow = length(group_genes),
    ncol = length(group_genes),
    dimnames = list(modified_group_names, modified_group_names)
  )
  
  for (i in seq_along(group_genes)) {
    for (j in seq_along(group_genes)) {
      overlap_matrix[i, j] <- szymkiewicz_simpson(
        group_genes[[group_names[i]]],
        group_genes[[group_names[j]]]
      )
    }
  }
  
  overlap_matrix
}

overlap_matrix2 <- extract_and_compare_genes(
  all_results_with_driver_presence,
  all_results_with_backbone_presence
)

melted_overlap_matrix2 <- reshape2::melt(overlap_matrix2)

# write.csv(
#   melted_overlap_matrix2,
#   "/corgi/sebas/tcell_multi/overlap_matrix_backbones_drivers/overlap_mtx_driv_bb_MF.csv",
#   row.names = FALSE
# )

ggplot(melted_overlap_matrix2, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0.5, limit = c(0, 1), space = "Lab",
    name = "Overlap\nCoefficient"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(vjust = 1, face = "bold"),
    axis.title  = element_blank(),
    legend.title.align = 0.5
  ) +
  coord_fixed()




