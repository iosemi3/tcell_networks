
#!/usr/bin/env Rscript
## ============================================================================
## 02_backbone_activation_effector_analysis.R
##
## Standalone script:
##   - consumes frozen 2023 degree-centrality matrix
##     (merged_graph_obj_degree_centrality.csv)
##   - derives backbone TFs algorithmically:
##       * TFs with no zero degree in any subtype
##       * ranked by mean degree centrality across subtypes
##       * top 8 retained as backbone.bois
##   - defines activation, stimulation, effector TF sets
##   - prepares Cytoscape gene lists (exports commented)
##   - plots backbone heatmap
##
## This script does NOT recompute networks or centrality.
## ============================================================================

library(dplyr)
library(reshape2)
library(ggplot2)

## -------------------- Paths / parameters ------------------------------------

# Path to centrality matrix
merged_cents <- read.csv(
  "/Users/sebas/Desktop/tcell network paper submission/csv_files/merged_graph_obj_degree_centrality.csv",
  header = TRUE,
  check.names = FALSE
)
# First column = rownames, second column = gene id; restore original structure
merged_cents <- merged_cents[, -1]
rownames(merged_cents) <- merged_cents[, 1]
merged_cents <- merged_cents[, -1]

# Ensure we’re working with numeric centralities
num_mat <- as.matrix(merged_cents)
if (!is.numeric(num_mat)) {
  stop("merged_cents contains non-numeric entries; please inspect the input file.")
}

## Optional sanity check: expected columns
required_cols <- c("naive", "naive.stim", "Th2")
missing_cols  <- setdiff(required_cols, colnames(merged_cents))
if (length(missing_cols) > 0) {
  warning("Missing expected columns in merged_cents: ",
          paste(missing_cols, collapse = ", "))
}

## -------------------- Derive backbone.bois from centrality matrix ----------

# TFs that are central (non-zero) across all subtypes
no_zero <- rowSums(num_mat == 0) == 0
central_tfs <- merged_cents[no_zero, , drop = FALSE]

# Rank by mean centrality across subtypes
central_tfs$mean_cent <- rowMeans(central_tfs, na.rm = TRUE)
central_tfs <- central_tfs[order(-central_tfs$mean_cent), ]

# Top 8 = backbone TFs
backbone.bois <- rownames(central_tfs)[1:8]
cat("Backbone TFs:\n")
print(backbone.bois)

## -------------------- Backbone / activation / stimulation / effector --------

# backbone TFs: rows corresponding to backbone.bois
backbones <- merged_cents[rownames(merged_cents) %in% backbone.bois, ]

# activation TFs: not central in naive (naive == 0)
activation <- merged_cents[merged_cents$naive == 0, ]

# stimulation TFs: not central in naive AND not central in naive.stim
stimulation <- merged_cents[
  merged_cents$naive == 0 & merged_cents$naive.stim == 0,
]

# effector TFs: in 'stimulation' and central in exactly one subtype
effector <- stimulation[rowSums(stimulation != 0) == 1, ]

## -------------------- Cytoscape exports (all write.csv commented) -----------

#### Activation: top 6 genes by naive.stim

sorted_stim <- activation[order(-activation$naive.stim), ]
sorted_stim <- sorted_stim[1:6, ]
sorted_activation <- rownames(sorted_stim)

# write.csv(
#   sorted_activation,
#   "/corgi/sebas/tcell_multi/cytoscape/activation_genes.csv",
#   row.names = FALSE,
#   quote = FALSE
# )

#### Stimulation: top 6 genes by Th2

sorted_stim <- stimulation[order(-stimulation$Th2), ]
sorted_stim <- sorted_stim[1:6, ]
sorted_stim <- rownames(sorted_stim)

# write.csv(
#   sorted_stim,
#   "/corgi/sebas/tcell_multi/cytoscape/stimulation_genes.csv",
#   row.names = FALSE,
#   quote = FALSE
# )

#### Effector: per-subtype top 3, plus per-subtype top 6 exports

sorted_effectors_all_st <- character()

# NOTE: the index 1:9 assumes the same column order as in the 2023 export.
# This preserves the original behavior exactly.
for (i in 1:9) {  ### overlapping of effectors, take top 3 for each subtype
  idx <- which(effector[, i] > 0)
  if (length(idx) == 0) next
  
  eff_vals  <- effector[idx, i]
  eff_genes <- rownames(effector)[idx][order(eff_vals, decreasing = TRUE)]
  effectors <- head(eff_genes, 3)
  
  sorted_effectors_all_st <- c(sorted_effectors_all_st, effectors)
}

# For overlapping effector TFs across subtypes
cat("Combined effector TFs (top 3 per subtype):\n")
print(sorted_effectors_all_st)

## Per-subtype top 6 effector genes (column indices follow original order)

# Th1
# write.csv(
#   head(rownames(effector)[which(effector[, 1] > 0)][
#     order(effector[effector[, 1] > 0, 1], decreasing = TRUE)
#   ], 6),
#   "/corgi/sebas/tcell_multi/cytoscape/effector_th1_genes.csv",
#   row.names = FALSE, quote = FALSE
# )

# Th2
# write.csv(
#   head(rownames(effector)[which(effector[, 2] > 0)][
#     order(effector[effector[, 2] > 0, 2], decreasing = TRUE)
#   ], 6),
#   "/corgi/sebas/tcell_multi/cytoscape/effector_th2_genes.csv",
#   row.names = FALSE, quote = FALSE
# )

# Th1/17
# write.csv(
#   head(rownames(effector)[which(effector[, 3] > 0)][
#     order(effector[effector[, 3] > 0, 3], decreasing = TRUE)
#   ], 6),
#   "/corgi/sebas/tcell_multi/cytoscape/effector_th1_17_genes.csv",
#   row.names = FALSE, quote = FALSE
# )

# Th17
# write.csv(
#   head(rownames(effector)[which(effector[, 4] > 0)][
#     order(effector[effector[, 4] > 0, 4], decreasing = TRUE)
#   ], 6),
#   "/corgi/sebas/tcell_multi/cytoscape/effector_th17_genes.csv",
#   row.names = FALSE, quote = FALSE
# )

# NKT
# write.csv(
#   head(rownames(effector)[which(effector[, 5] > 0)][
#     order(effector[effector[, 5] > 0, 5], decreasing = TRUE)
#   ], 6),
#   "/corgi/sebas/tcell_multi/cytoscape/effector_NKT_genes.csv",
#   row.names = FALSE, quote = FALSE
# )

# naive.stim
# write.csv(
#   head(rownames(effector)[which(effector[, 6] > 0)][
#     order(effector[effector[, 6] > 0, 6], decreasing = TRUE)
#   ], 6),
#   "/corgi/sebas/tcell_multi/cytoscape/effector_naive.stim_genes.csv",
#   row.names = FALSE, quote = FALSE
# )

# Treg.mem
# write.csv(
#   head(rownames(effector)[which(effector[, 7] > 0)][
#     order(effector[effector[, 7] > 0, 7], decreasing = TRUE)
#   ], 6),
#   "/corgi/sebas/tcell_multi/cytoscape/effector_Treg.mem_genes.csv",
#   row.names = FALSE, quote = FALSE
# )

# Treg.naive
# write.csv(
#   head(rownames(effector)[which(effector[, 8] > 0)][
#     order(effector[effector[, 8] > 0, 8], decreasing = TRUE)
#   ], 6),
#   "/corgi/sebas/tcell_multi/cytoscape/effector_Treg.naive_genes.csv",
#   row.names = FALSE, quote = FALSE
# )

# naive
# write.csv(
#   head(rownames(effector)[which(effector[, 9] > 0)][
#     order(effector[effector[, 9] > 0, 9], decreasing = TRUE)
#   ], 6),
#   "/corgi/sebas/tcell_multi/cytoscape/effector_naive_genes.csv",
#   row.names = FALSE, quote = FALSE
# )

## -------------------- Backbone heatmap --------------------------------------

backbones$gene <- rownames(backbones)

df_top <- melt(backbones, id.vars = "gene")

ggplot(df_top, aes(x = variable, y = gene)) +
  geom_tile(aes(fill = value), color = "black", size = 0.5) +
  scale_fill_gradient(low = "cornsilk", high = "brown3") +
  labs(
    y    = "Backbone TF",
    x    = "Subtype",
    fill = "Degree Centrality"
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(
    panel.grid  = element_line(size = 1, linewidth = 0.5, colour = "black"),
    axis.text   = element_text(face = "bold", size = 14),
    axis.title  = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold")
  )
