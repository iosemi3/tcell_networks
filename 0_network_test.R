#!/usr/bin/env Rscript
## ============================================================================
## network_validation.R
## Validation of xgboost GRN centrality against random graph models
## - Uses tbl_graphs.RDS (xgboost, no cell-cycle regression)
## - Generates:
##    <name>Wdegdist_real / <name>Wdegdist_sim  (WS model)
##    <name>Wdegcent_real / <name>Wdegcent_sim
##    <name>degdist_real / <name>degdist_sim    (BA model)
##    <name>degcent_real / <name>degcent_sim
## - Performs permutation tests, Wilcoxon tests, t-tests
## - All write.csv() are commented out (restore if export is needed)
## ============================================================================

library(igraph)
library(ggplot2)

set.seed(666)

## ------------------------- Load graph objects -------------------------------

objects <- readRDS("/Users/sebas/Desktop/tcell network paper submission/rds_files/tbl_graphs.RDS")
# 'objects' is expected to be a named list of graph-like objects (igraph or tbl_graph)

## ------------------------- Helper functions ---------------------------------

plot_density_overlay <- function(data1, data2, main = "Density Plot") {
  plot(density(data1), main = main, col = "green", lwd = 2)
  lines(density(data2), col = "red", lwd = 2)
}

permutation_test_and_plot <- function(real_deg_dist, sim_deg_dist, n_permutations = 1000) {
  set.seed(666)
  real_var_name <- deparse(substitute(real_deg_dist))
  
  observed_diff <- mean(real_deg_dist) - mean(sim_deg_dist)
  pooled_data   <- c(real_deg_dist, sim_deg_dist)
  
  permuted_diffs <- numeric(n_permutations)
  for (perm_idx in seq_len(n_permutations)) {
    permuted_data  <- sample(pooled_data)
    permuted_real  <- permuted_data[seq_along(real_deg_dist)]
    permuted_sim   <- permuted_data[(length(real_deg_dist) + 1):length(pooled_data)]
    permuted_diffs[perm_idx] <- mean(permuted_real) - mean(permuted_sim)
  }
  
  p_value <- mean(abs(permuted_diffs) >= abs(observed_diff))
  
  common_xlim <- range(c(real_deg_dist, sim_deg_dist))
  
  hist(sim_deg_dist,
       breaks = 20, col = rgb(1, 0, 0, 0.5),
       xlim = common_xlim, freq = FALSE,
       xlab = "Degree", ylab = "Density", main = real_var_name)
  
  hist(real_deg_dist,
       breaks = 20, col = rgb(0, 1, 0, 0.7),
       add = TRUE, freq = FALSE)
  
  abline(v = mean(real_deg_dist), col = "green", lty = 2, lwd = 2)
  abline(v = mean(sim_deg_dist),  col = "red",   lty = 2, lwd = 2)
  
  legend_text <- c(
    "Real", "Simulated", "Real Avg", "Sim Avg",
    paste0("p-value: ", round(p_value, 4)),
    paste0("Permutations: ", n_permutations)
  )
  legend("topright",
         legend = legend_text,
         col    = c(rgb(0,1,0,0.5), rgb(1,0,0,0.5), "green", "red", "black", "black"),
         pch    = c(15, 15, NA, NA, NA, NA),
         lty    = c(NA, NA, 2, 2, NA, NA),
         text.font = 2)
  
  return(p_value)
}

wilcoxon_test_and_plot <- function(real_deg_dist, sim_deg_dist) {
  test_result <- wilcox.test(real_deg_dist, sim_deg_dist, alternative = "two.sided")
  p_value <- test_result$p.value
  real_var_name <- deparse(substitute(real_deg_dist))
  
  common_xlim <- range(c(real_deg_dist, sim_deg_dist))
  
  hist(sim_deg_dist,
       breaks = 10, col = rgb(1, 0, 0, 0.5),
       xlim = common_xlim, freq = FALSE,
       xlab = "Degree", ylab = "Density", main = real_var_name)
  
  hist(real_deg_dist,
       breaks = 10, col = rgb(0, 1, 0, 0.7),
       add = TRUE, freq = FALSE)
  
  abline(v = median(real_deg_dist), col = "green", lty = 2, lwd = 2)
  abline(v = median(sim_deg_dist),  col = "red",   lty = 2, lwd = 2)
  
  legend_text <- c(
    "Real", "Simulated", "Real Median", "Sim Median",
    paste0("p-value: ", signif(p_value, 4))
  )
  legend("topright",
         legend = legend_text,
         col    = c(rgb(0,1,0,0.5), rgb(1,0,0,0.5), "green", "red", "black"),
         pch    = c(15, 15, NA, NA, NA),
         lty    = c(NA, NA, 2, 2, NA),
         text.font = 2)
  
  return(p_value)
}

t_test_and_plot <- function(real_deg_dist, sim_deg_dist) {
  test_result <- t.test(real_deg_dist, sim_deg_dist, alternative = "two.sided")
  p_value <- test_result$p.value
  
  real_mean <- formatC(mean(real_deg_dist), format = "f", digits = 3)
  sim_mean  <- formatC(mean(sim_deg_dist),  format = "f", digits = 3)
  
  real_var_name <- deparse(substitute(real_deg_dist))
  common_xlim <- range(c(real_deg_dist, sim_deg_dist))
  
  hist(sim_deg_dist,
       breaks = 5, col = rgb(1,0,0,0.5),
       xlim = common_xlim, freq = FALSE,
       xlab = "Degree", ylab = "Density", main = real_var_name)
  
  hist(real_deg_dist,
       breaks = 20, col = rgb(0,1,0,0.7),
       add = TRUE, freq = FALSE)
  
  abline(v = mean(real_deg_dist), col = "green", lty = 2, lwd = 2)
  abline(v = mean(sim_deg_dist),  col = "red",   lty = 2, lwd = 2)
  
  legend_text <- c(
    "Real", "Simulated",
    paste("Real Mean:", real_mean),
    paste("Sim Mean:",  sim_mean),
    paste0("p-value: ", signif(p_value, 4))
  )
  legend("topright",
         legend = legend_text,
         col    = c(rgb(0,1,0,0.5), rgb(1,0,0,0.5), "green", "red", "black"),
         pch    = c(15, 15, NA, NA, NA),
         lty    = c(NA, NA, 2, 2, NA),
         text.font = 2)
  
  return(p_value)
}

perform_hypothesis_test <- function(gl_centralities, g_centralities_real, p_value_threshold = 0.05) {
  global_mean_gl <- mean(gl_centralities)
  
  t_test_result <- t.test(g_centralities_real,
                          mu  = global_mean_gl,
                          alternative = "two.sided",
                          conf.level  = 0.95)
  
  df <- data.frame(Value = g_centralities_real)
  hist_plot <- ggplot(df, aes(x = Value)) +
    geom_histogram(binwidth = 0.02, fill = "blue", color = "black") +
    geom_vline(xintercept = global_mean_gl, color = "red",   linetype = "dashed", size = 1) +
    geom_vline(xintercept = quantile(g_centralities_real, 1 - p_value_threshold),
               color = "green", linetype = "dashed", size = 1) +
    labs(x = "Eigenvector Centrality") +
    theme_minimal()
  
  print(hist_plot)
  
  if (t_test_result$p.value <= p_value_threshold) {
    cat("Reject H0: means differ (p =", t_test_result$p.value, ")\n")
  } else {
    cat("Fail to reject H0: means not significantly different (p =", t_test_result$p.value, ")\n")
  }
  
  invisible(t_test_result$p.value)
}

## -------------------- NAIVE network: Watts–Strogatz -------------------------
## NOTE: assumes naive networks have "naive" in the name.
## If this assumption is wrong, edit the grepl("naive", ...) filters.

for (stp in names(objects)) {
  if (!grepl("naive", stp, ignore.case = TRUE)) next
  
  i <- stp  # keep your naming
  current_object <- objects[[stp]]
  g <- as.igraph(current_object)
  
  degree_distribution_values_g <- degree_distribution(g)
  degree_g <- degree(g)
  
  namecnt <- sub("(.*)_graph_obj", "\\1", stp)
  message("NAIVE / WS model: ", namecnt)
  
  name1 <- paste0(namecnt, "Wdegdist_real")
  name2 <- paste0(namecnt, "Wdegcent_real")
  
  assign(name1, degree_distribution_values_g)
  assign(name2, degree_g)
  
  # simulate 10,000 random WS graphs
  n   <- vcount(g)
  avg_k <- mean(neighborhood.size(g, order = 1, mode = "all"))
  avg_edge_density <- edge_density(g)
  
  gl <- vector("list", 10000)
  for (rep_idx in seq_len(10000)) {
    random_network <- watts.strogatz.game(1, n, avg_k, avg_edge_density)
    V(random_network)$name <- V(g)$name
    gl[[rep_idx]] <- random_network
  }
  
  gl_degdists <- unlist(lapply(gl, degree_distribution, mode = "total"))
  name2 <- paste0(namecnt, "Wdegdist_sim")
  assign(name2, gl_degdists)
  message("Stored: ", name2)
  
  gl_degcents <- unlist(lapply(gl, degree, mode = "total"))
  name2 <- paste0(namecnt, "Wdegcent_sim")
  assign(name2, gl_degcents)
  message("Stored: ", name2)
}

## Example usage (replace with your subtype-specific names):
# plot_density_overlay(naiveWdegdist_real, naiveWdegdist_sim)
# p_value_perm <- permutation_test_and_plot(naiveWdegdist_real, naiveWdegdist_sim, 1000)

## -------------------- ALL OTHER subtypes: Barabási–Albert -------------------

for (stp in names(objects)) {
  if (grepl("naive", stp, ignore.case = TRUE)) next
  
  i <- stp  # keep your naming
  current_object <- objects[[stp]]
  g <- as.igraph(current_object)
  
  degree_distribution_values_g <- degree_distribution(g)
  degree_g <- degree(g)
  
  namecnt <- sub("(.*)_graph_obj", "\\1", stp)
  message("BA model: ", namecnt)
  
  name1 <- paste0(namecnt, "degdist_real")
  name2 <- paste0(namecnt, "degcent_real")
  
  assign(name1, degree_distribution_values_g)
  assign(name2, degree_g)
  
  gl <- vector("list", 10000)
  for (rep_idx in seq_len(10000)) {
    m_param <- max(1L, floor(ecount(g) / gorder(g)))
    gl[[rep_idx]] <- barabasi.game(n = gorder(g), m = m_param, directed = TRUE)
    V(gl[[rep_idx]])$name <- V(g)$name
  }
  
  gl_degdists <- unlist(lapply(gl, degree_distribution, mode = "total"))
  name2 <- paste0(namecnt, "degdist_sim")
  assign(name2, gl_degdists)
  message("Stored: ", namecnt, "degdist_sim")
  
  gl_degcents <- unlist(lapply(gl, degree, mode = "total"))
  name2 <- paste0(namecnt, "degcent_sim")
  assign(name2, gl_degcents)
  message("Stored: ", namecnt, "degcent_sim")
}

## -------------------- Example tests / plots (manually run) ------------------
## DEGREE distributions:
## wilcoxon_test_and_plot(naive_stimdegdist_real, naive_stimdegdist_sim)
## wilcoxon_test_and_plot(Th17degdist_real, Th17degdist_sim)
##
## t-test:
## t_test_and_plot(Th17degdist_real, Th17degdist_sim)
##
## CENTRALITIES:
## wilcoxon_test_and_plot(Th17degcent_real, Th17degcent_sim)
## t_test_and_plot(Th17degcent_real, Th17degcent_sim)
##
## One-sample hypothesis test on eigenvector centrality:
## perform_hypothesis_test(
##   gl_centralities      = Th17centralities_sim,
##   g_centralities_real  = Th17centralities_real,
##   p_value_threshold    = 0.05
## )
##
## NOTE: no write.csv() here; if you want to export summaries, add them
## explicitly in a separate block and keep them commented in the repo.
