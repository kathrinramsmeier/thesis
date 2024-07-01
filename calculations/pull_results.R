
sessions <- c(
  "C190127", "C190128", "C190203", "C190204", 
  "C190206", "C190207", "C190213", "C190218", 
  "C190219", "C190220", "C190222", "C190223",
  "C190225", "C190307", "C190313", "C190320", 
  "C190322", "C190327", "C190416", "H190608", 
  "H190612", "H190620", "H190621", "H190625", 
  "H190626", "H190627", "H190629"
)

# frequency bands
frequency_bands <- list(
  theta = c(4, 8),
  alpha = c(8, 12),
  beta = c(12, 30),
  gamma = c(30, 200)
)

num_sessions <- length(sessions)
num_channels <- 15
num_frequency_bands <- length(frequency_bands)



# Power Analysis ----------------------------------------------------------

mean_bands_power_unprimed_all_sessions <- array(dim = c(num_sessions, num_channels, num_frequency_bands))
mean_bands_power_primed_all_sessions <- array(dim = c(num_sessions, num_channels, num_frequency_bands))

for (i in 1:num_sessions) {
  
  session <- sessions[i]
  
  setwd("~/data_results/power/unprimed")
  mean_bands_power_unprimed_all_sessions[i, , ] <- readRDS(paste0(session, "_mean_bands_power_unprimed.Rds"))
  
  setwd("~/data_results/power/primed")
  mean_bands_power_primed_all_sessions[i, , ] <- readRDS(paste0(session, "_mean_bands_power_primed.Rds"))
  
}

# calculate the average mean bands power over all sessions
mean_bands_power_unprimed_all_sessions_avg <- apply(mean_bands_power_unprimed_all_sessions, c(2, 3), sum) / num_sessions
mean_bands_power_primed_all_sessions_avg <- apply(mean_bands_power_primed_all_sessions, c(2, 3), sum) / num_sessions

colnames(mean_bands_power_unprimed_all_sessions_avg) <- names(frequency_bands)
colnames(mean_bands_power_primed_all_sessions_avg) <- names(frequency_bands)
rownames(mean_bands_power_unprimed_all_sessions_avg) <- 1:num_channels
rownames(mean_bands_power_primed_all_sessions_avg) <- 1:num_channels

plot_power(
  mean_bands_power_unprimed = mean_bands_power_unprimed_all_sessions_avg, 
  mean_bands_power_primed = mean_bands_power_primed_all_sessions_avg,
  frequency_bands = frequency_bands
)

# grouping the results
mean_bands_power_unprimed_all_sessions_grouped <- array(0, dim = c(27, 3, 4))
group_indices <- list(1:5, 6:10, 11:15)
for (i in 1:3) {
  mean_bands_power_unprimed_all_sessions_grouped[, i, ] <- apply(mean_bands_power_unprimed_all_sessions[, group_indices[[i]], ], c(1, 3), mean)
}
mean_bands_power_primed_all_sessions_grouped <- array(0, dim = c(27, 3, 4))
group_indices <- list(1:5, 6:10, 11:15)
for (i in 1:3) {
  mean_bands_power_primed_all_sessions_grouped[, i, ] <- apply(mean_bands_power_primed_all_sessions[, group_indices[[i]], ], c(1, 3), mean)
}
mean_bands_power_unprimed_all_sessions_avg_grouped <- apply(mean_bands_power_unprimed_all_sessions_grouped, c(2, 3), sum) / num_sessions
mean_bands_power_primed_all_sessions_avg_grouped <- apply(mean_bands_power_primed_all_sessions_grouped, c(2, 3), sum) / num_sessions
colnames(mean_bands_power_unprimed_all_sessions_avg_grouped) <- names(frequency_bands)
colnames(mean_bands_power_primed_all_sessions_avg_grouped) <- names(frequency_bands)

rm(list = c("mean_bands_power_unprimed_all_sessions", "mean_bands_power_primed_all_sessions"))

setwd("~/data_results/combined_results/power")
saveRDS(mean_bands_power_unprimed_all_sessions_avg, "mean_bands_power_unprimed_all_sessions_avg.Rds")
saveRDS(mean_bands_power_primed_all_sessions_avg, "mean_bands_power_primed_all_sessions_avg.Rds")
saveRDS(mean_bands_power_unprimed_all_sessions_avg_grouped, "mean_bands_power_unprimed_all_sessions_avg_grouped.Rds")
saveRDS(mean_bands_power_primed_all_sessions_avg_grouped, "mean_bands_power_primed_all_sessions_avg_grouped.Rds")



# Phase Analysis ----------------------------------------------------------

coherence_hat_unprimed_all_sessions <- 
  coherence_hat_primed_all_sessions <-
  PLV_hat_unprimed_all_sessions <-
  PLV_hat_primed_all_sessions <-
  PLI_hat_unprimed_all_sessions <-
  PLI_hat_primed_all_sessions <-
  PPC_hat_unprimed_all_sessions <-
  PPC_hat_primed_all_sessions <-
  array(dim = c(num_sessions, num_channels, num_channels, num_frequency_bands))

# load the results from all sessions
for (i in 1:num_sessions) {
  
  session <- sessions[i]
  
  setwd("~/data_results/coherence/unprimed")
  coherence_hat_unprimed_all_sessions[i, , , ] <- readRDS(paste0(session, "_avg_coherence_unprimed.Rds"))
  
  setwd("~/data_results/coherence/primed")
  coherence_hat_primed_all_sessions[i, , , ] <- readRDS(paste0(session, "_avg_coherence_primed.Rds"))
  
  setwd("~/data_results/PLV/unprimed")
  PLV_hat_unprimed_all_sessions[i, , , ] <- readRDS(paste0(session, "_PLV_unprimed.Rds"))

  setwd("~/data_results/PLI/unprimed")
  PLI_hat_unprimed_all_sessions[i, , , ] <- readRDS(paste0(session, "_PLI_unprimed.Rds"))
  
  setwd("~/data_results/PPC/unprimed")
  PPC_hat_unprimed_all_sessions[i, , , ] <- readRDS(paste0(session, "_PPC_unprimed.Rds"))

  setwd("~/data_results/PLV/primed")
  PLV_hat_primed_all_sessions[i, , , ] <- readRDS(paste0(session, "_PLV_primed.Rds"))

  setwd("~/data_results/PLI/primed")
  PLI_hat_primed_all_sessions[i, , , ] <- readRDS(paste0(session, "_PLI_primed.Rds"))
  
  setwd("~/data_results/PPC/primed")
  PPC_hat_primed_all_sessions[i, , , ] <- readRDS(paste0(session, "_PPC_primed.Rds"))
  
}

# calculate the average over all sessions
coherence_hat_unprimed_all_sessions_avg <- apply(coherence_hat_unprimed_all_sessions, c(2, 3, 4), sum) / num_sessions
coherence_hat_primed_all_sessions_avg <- apply(coherence_hat_primed_all_sessions, c(2, 3, 4), sum) / num_sessions
PLV_hat_unprimed_all_sessions_avg <- apply(PLV_hat_unprimed_all_sessions, c(2, 3, 4), sum) / num_sessions
PLV_hat_primed_all_sessions_avg <- apply(PLV_hat_primed_all_sessions, c(2, 3, 4), sum) / num_sessions
PLI_hat_unprimed_all_sessions_avg <- apply(PLI_hat_unprimed_all_sessions, c(2, 3, 4), sum) / num_sessions
PLI_hat_primed_all_sessions_avg <- apply(PLI_hat_primed_all_sessions, c(2, 3, 4), sum) / num_sessions
PPC_hat_unprimed_all_sessions_avg <- apply(PPC_hat_unprimed_all_sessions, c(2, 3, 4), sum) / num_sessions
PPC_hat_primed_all_sessions_avg <- apply(PPC_hat_primed_all_sessions, c(2, 3, 4), sum) / num_sessions

rm(list = c(
  "coherence_hat_unprimed_all_sessions", "coherence_hat_primed_all_sessions",
  "PLV_hat_unprimed_all_sessions", "PLV_hat_primed_all_sessions",
  "PLI_hat_unprimed_all_sessions", "PLI_hat_primed_all_sessions",
  "PPC_hat_unprimed_all_sessions", "PPC_hat_primed_all_sessions"
))

# plot the results
plot_heatmap(
  result_array = coherence_hat_unprimed_all_sessions_avg, 
  frequency_bands = frequency_bands, 
  label_name = "phase coherence"
)
plot_heatmap(
  result_array = coherence_hat_primed_all_sessions_avg, 
  frequency_bands = frequency_bands, 
  label_name = "phase coherence"
)
plot_diff_heatmap(
  result_array_unprimed = coherence_hat_unprimed_all_sessions_avg,
  results_array_primed =  coherence_hat_primed_all_sessions_avg, 
  frequency_bands = frequency_bands,
  label_name = expression(delta[coherence])
)
plot_heatmap(
  result_array = PLV_hat_unprimed_all_sessions_avg, 
  frequency_bands = frequency_bands, 
  label_name = "PLV"
)
plot_heatmap(
  result_array = PLV_hat_primed_all_sessions_avg, 
  frequency_bands = frequency_bands, 
  label_name = "PLV"
)
plot_diff_heatmap(
  result_array_unprimed = PLV_hat_unprimed_all_sessions_avg,
  results_array_primed =  PLV_hat_primed_all_sessions_avg, 
  frequency_bands = frequency_bands,
  label_name = expression(delta[PLV])
)
plot_heatmap(
  result_array = PLI_hat_unprimed_all_sessions_avg, 
  frequency_bands = frequency_bands, 
  label_name = "PLI"
)
plot_heatmap(
  result_array = PLI_hat_primed_all_sessions_avg, 
  frequency_bands = frequency_bands, 
  label_name = "PLI"
)
plot_diff_heatmap(
  result_array_unprimed = PLI_hat_unprimed_all_sessions_avg,
  results_array_primed =  PLI_hat_primed_all_sessions_avg, 
  frequency_bands = frequency_bands,
  label_name = expression(delta[PLI])
)
plot_heatmap(
  result_array = PPC_hat_unprimed_all_sessions_avg, 
  frequency_bands = frequency_bands, 
  label_name = "PPC"
)
plot_heatmap(
  result_array = PPC_hat_primed_all_sessions_avg, 
  frequency_bands = frequency_bands, 
  label_name = "PPC"
)
plot_diff_heatmap(
  result_array_unprimed = PPC_hat_unprimed_all_sessions_avg,
  results_array_primed = PPC_hat_primed_all_sessions_avg,
  frequency_bands = frequency_bands,
  label_name = expression(delta[PPC])
)

# save the results
setwd("~/data_results/combined_results/coherence")
saveRDS(coherence_hat_unprimed_all_sessions_avg, "coherence_hat_unprimed_all_sessions_avg.Rds")
saveRDS(coherence_hat_primed_all_sessions_avg, "coherence_hat_primed_all_sessions_avg.Rds")
setwd("~/data_results/combined_results/PLV")
saveRDS(PLV_hat_unprimed_all_sessions_avg, "PLV_hat_unprimed_all_sessions_avg.Rds")
saveRDS(PLV_hat_primed_all_sessions_avg, "PLV_hat_primed_all_sessions_avg.Rds")
setwd("~/data_results/combined_results/PLI")
saveRDS(PLI_hat_unprimed_all_sessions_avg, "PLI_hat_unprimed_all_sessions_avg.Rds")
saveRDS(PLI_hat_primed_all_sessions_avg, "PLI_hat_primed_all_sessions_avg.Rds")
setwd("~/data_results/combined_results/PPC")
saveRDS(PPC_hat_unprimed_all_sessions_avg, "PPC_hat_unprimed_all_sessions_avg.Rds")
saveRDS(PPC_hat_primed_all_sessions_avg, "PPC_hat_primed_all_sessions_avg.Rds")



# Granger Causality -------------------------------------------------------

GC_p_values_unprimed_all_sessions <- matrix(nrow = 210, ncol = num_sessions)
GC_p_values_primed_all_sessions <- matrix(nrow = 210, ncol = num_sessions)
GC_percentage_sign_layer_levels_unprimed_all_sessions <- matrix(nrow = 9, ncol = num_sessions)
GC_percentage_sign_layer_levels_primed_all_sessions <- matrix(nrow = 9, ncol = num_sessions)

for (i in 1:num_sessions) {
  
  session <- sessions[i]
  
  setwd("~/data_results/granger_causality/unprimed")
  
  p_values_GC <- readRDS(paste0(session, "_GC_p_values_unprimed.Rds"))$p_values
  p_values_GC <- p_values_GC[, 3:ncol(p_values_GC)]
  GC_p_values_unprimed_all_sessions[, i] <- rowSums(p_values_GC < 0.05) / ncol(p_values_GC)
  
  GC_percentage_sign_layer_levels_unprimed_all_sessions[, i] <- as.numeric(readRDS(paste0(session, "_GC_percentage_sign_layer_levels_unprimed.Rds"))[, 3])
  
  setwd("~/data_results/granger_causality/non-stationary_not_removed/primed")
  
  p_values_GC <- readRDS(paste0(session, "_GC_p_values_primed.Rds"))$p_values
  p_values_GC <- p_values_GC[, 3:ncol(p_values_GC)]
  GC_p_values_primed_all_sessions[, i] <- rowSums(p_values_GC < 0.05) / ncol(p_values_GC)
  
  GC_percentage_sign_layer_levels_primed_all_sessions[, i] <- as.numeric(readRDS(paste0(session, "_GC_percentage_sign_layer_levels_primed.Rds"))[, 3])
  
}

# calculate the average over all sessions
GC_percentage_sign_layer_levels_unprimed_all_sessions_avg <- apply(GC_percentage_sign_layer_levels_unprimed_all_sessions, 1, mean)
GC_percentage_sign_layer_levels_primed_all_sessions_avg <- apply(GC_percentage_sign_layer_levels_primed_all_sessions, 1, mean)
GC_percentage_sign_unprimed_all_sessions_avg <- apply(GC_p_values_unprimed_all_sessions, 1, mean)
GC_percentage_sign_primed_all_sessions_avg <- apply(GC_p_values_primed_all_sessions, 1, mean)

rm(list = c(
  "GC_percentage_sign_layer_levels_unprimed_all_sessions", 
  "GC_percentage_sign_layer_levels_primed_all_sessions",
  "GC_p_values_unprimed_all_sessions",
  "GC_p_values_primed_all_sessions"
))

GC_percentage_sign_layer_levels_unprimed_all_sessions_avg <- cbind(
  c("upper (1-5)", "upper (1-5)", "upper (1-5)", "middle (6-10)", "middle (6-10)", "middle (6-10)", "deep (11-15)", "deep (11-15)", "deep (11-15)"),
  c("upper (1-5)", "middle (6-10)", "deep (11-15)", "upper (1-5)", "middle (6-10)", "deep (11-15)", "upper (1-5)", "middle (6-10)", "deep (11-15)"),
  GC_percentage_sign_layer_levels_unprimed_all_sessions_avg
)
GC_percentage_sign_layer_levels_primed_all_sessions_avg <- cbind(
  c("upper (1-5)", "upper (1-5)", "upper (1-5)", "middle (6-10)", "middle (6-10)", "middle (6-10)", "deep (11-15)", "deep (11-15)", "deep (11-15)"),
  c("upper (1-5)", "middle (6-10)", "deep (11-15)", "upper (1-5)", "middle (6-10)", "deep (11-15)", "upper (1-5)", "middle (6-10)", "deep (11-15)"),
  GC_percentage_sign_layer_levels_primed_all_sessions_avg
)
colnames(GC_percentage_sign_layer_levels_unprimed_all_sessions_avg) <- colnames(GC_percentage_sign_layer_levels_primed_all_sessions_avg) <- 
  c("influencing_channel_level", "influenced_channel_level", "percentage_sign")

combinations <- matrix(nrow = 15 * 15 - 15, ncol = 2)
k <- 1
for (i in 1:15) {
  for (j in 1:15) {
    # exclude identical pairs
    if (i != j) {
      combinations[k, 1] <- i
      combinations[k, 2] <- j
      k <- k + 1
    }
  }
}
GC_percentage_sign_all_sessions_avg <- cbind(
  combinations, 
  GC_percentage_sign_unprimed_all_sessions_avg,
  GC_percentage_sign_primed_all_sessions_avg
)
colnames(GC_percentage_sign_all_sessions_avg) <- c("influencing_channel_level", "influenced_channel_level", "percentage_sign_unprimed", "percentage_sign_primed")

plot_layer_influences(
  percentage_sign = GC_percentage_sign_layer_levels_unprimed_all_sessions_avg,
  un_primed = "Unprimed"
)

plot_layer_influences(
  percentage_sign = GC_percentage_sign_layer_levels_primed_all_sessions_avg,
  un_primed = "Primed"
)

setwd("~/data_results/combined_results/GC")
saveRDS(GC_percentage_sign_layer_levels_unprimed_all_sessions_avg, "GC_percentage_sign_layer_levels_unprimed_all_sessions_avg.Rds")
saveRDS(GC_percentage_sign_layer_levels_primed_all_sessions_avg, "GC_percentage_sign_layer_levels_primed_all_sessions_avg.Rds")
saveRDS(GC_percentage_sign_all_sessions_avg, "GC_percentage_sign_all_sessions_avg.Rds")
