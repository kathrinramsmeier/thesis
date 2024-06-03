
sessions <- c(
  "C190127"#, "C190128", 
  # "C190201",
  #"C190203", "C190204", "C190206"#, "C190207",
  # "C190208",
  # "C190213", "C190218", "C190219", "C190220", "C190222", "C190223",
  # "C190225", "C190307", "C190313", "C190320", "C190322", "C190327", "C190416",
  # "H190608", "H190612", "H190620", "H190621", "H190625", "H190626", "H190627", "H190629"
)

# frequency bands
# https://ncbi.nlm.nih.gov/pmc/articles/PMC3122299/#:~:text=Neural%20oscillations%20are%20electrical%20activities,gamma%20(%3E80%20Hz).
frequency_bands <- list(
  delta = c(0, 4),
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
  
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/power/unprimed")
  mean_bands_power_unprimed_all_sessions[i, , ] <- readRDS(paste0(session, "_mean_bands_power_unprimed.Rds"))
  
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/power/primed")
  mean_bands_power_primed_all_sessions[i, , ] <- readRDS(paste0(session, "_mean_bands_power_primed.Rds"))
  
}

# calculate the average mean bands power over all sessions
mean_bands_power_unprimed_all_sessions_avg <- apply(mean_bands_power_unprimed_all_sessions, c(2, 3), sum) / num_sessions
mean_bands_power_primed_all_sessions_avg <- apply(mean_bands_power_primed_all_sessions, c(2, 3), sum) / num_sessions

rm(list = c("mean_bands_power_unprimed_all_sessions", "mean_bands_power_primed_all_sessions"))

colnames(mean_bands_power_unprimed_all_sessions_avg) <- names(frequency_bands)
colnames(mean_bands_power_primed_all_sessions_avg) <- names(frequency_bands)
rownames(mean_bands_power_unprimed_all_sessions_avg) <- 1:num_channels
rownames(mean_bands_power_primed_all_sessions_avg) <- 1:num_channels

plot_power_heatmaps(
  mean_bands_power_unprimed = mean_bands_power_unprimed_all_sessions_avg, 
  mean_bands_power_primed = mean_bands_power_primed_all_sessions_avg, 
  un_primed = "unprimed"
)

plot_power_heatmaps(
  mean_bands_power_unprimed = mean_bands_power_unprimed_all_sessions_avg, 
  mean_bands_power_primed = mean_bands_power_primed_all_sessions_avg, 
  un_primed = "primed"
)

plot_power(
  mean_bands_power_unprimed = mean_bands_power_unprimed_all_sessions_avg, 
  mean_bands_power_primed = mean_bands_power_primed_all_sessions_avg,
  frequency_bands = frequency_bands
) 

plot_power_primed
plot_power_unprimed
plot_power_unprimed_primed # primed is dark blue, unprimed light blue

setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/combined_results/power")
saveRDS(mean_bands_power_unprimed_all_sessions_avg, "mean_bands_power_unprimed_all_sessions_avg.Rds")
saveRDS(mean_bands_power_primed_all_sessions_avg, "mean_bands_power_primed_all_sessions_avg.Rds")



# Coherence  --------------------------------------------

coherence_hat_unprimed_all_sessions <- coherence_hat_primed_all_sessions <- array(dim = c(num_sessions, num_channels, num_channels, num_frequency_bands))

for (i in 1:num_sessions) {
  
  session <- sessions[i]
  
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/coherence/unprimed")
  coherence_hat_unprimed_all_sessions[i, , , ] <- readRDS(paste0(session, "_avg_coherence_unprimed.Rds"))
  
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/coherence/primed")
  coherence_hat_primed_all_sessions[i, , , ] <- readRDS(paste0(session, "_avg_coherence_primed.Rds"))
  
}

# calculate the average over all sessions
coherence_hat_unprimed_all_sessions_avg <- apply(coherence_hat_unprimed_all_sessions, c(2, 3, 4), sum) / num_sessions
coherence_hat_primed_all_sessions_avg <- apply(coherence_hat_primed_all_sessions, c(2, 3, 4), sum) / num_sessions

rm(list = c("coherence_hat_unprimed_all_sessions", "coherence_hat_primed_all_sessions"))

plot_heatmap(
  result_array = coherence_hat_unprimed_all_sessions_avg, 
  frequency_bands = frequency_bands, 
  un_primed = "Unprimed"
)
plot_heatmap(
  result_array = coherence_hat_primed_all_sessions_avg, 
  frequency_bands = frequency_bands, 
  un_primed = "Primed"
)
plot_diff_heatmap(
  result_array_unprimed = coherence_hat_unprimed_all_sessions_avg,
  results_array_primed =  coherence_hat_primed_all_sessions_avg, 
  frequency_bands = frequency_bands
)

setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/combined_results/coherence")
saveRDS(coherence_hat_unprimed_all_sessions_avg, "coherence_hat_unprimed_all_sessions_avg.Rds")
saveRDS(coherence_hat_primed_all_sessions_avg, "coherence_hat_primed_all_sessions_avg.Rds")



# Phase Analysis ----------------------------------------------------------

# PLV_hat_unprimed_all_sessions <- 
  # PLV_hat_primed_all_sessions <- 
  PLI_hat_unprimed_all_sessions <- 
  PLI_hat_primed_all_sessions <- array(dim = c(num_sessions, num_channels, num_channels, num_frequency_bands))

for (i in 1:num_sessions) {
  
  session <- sessions[i]
  
  # setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/PLV/unprimed")
  # PLV_hat_unprimed_all_sessions[i, , , ] <- readRDS(paste0(session, "_PLV_unprimed.Rds"))
  
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/PLI/unprimed")
  PLI_hat_unprimed_all_sessions[i, , , ] <- readRDS(paste0(session, "_PLI_unprimed.Rds"))
  
  # setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/PLV/primed")
  # PLV_hat_primed_all_sessions[i, , , ] <- readRDS(paste0(session, "_PLV_primed.Rds"))
  
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/PLI/primed")
  PLI_hat_primed_all_sessions[i, , , ] <- readRDS(paste0(session, "_PLI_primed.Rds"))
  
}

# calculate the average over all sessions
PLV_hat_unprimed_all_sessions_avg <- apply(PLV_hat_unprimed_all_sessions, c(2, 3, 4), sum) / num_sessions
PLV_hat_primed_all_sessions_avg <- apply(PLV_hat_primed_all_sessions, c(2, 3, 4), sum) / num_sessions
PLI_hat_unprimed_all_sessions_avg <- apply(PLI_hat_unprimed_all_sessions, c(2, 3, 4), sum) / num_sessions
PLI_hat_primed_all_sessions_avg <- apply(PLI_hat_primed_all_sessions, c(2, 3, 4), sum) / num_sessions

rm(list = c(
  "coherence_hat_unprimed_all_sessions", "coherence_hat_primed_all_sessions",
  "PLV_hat_unprimed_all_sessions", "PLV_hat_primed_all_sessions",
  "PLI_hat_unprimed_all_sessions", "PLI_hat_primed_all_sessions"
))

plot_heatmap(
  result_array = PLV_hat_unprimed_all_sessions_avg, 
  frequency_bands = frequency_bands, 
  un_primed = "Unprimed"
)
plot_heatmap(
  result_array = PLV_hat_primed_all_sessions_avg, 
  frequency_bands = frequency_bands, 
  un_primed = "Primed"
)
plot_heatmap(
  result_array = PLI_hat_unprimed_all_sessions_avg, 
  frequency_bands = frequency_bands, 
  un_primed = "Unprimed"
)
plot_heatmap(
  result_array = PLI_hat_primed_all_sessions_avg, 
  frequency_bands = frequency_bands, 
  un_primed = "Primed"
)

# Granger Causality -------------------------------------------------------

# GC_p_values_unprimed_all_sessions <- array(dim = c(num_sessions, , ))
# GC_p_values_primed_all_sessions <- array(dim = c(num_sessions, , ))
GC_percentage_sign_layer_levels_unprimed_all_sessions <- matrix(nrow = 9, ncol = num_sessions)
GC_percentage_sign_layer_levels_primed_all_sessions <- matrix(nrow = 9, ncol = num_sessions)

for (i in 1:num_sessions) {
  
  session <- sessions[i]
  
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/granger_causality/unprimed")
  # GC_p_values_unprimed_all_sessions[i, , ] <- readRDS(paste0(session, "_GC_p_values_unprimed.Rds"))
  GC_percentage_sign_layer_levels_unprimed_all_sessions[, i] <- as.numeric(readRDS(paste0(session, "_GC_percentage_sign_layer_levels_unprimed.Rds"))[, 3])
  
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/granger_causality/primed")
  # GC_p_values_primed_all_sessions[i, , ] <- readRDS(paste0(session, "_GC_p_values_primed.Rds"))
  GC_percentage_sign_layer_levels_primed_all_sessions[, i] <- as.numeric(readRDS(paste0(session, "_GC_percentage_sign_layer_levels_primed.Rds"))[, 3])
  
}

# calculate the average over all sessions
GC_percentage_sign_layer_levels_unprimed_all_sessions_avg <- apply(GC_percentage_sign_layer_levels_unprimed_all_sessions, 1, mean)
GC_percentage_sign_layer_levels_primed_all_sessions_avg <- apply(GC_percentage_sign_layer_levels_primed_all_sessions, 1, mean)

rm(list = c("GC_percentage_sign_layer_levels_unprimed_all_sessions", "GC_percentage_sign_layer_levels_primed_all_sessions"))

GC_percentage_sign_layer_levels_unprimed_all_sessions_avg <- cbind(
  c("upper", "upper", "upper", "middle", "middle", "middle", "deep", "deep", "deep"),
  c("upper", "middle", "deep", "upper", "middle", "deep", "upper", "middle", "deep"),
  GC_percentage_sign_layer_levels_unprimed_all_sessions_avg
)
GC_percentage_sign_layer_levels_primed_all_sessions_avg <- cbind(
  c("upper", "upper", "upper", "middle", "middle", "middle", "deep", "deep", "deep"),
  c("upper", "middle", "deep", "upper", "middle", "deep", "upper", "middle", "deep"),
  GC_percentage_sign_layer_levels_primed_all_sessions_avg
)
colnames(GC_percentage_sign_layer_levels_unprimed_all_sessions_avg) <- colnames(GC_percentage_sign_layer_levels_primed_all_sessions_avg) <- 
  c("influencing_channel_level", "influenced_channel_level", "percentage_sign")

plot_layer_influences(
  percentage_sign = GC_percentage_sign_layer_levels_unprimed_all_sessions_avg,
  un_primed = "Unprimed"
)

plot_layer_influences(
  percentage_sign = GC_percentage_sign_layer_levels_primed_all_sessions_avg,
  un_primed = "Primed"
)

setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/combined_results/GC")
saveRDS(GC_percentage_sign_layer_levels_unprimed_all_sessions_avg, "GC_percentage_sign_layer_levels_unprimed_all_sessions_avg.Rds")
saveRDS(GC_percentage_sign_layer_levels_primed_all_sessions_avg, "GC_percentage_sign_layer_levels_primed_all_sessions_avg.Rds")
