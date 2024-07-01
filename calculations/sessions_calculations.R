
sessions <- c(
  "C190127", "C190128", "C190203", "C190204", 
  "C190206", "C190207", "C190213", "C190218", 
  "C190219", "C190220", "C190222", "C190223",
  "C190225", "C190307", "C190313", "C190320", 
  "C190322", "C190327", "C190416", "H190608", 
  "H190612", "H190620", "H190621", "H190625", 
  "H190626", "H190627", "H190629"
)

# session <- "C190127"

# frequency bands
frequency_bands <- list(
  theta = c(4, 8),
  alpha = c(8, 12),
  beta = c(12, 30),
  gamma = c(30, 200)
)

sampling_rate <- 1017.253
nyquist_freq <- sampling_rate / 2

upper_channels <- 1:5
middle_channels <- 6:10
deep_channels <- 11:15



# Data Preprocessing ------------------------------------------------------

session_data <- load_data(session = session) # if error do: detach("package:reticulate")
session_data <- data_preprocessing(
  session_data = session_data, 
  same_lengths = FALSE,
  ms_clipping = 10
)

# get the data sets out of the list into the global environment
list2env(session_data, envir = .GlobalEnv)
rm(list = "session_data")



# Stationarity ------------------------------------------------------------

# check stationarity
adf_result_baseline <- adf_stationarity_test(LFP = LFP)
adf_result_baseline

# find best filter order for notch filtering
LFP_stationary <- transform_stationary(
  LFP = LFP, 
  detrending = TRUE, 
  normalising = TRUE, 
  ensemble_adjusting = FALSE, 
  notch_filtering = FALSE,
  notch_filter_order = notch_filter_order,
  sampling_rate = sampling_rate,
  print_remaining_non_stationary = FALSE
)
notch_filter_order <- find_best_filter_order(
  LFP = LFP_stationary, 
  max_filter_order = 20, 
  nyquist_freq = nyquist_freq
)

# save best filter orders
setwd("~/data_results/notch_filter_order")
saveRDS(notch_filter_order, paste0(session, "_notch_filter_order.Rds"))

# detrending, normalising and notch filtering
LFP_stationary <- transform_stationary(
  LFP = LFP, 
  detrending = TRUE, 
  normalising = TRUE, 
  ensemble_adjusting = FALSE, 
  notch_filtering = TRUE,
  notch_filter_order = notch_filter_order,
  sampling_rate = sampling_rate,
  print_remaining_non_stationary = FALSE
)

# check stationarity
adf_result_notch_filtered <- adf_stationarity_test(LFP = LFP_stationary)

# filter all primed and unprimed trials (indices)
unprimed_ind <- task$trial_number_count[task$block_trial_count %in% c(1, 2)]
primed_ind <- task$trial_number_count[!(task$block_trial_count %in% c(1, 2))]

# save stationary LFPs
setwd("~/data_results/stationarity")
saveRDS(LFP_stationary, paste0(session, "_LFP_stationary.Rds"))

# save (un)primed indices
un_primed_indices <- list(
  "unprimed_ind" = unprimed_ind,
  "primed_ind" = primed_ind,
  "unprimed_ind_hra" = unprimed_ind_hra,
  "primed_ind_hra" = primed_ind_hra
)
setwd("~/data_results/stationarity")
saveRDS(un_primed_indices, paste0(session, "_un_primed_indices.Rds"))



# Power Analysis ----------------------------------------------------------

for (i in 1:length(sessions)) {
  
  session <- sessions[i]
  
  # load stationary LFPs
  setwd("~/data_results/stationarity")
  LFP_stationary <- readRDS(paste0(session, "_LFP_stationary.Rds"))
  
  # load (un)primed indices
  setwd("~/data_results/stationarity")
  un_primed_indices <- readRDS(paste0(session, "_un_primed_indices.Rds"))
  unprimed_ind <- un_primed_indices$unprimed_ind
  primed_ind <- un_primed_indices$primed_ind
  
  # mean band power in every channel for unprimed trials
  mean_bands_power_unprimed <- calculate_mean_bands_power(
    LFP = LFP_stationary, 
    un_primed_ind = unprimed_ind,
    hanning_windowing = TRUE,
    frequency_bands = frequency_bands
  )
  
  # mean band power in every channel for primed trials
  mean_bands_power_primed <- calculate_mean_bands_power(
    LFP = LFP_stationary, 
    un_primed_ind = primed_ind,
    hanning_windowing = TRUE,
    frequency_bands = frequency_bands
  )
  
  # save power analysis results
  setwd("~/data_results/power_hanning_windowed/unprimed")
  saveRDS(mean_bands_power_unprimed, paste0(session, "_mean_bands_power_unprimed.Rds"))
  setwd("~/data_results/power_hanning_windowed/primed")
  saveRDS(mean_bands_power_primed, paste0(session, "_mean_bands_power_primed.Rds"))
  
  print(paste("Session", i, "of", length(sessions), "completed."))
  
}



# Phase Analysis ----------------------------------------------------------


for (i in 1:length(sessions)) {
  
  session <- sessions[i]
  
  # load stationary LFPs
  setwd("~/data_results/stationarity")
  LFP_stationary <- readRDS(paste0(session, "_LFP_stationary.Rds"))
  
  # load (un)primed indices
  setwd("~/data_results/stationarity")
  un_primed_indices <- readRDS(paste0(session, "_un_primed_indices.Rds"))
  unprimed_ind <- un_primed_indices$unprimed_ind
  primed_ind <- un_primed_indices$primed_ind
  
  # calculate coherence_hat between the channels for all frequency bands for unprimed trials
  coherence_hat_unprimed <- calculate_avg_coherence(
    LFP = LFP_stationary,
    un_primed_ind = unprimed_ind,
    frequency_bands = frequency_bands,
    sampling_rate = sampling_rate, 
    filter_order = 2
  )
  
  setwd("~/data_results/coherence/unprimed")
  saveRDS(coherence_hat_unprimed, paste0(session, "_avg_coherence_unprimed.Rds"))
  
  # calculate coherence_hat between the channels for all frequency bands for primed trials
  coherence_hat_primed <- calculate_avg_coherence(
    LFP = LFP_stationary,
    un_primed_ind = primed_ind,
    frequency_bands = frequency_bands,
    sampling_rate = sampling_rate, 
    filter_order = 2
  )
  
  setwd("~/data_results/coherence/primed")
  saveRDS(coherence_hat_primed, paste0(session, "_avg_coherence_primed.Rds"))
  
  # calculate PLV between the channels for all frequency bands for unprimed trials
  PLV_unprimed <- calculate_PLV_PLI_PPC_hat(
    LFP = LFP_stationary,
    un_primed_ind = unprimed_ind,
    frequency_bands = frequency_bands,
    sampling_rate = sampling_rate, 
    filter_order = 2,
    method = "PLV"
  )
  
  setwd("~/data_results/PLV/unprimed")
  saveRDS(PLV_unprimed, paste0(session, "_PLV_unprimed_unprimed.Rds"))
  
  # calculate PLV between the channels for all frequency bands for primed trials
  PLV_primed <- calculate_PLV_PLI_PPC_hat(
    LFP = LFP_stationary,
    un_primed_ind = primed_ind,
    frequency_bands = frequency_bands,
    sampling_rate = sampling_rate, 
    filter_order = 2,
    method = "PLV"
  )
  
  setwd("~/data_results/PLV/primed")
  saveRDS(PLV_primed, paste0(session, "_PLV_unprimed_primed.Rds"))
  
  # calculate PLI between the channels for all frequency bands for unprimed trials
  PLI_unprimed <- calculate_PLV_PLI_PPC_hat(
    LFP = LFP_stationary,
    un_primed_ind = unprimed_ind,
    frequency_bands = frequency_bands,
    sampling_rate = sampling_rate, 
    filter_order = 2,
    method = "PLI"
  )
  
  setwd("~/data_results/PLI/unprimed")
  saveRDS(PLI_unprimed, paste0(session, "_PLI_unprimed.Rds"))
  
  # calculate PLI between the channels for all frequency bands for primed trials
  PLI_primed <- calculate_PLV_PLI_PPC_hat(
    LFP = LFP_stationary,
    un_primed_ind = primed_ind,
    frequency_bands = frequency_bands,
    sampling_rate = sampling_rate, 
    filter_order = 2,
    method = "PLI"
  )
  
  setwd("~/data_results/PLI/primed")
  saveRDS(PLI_primed, paste0(session, "_PLI_primed.Rds"))
  
  # calculate PPC between the channels for all frequency bands for unprimed trials
  PPC_unprimed <- calculate_PLV_PLI_PPC_hat(
    LFP = LFP_stationary,
    un_primed_ind = unprimed_ind,
    frequency_bands = frequency_bands[2:5],
    sampling_rate = sampling_rate, 
    filter_order = 2,
    method = "PPC"
  )
  
  setwd("~/data_results/PPC/unprimed")
  saveRDS(PPC_unprimed, paste0(session, "_PPC_unprimed.Rds"))
  
  # calculate PPC between the channels for all frequency bands for primed trials
  PPC_primed <- calculate_PLV_PLI_PPC_hat(
    LFP = LFP_stationary,
    un_primed_ind = primed_ind,
    frequency_bands = frequency_bands[2:5],
    sampling_rate = sampling_rate, 
    filter_order = 2,
    method = "PPC"
  )
  
  setwd("~/data_results/PPC/primed")
  saveRDS(PPC_primed, paste0(session, "_PPC_primed.Rds"))
  
  print(paste("Session", i, "of", length(sessions), "Sessions completed."))
  
}



# Granger Causality -------------------------------------------------------

for (i in 1:length(sessions)) {
  
  session <- sessions[i]
  
  # load stationary LFPs
  setwd("~/data_results/stationarity")
  LFP_stationary_hra <- readRDS(paste0(session, "_LFP_stationary.Rds"))
  
  # load (un)primed indices
  setwd("~/data_results/stationarity")
  un_primed_indices <- readRDS(paste0(session, "_un_primed_indices.Rds"))
  unprimed_ind_hra <- un_primed_indices$unprimed_ind
  primed_ind_hra <- un_primed_indices$primed_ind
  
  # GC for unprimed and primed trials
  GC_p_values_unprimed <- calculate_GC_p_values(
    LFP = LFP_stationary_hra,
    un_primed_ind = unprimed_ind_hra,
    lag_order = 2
  )
  GC_percentage_sign_layer_levels_unprimed <- transform_layer_levels(GC_p_values = GC_p_values_unprimed)
  GC_p_values_primed <- calculate_GC_p_values(
    LFP = LFP_stationary_hra,
    un_primed_ind = primed_ind_hra,
    lag_order = 2
  )
  GC_percentage_sign_layer_levels_primed <- transform_layer_levels(GC_p_values = GC_p_values_primed)
  
  # save GC results
  setwd("~/data_results/granger_causality/unprimed")
  saveRDS(GC_p_values_unprimed, paste0(session, "_GC_p_values_unprimed.Rds"))
  saveRDS(GC_percentage_sign_layer_levels_unprimed, paste0(session, "_GC_percentage_sign_layer_levels_unprimed.Rds"))
  setwd("~/data_results/granger_causality/primed")
  saveRDS(GC_p_values_primed, paste0(session, "_GC_p_values_primed.Rds"))
  saveRDS(GC_percentage_sign_layer_levels_primed, paste0(session, "_GC_percentage_sign_layer_levels_primed.Rds"))
  
}

# checks
any(is.na(GC_p_values_unprimed$p_values))
mean(GC_p_values_unprimed$adj_sum_square_err < 0.3)
mean(GC_p_values_unprimed$consistency < 0.8)
GC_p_values_unprimed$portmanteu_test_res
GC_p_values_unprimed$portmanteu_test_res
plot_residuals(GC_p_values_unprimed$res[[38]])

# checks
any(is.na(GC_p_values_primed$p_values))
mean(GC_p_values_primed$adj_sum_square_err < 0.3)
mean(GC_p_values_primed$consistency < 0.8)
GC_p_values_primed$portmanteu_test_res
GC_p_values_primed$portmanteu_test_res
plot_residuals(GC_p_values_primed$res[[30]])
