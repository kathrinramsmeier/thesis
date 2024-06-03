
# delete data from previous session to get storage space
rm(list = c(
  "LFP", "probe", "time", "behaviour", "task", "recordinginfo", "probe",
  "LFP_stationary", "LFP_stationary_hra"
))

session <- "C190127"

sessions <- c(
  # "C190127"#, 
  # "C190128", 
  # # "C190201",
  # "C190203", "C190204", 
  # "C190206", 
  "C190207",
  # # "C190208",
  "C190213", "C190218", "C190219", "C190220", "C190222", "C190223",
  "C190225", "C190307", "C190313", "C190320", "C190322", "C190327", "C190416",
  "H190608", "H190612", "H190620", "H190621", "H190625", "H190626", "H190627", "H190629"
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

# sampling rate
sampling_rate <- recordinginfo$data_sampling_rate_hz

# Nyquist frequency: maximum frequency that can be represented in the digital signal
nyquist_freq <- sampling_rate / 2

# get channel information
upper_channels <- probe$contact[probe$layer_category == "upper"]
middle_channels <- probe$contact[probe$layer_category == "middle"]
middle_channels[c(1, 3)] <- c(6, 8)
deep_channels <- probe$contact[probe$layer_category == "deep"]



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
setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/notch_filter_order")
saveRDS(notch_filter_order, paste0(session, "_notch_filter_order.Rds"))

# # load best filter orders
# setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/notch_filter_order")
# notch_filter_order <- readRDS(paste0(session, "_notch_filter_order.Rds"))

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
adf_result_notch_filtered

# remove remaining non-stationary trials (for GC analysis)
LFP_stationary_hra <- exclude_remaining_non_stationary(LFP = LFP_stationary)
length(LFP_stationary) - length(LFP_stationary_hra)

# filter all primed and unprimed trials (indices)
unprimed_ind <- task$trial_number_count[task$block_trial_count %in% c(1, 2)]
primed_ind <- task$trial_number_count[!(task$block_trial_count %in% c(1, 2))]

# adapt unprimed and primed ind to removed trials
unprimed_ind_hra <- unprimed_ind[which(unprimed_ind %in% c(1:length(LFP_stationary_hra)))]
primed_ind_hra <- primed_ind[which(primed_ind %in% c(1:length(LFP_stationary_hra)))]
length(unprimed_ind_hra)
length(primed_ind_hra)

# save stationary LFPs
setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/stationarity")
saveRDS(LFP_stationary, paste0(session, "_LFP_stationary.Rds"))
saveRDS(LFP_stationary_hra, paste0(session, "_LFP_stationary_hra.Rds"))

# # load stationary LFPs
# setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/stationarity")
# LFP_stationary <- readRDS(paste0(session, "_LFP_stationary.Rds"))
# LFP_stationary_hra <- readRDS(paste0(session, "_LFP_stationary_hra.Rds"))

# save (un)primed indices
un_primed_indices <- list(
  "unprimed_ind" = unprimed_ind,
  "primed_ind" = primed_ind,
  "unprimed_ind_hra" = unprimed_ind_hra,
  "primed_ind_hra" = primed_ind_hra
)
setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/stationarity")
saveRDS(un_primed_indices, paste0(session, "_un_primed_indices.Rds"))

# # load (un)primed indices
# setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/stationarity")
# un_primed_indices <- readRDS(paste0(session, "_un_primed_indices.Rds"))
# unprimed_ind <- un_primed_indices$unprimed_ind
# primed_ind <- un_primed_indices$primed_ind
# unprimed_ind_hra <- un_primed_indices$unprimed_ind_hra
# primed_ind_hra <- un_primed_indices$primed_ind_hra



# Power Analysis ----------------------------------------------------------

for (i in 1:length(sessions)) {
  
  session <- sessions[i]
  
  # load stationary LFPs
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/stationarity")
  LFP_stationary <- readRDS(paste0(session, "_LFP_stationary.Rds"))
  
  # load (un)primed indices
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/stationarity")
  un_primed_indices <- readRDS(paste0(session, "_un_primed_indices.Rds"))
  unprimed_ind <- un_primed_indices$unprimed_ind
  primed_ind <- un_primed_indices$primed_ind
  
  sampling_rate <- 1017.253
  
  # Nyquist frequency: maximum frequency that can be represented in the digital signal
  nyquist_freq <- sampling_rate / 2
  
  # mean band power in every channel for unprimed trials
  mean_bands_power_unprimed <- calculate_mean_bands_power(
    LFP = LFP_stationary, 
    un_primed_ind = unprimed_ind,
    hanning_windowing = FALSE,
    frequency_bands = frequency_bands
  )
  
  # mean band power in every channel for primed trials
  mean_bands_power_primed <- calculate_mean_bands_power(
    LFP = LFP_stationary, 
    un_primed_ind = primed_ind,
    hanning_windowing = FALSE,
    frequency_bands = frequency_bands
  )
  
  # save power analysis results
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/power/unprimed")
  saveRDS(mean_bands_power_unprimed, paste0(session, "_mean_bands_power_unprimed.Rds"))
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/power/primed")
  saveRDS(mean_bands_power_primed, paste0(session, "_mean_bands_power_primed.Rds"))
  
  print(paste("Session", i, "of", length(sessions), "completed."))
  
}



# Coherence Analysis --------------------------------------------------------

for (i in 1:length(sessions)) {
  
  session <- sessions[i]
  
  # load stationary LFPs
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/stationarity")
  LFP_stationary <- readRDS(paste0(session, "_LFP_stationary.Rds"))
  
  # load (un)primed indices
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/stationarity")
  un_primed_indices <- readRDS(paste0(session, "_un_primed_indices.Rds"))
  unprimed_ind <- un_primed_indices$unprimed_ind
  primed_ind <- un_primed_indices$primed_ind
  
  sampling_rate <- 1017.253
  
  # calculate coherence_hat between the channels for all frequency bands for unprimed trials
  coherence_hat_unprimed <- calculate_avg_coherence(
    LFP = LFP_stationary,
    un_primed_ind = unprimed_ind,
    frequency_bands = frequency_bands,
    sampling_rate = sampling_rate, 
    filter_order = 2
  )
  
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/coherence/unprimed")
  saveRDS(coherence_hat_unprimed, paste0(session, "_avg_coherence_unprimed.Rds"))
  
  # calculate coherence_hat between the channels for all frequency bands for primed trials
  coherence_hat_primed <- calculate_avg_coherence(
    LFP = LFP_stationary,
    un_primed_ind = primed_ind,
    frequency_bands = frequency_bands,
    sampling_rate = sampling_rate, 
    filter_order = 2
  )
  
  setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/coherence/primed")
  saveRDS(coherence_hat_primed, paste0(session, "_avg_coherence_primed.Rds"))
  
  print(paste("Session", i, "of", length(sessions), "Sessions completed."))
  
}



# Phase Analysis ----------------------------------------------------------

# calculate PLV between the channels for all frequency bands for unprimed trials
PLV_unprimed <- calculate_PLV_PLI_PPC_hat(
  LFP = LFP_stationary,
  un_primed_ind = unprimed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  method = "PLV"
)

setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/PLV/unprimed")
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

setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/PLV/primed")
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

setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/PLI/unprimed")
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

setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/PLI/primed")
saveRDS(PLI_primed, paste0(session, "_PLI_primed.Rds"))

# calculate PPC between the channels for all frequency bands for unprimed trials
PPC_unprimed <- calculate_PLV_PLI_PPC_hat(
  LFP = LFP_stationary,
  un_primed_ind = unprimed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  method = "PPC"
)

setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/PPC/unprimed")
saveRDS(PPC_unprimed, paste0(session, "_PPC_unprimed_unprimed.Rds"))

# calculate PPC between the channels for all frequency bands for primed trials
PPC_primed <- calculate_PLV_PLI_PPC_hat(
  LFP = LFP_stationary,
  un_primed_ind = primed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  method = "PPC"
)

setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/PPC/primed")
saveRDS(PPC_primed, paste0(session, "_PPC_unprimed_primed.Rds"))



# Granger Causality -------------------------------------------------------

# load stationary LFPs
setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/stationarity")
LFP_stationary_hra <- readRDS(paste0(session, "_LFP_stationary_hra.Rds"))

# load (un)primed indices
setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/stationarity")
un_primed_indices <- readRDS(paste0(session, "_un_primed_indices.Rds"))
unprimed_ind_hra <- un_primed_indices$unprimed_ind_hra
primed_ind_hra <- un_primed_indices$primed_ind_hra

GC_p_values_unprimed <- calculate_GC_p_values(
  LFP = LFP_stationary_hra,
  un_primed_ind = unprimed_ind_hra,
  lag_order = 4
)

# checks
any(is.na(GC_p_values_unprimed$p_values))
mean(GC_p_values_unprimed$adj_sum_square_err < 0.3)
mean(GC_p_values_unprimed$consistency < 0.8)
GC_p_values_unprimed$portmanteu_test_res
GC_p_values_unprimed$portmanteu_test_res
plot_residuals(GC_p_values_unprimed$res[[38]])

GC_percentage_sign_layer_levels_unprimed <- transform_layer_levels(GC_p_values = GC_p_values_unprimed)

GC_p_values_primed <- calculate_GC_p_values(
  LFP = LFP_stationary_hra,
  un_primed_ind = primed_ind_hra,
  lag_order = 4
)

# checks
any(is.na(GC_p_values_primed$p_values))
mean(GC_p_values_primed$adj_sum_square_err < 0.3)
mean(GC_p_values_primed$consistency < 0.8)
GC_p_values_primed$portmanteu_test_res
GC_p_values_primed$portmanteu_test_res
plot_residuals(GC_p_values_primed$res[[30]])

GC_percentage_sign_layer_levels_primed <- transform_layer_levels(GC_p_values = GC_p_values_primed)

# save GC results
setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/granger_causality/unprimed")
saveRDS(GC_p_values_unprimed, paste0(session, "_GC_p_values_unprimed.Rds"))
saveRDS(GC_percentage_sign_layer_levels_unprimed, paste0(session, "_GC_percentage_sign_layer_levels_unprimed.Rds"))
setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/granger_causality/primed")
saveRDS(GC_p_values_primed, paste0(session, "_GC_p_values_primed.Rds"))
saveRDS(GC_percentage_sign_layer_levels_primed, paste0(session, "_GC_percentage_sign_layer_levels_primed.Rds"))

plot_layer_influences(GC_percentage_sign_layer_levels_unprimed, "Unprimed")
plot_layer_influences(GC_percentage_sign_layer_levels_primed, "Primed")
