
session <- "C190127"



# Data Preprocessing ------------------------------------------------------

# delete data from previous session to get storage space
rm(list = c(
  "LFP", "probe", "time", "behaviour", "task", "recordinginfo", "probe",
  "LFP_stationary", "LFP_stationary_hra"
))

session_data <- load_data(session = session)
session_data <- data_preprocessing(session_data = session_data, same_lengths = FALSE)

# get the data sets out of the list into the global environment
list2env(session_data, envir = .GlobalEnv)
rm(list = "session_data")



# Stationarity ------------------------------------------------------------

# find best filter order for notch filter
LFP_stationary <- transform_stationary(
  LFP = LFP, 
  detrending = TRUE, 
  normalising = TRUE, 
  ensemble_adjusting = TRUE, 
  notch_filtering = FALSE,
  notch_filter_order = notch_filter_order,
  sampling_rate = sampling_rate,
  print_remaining_non_stationary = TRUE,
  hard_remove_all = FALSE
)
notch_filter_order <- find_best_filter_order(LFP = LFP_stationary)

# save best filter orders
setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/notch_filter_order")
saveRDS(notch_filter_order, paste0(session, "_notch_filter_order.Rds"))

# # load best filter orders
# setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/notch_filter_order")
# notch_filter_order <- readRDS(paste0(session, "_notch_filter_order.Rds"))

LFP_stationary <- transform_stationary(
  LFP = LFP, 
  detrending = TRUE, 
  normalising = TRUE, 
  ensemble_adjusting = TRUE, 
  notch_filtering = TRUE,
  notch_filter_order = notch_filter_order,
  sampling_rate = sampling_rate,
  print_remaining_non_stationary = TRUE,
  hard_remove_all = FALSE
)
LFP_stationary[1:3, 1:3, 1:3]

LFP_stationary_hra <- transform_stationary(
  LFP = LFP, 
  detrending = TRUE, 
  normalising = TRUE, 
  ensemble_adjusting = TRUE, 
  notch_filtering = TRUE,
  notch_filter_order = notch_filter_order,
  sampling_rate = sampling_rate,
  print_remaining_non_stationary = FALSE,
  hard_remove_all = TRUE
)
dim(LFP_stationary_hra)

# filter all primed and unprimed trials (indices)
primed_ind <- task$trial_number_count[task$block_trial_count %in% c(1, 2)]
unprimed_ind <- task$trial_number_count[!(task$block_trial_count %in% c(1, 2))]

# adapt unprimed and primed ind to removed trials
unprimed_ind_hra <- unprimed_ind[which(unprimed_ind %in% c(1:dim(LFP_stationary_hra)[3]))]
primed_ind_hra <- primed_ind[which(primed_ind %in% c(1:dim(LFP_stationary_hra)[3]))]
length(unprimed_ind_hra)
length(primed_ind_hra)

# save stationary LFPs
setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/stationarity")
saveRDS(LFP_stationary, paste0(session, "_LFP_stationary.Rds"))
saveRDS(LFP_stationary_hra, paste0(session, "_LFP_stationary_hra.Rds"))

# load stationary LFPs
setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/stationarity")
LFP_stationary <- readRDS(paste0(session, "_LFP_stationary.Rds"))
LFP_stationary_hra <- readRDS(paste0(session, "_LFP_stationary_hra.Rds"))

# sampling rate
sampling_rate <- recordinginfo$data_sampling_rate_hz

# Nyquist frequency: maximum frequency that can be represented in the digital signal
nyquist_freq <- sampling_rate / 2

# get channel information
upper_channels <- probe$contact[probe$layer_category == "upper"]
middle_channels <- probe$contact[probe$layer_category == "middle"]
middle_channels[c(1, 3)] <- c(6, 8)
deep_channels <- probe$contact[probe$layer_category == "deep"]



# Power Analysis ----------------------------------------------------------

# mean band power in every channel for unprimed trials
mean_bands_power_hanning_windowed_unprimed <- calculate_mean_bands_power(
  LFP = LFP_stationary, 
  un_primed_ind = unprimed_ind,
  hanning_windowed = TRUE,
  frequency_bands = frequency_bands
)
mean_bands_power_unprimed <- calculate_mean_bands_power(
  LFP = LFP_stationary, 
  un_primed_ind = unprimed_ind,
  hanning_windowed = FALSE,
  frequency_bands = frequency_bands
)

# mean band power in every channel for primed trials
mean_bands_power_hanning_windowed_primed <- calculate_mean_bands_power(
  LFP = LFP_stationary, 
  un_primed_ind = primed_ind,
  hanning_windowed = TRUE,
  frequency_bands = frequency_bands
)
mean_bands_power_primed <- calculate_mean_bands_power(
  LFP = LFP_stationary, 
  un_primed_ind = primed_ind,
  hanning_windowed = FALSE,
  frequency_bands = frequency_bands
)

# check
plot_power_heatmaps(
  mean_bands_power_unprimed = mean_bands_power_hanning_windowed_unprimed, 
  mean_bands_power_primed = mean_bands_power_hanning_windowed_primed, 
  un_primed = "unprimed"
)
plot_power_heatmaps(
  mean_bands_power_unprimed = mean_bands_power_hanning_windowed_unprimed, 
  mean_bands_power_primed = mean_bands_power_hanning_windowed_primed, 
  un_primed = "primed"
)
plot_power_heatmaps(
  mean_bands_power_unprimed = mean_bands_power_unprimed, 
  mean_bands_power_primed = mean_bands_power_primed, 
  un_primed = "unprimed"
)
plot_power_heatmaps(
  mean_bands_power_unprimed = mean_bands_power_unprimed, 
  mean_bands_power_primed = mean_bands_power_primed, 
  un_primed = "primed"
)

# save power analysis results
setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/power/unprimed")
saveRDS(mean_bands_power_hanning_windowed_unprimed, paste0(session, "_mean_bands_power_hanning_windowed_unprimed.Rds"))
saveRDS(mean_bands_power_unprimed, paste0(session, "_mean_bands_power_unprimed.Rds"))
setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/power/primed")
saveRDS(mean_bands_power_hanning_windowed_primed, paste0(session, "_mean_bands_power_hanning_windowed_primed.Rds"))
saveRDS(mean_bands_power_primed, paste0(session, "_mean_bands_power_primed.Rds"))



# Phase Analysis ----------------------------------------------------------

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

# calculate PPC - PLV between the channels for all frequency bands for unprimed trials
PPC_PLV_unprimed <- calculate_PPC(
  LFP = LFP_stationary,
  un_primed_ind = unprimed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  PLV_PLI = "PLV"
)

setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/PPC_PLV/unprimed")
saveRDS(PPC_PLV_unprimed, paste0(session, "_PPC_PLV_unprimed_unprimed.Rds"))

# calculate PPC - PLV between the channels for all frequency bands for primed trials
PPC_PLV_primed <- calculate_PPC(
  LFP = LFP_stationary,
  un_primed_ind = primed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  PLV_PLI = "PLV"
)

setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/PPC_PLV/primed")
saveRDS(PPC_PLV_primed, paste0(session, "_PPC_PLV_unprimed_primed.Rds"))

# calculate PPC - PLI between the channels for all frequency bands for unprimed trials
PPC_PLI_unprimed <- calculate_PPC(
  LFP = LFP_stationary,
  un_primed_ind = unprimed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  PLV_PLI = "PLI"
)

setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/PPC_PLI/unprimed")
saveRDS(PPC_PLI_unprimed, paste0(session, "_PPC_PLI_unprimed_unprimed.Rds"))

# calculate PPC - PLI between the channels for all frequency bands for primed trials
PPC_PLI_primed <- calculate_PPC(
  LFP = LFP_stationary,
  un_primed_ind = primed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  PLV_PLI = "PLI"
)

setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/PPC_PLI/primed")
saveRDS(PPC_PLI_primed, paste0(session, "_PPC_PLI_unprimed_primed.Rds"))



# Granger Causality -------------------------------------------------------

GC_p_values_unprimed <- calculate_GC_p_values(
  LFP = LFP_stationary_hra,
  un_primed_ind = unprimed_ind_hra,
  lag_order = 5
)

# check
any(is.na(GC_p_values_unprimed))

GC_percentage_sign_layer_levels_unprimed <- transform_layer_levels(GC_p_values = GC_p_values_unprimed)

GC_p_values_primed <- calculate_GC_p_values(
  LFP = LFP_stationary_hra,
  un_primed_ind = primed_ind_hra,
  lag_order = 5
)

# check
any(is.na(GC_p_values_primed))

GC_percentage_sign_layer_levels_primed <- transform_layer_levels(GC_p_values = GC_p_values_primed)

# save GC results
setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/granger_causality/unprimed")
saveRDS(GC_p_values_unprimed, paste0(session, "_GC_p_values_unprimed.Rds"))
saveRDS(GC_percentage_sign_layer_levels_unprimed, paste0(session, "_GC_percentage_sign_layer_levels_unprimed.Rds"))
setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/granger_causality/primed")
saveRDS(GC_p_values_primed, paste0(session, "_GC_p_values_primed.Rds"))
saveRDS(GC_percentage_sign_layer_levels_primed, paste0(session, "_GC_percentage_sign_layer_levels_primed.Rds"))


