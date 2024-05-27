
# Calculation of Phase Locking Value (PLV), Phase Lag Index (PLI) and Pairwise Phase Consistency (PPC).
# Coherence and PLV are equivalent in the Gaussian case (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3674231/#:~:text=PLV%20can%20therefore%20be%20viewed,each%20scaling%20of%20the%20wavelet.).
# but here: individual LFP trials are not normally distributed.



# look at the 6th and 14th channel of trial 333 as an example first
trial <- 27
channel <- 6
# sample_LFP_1 <- LFP_stationary[channel, , trial]
sample_LFP_1 <- as.vector(LFP_stationary[[trial]][[channel]])
channel <- 14
# sample_LFP_2 <- LFP_stationary[channel, , trial]
sample_LFP_2 <- as.vector(LFP_stationary[[trial]][[channel]])
plot(sample_LFP_1, type = "l")
lines(sample_LFP_2, type = "l", col = "red")



# Filter the LFP in the different frequency bands -------------------------

freqband_filtering <- function(LFP_trial, frequency_bands, sampling_rate, filter_order) {
  
  LFP_freqband_filtered <- matrix(nrow = length(LFP_trial), ncol = length(frequency_bands))
  colnames(LFP_freqband_filtered) <- names(frequency_bands)
  
  nyquist_freq <- sampling_rate / 2
  
  # avoid numerical problem (only needed when all frequency bands are examined)
  # frequency_bands[[1]][1] <- 0.01
  
  # loop over all frequency bands
  for (i in 1:length(frequency_bands)) {
    
    freq_band <- c(frequency_bands[[i]][1], frequency_bands[[i]][2])
    
    # butterworth filtering
    butterworth <- gsignal::butter(n = filter_order, w = freq_band / nyquist_freq, type = "pass")
    
    # forward and reverse filtering (correction for phase distortion)
    LFP_freqband_filtered[, i] <- gsignal::filtfilt(filt = butterworth$b, a = butterworth$a, x = LFP_trial)
    
  }
  
  return(LFP_freqband_filtered)
  
}

# try on a sample LFP signal
sample_LFP_1_freqband_filtered <- freqband_filtering(
  LFP_trial = sample_LFP_1, 
  frequency_bands = frequency_bands, 
  sampling_rate = sampling_rate,
  filter_order = 2
) 
sample_LFP_2_freqband_filtered <- freqband_filtering(
  LFP_trial = sample_LFP_2, 
  frequency_bands = frequency_bands, 
  sampling_rate = sampling_rate,
  filter_order = 2
)



# PLV, PLI and PPC --------------------------------------------------------

# compute PLV and PLI across the frequency bands
# formulas from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3674231/#:~:text=PLV%20can%20therefore%20be%20viewed,each%20scaling%20of%20the%20wavelet.
calculate_PLV_PLI_PPC <- function(LFP_trial_1_freqband_filtered, LFP_trial_2_freqband_filtered, method) {
  
  stopifnot(method %in% c("PLV", "PLI", "PPC"))
  
  result <- rep(NA, dim(LFP_trial_1_freqband_filtered)[2])
  names(result) <- colnames(LFP_trial_1_freqband_filtered)
  
  # loop over all frequency bands
  for (i in 1:dim(LFP_trial_1_freqband_filtered)[2]) {
    
    # Hilbert transform to get the analytical signal
    LFP_1_freqband_ht <- gsignal::hilbert(LFP_trial_1_freqband_filtered[, i])
    LFP_2_freqband_ht <- gsignal::hilbert(LFP_trial_2_freqband_filtered[, i])
    
    # compute relative phase = phase difference
    phase_diffs <- pracma::angle(LFP_1_freqband_ht * Conj(LFP_2_freqband_ht))  # arg((LFP_1_band_ht * Conj(LFP_2_band_ht)) / (abs(LFP_1_band_ht) * abs(LFP_2_band_ht)))
    
    if (method == "PLV") {
      
      # compute the PLV for each frequency band
      result[i] <- abs(mean(exp(1i * phase_diffs)))
      
    } else if (method == "PLI") {
      
      # compute the PLI for each frequency band
      result[i] <- abs(mean(sign(phase_diffs)))
      
    } else if (method == "PPC") { # check this again!!!!!

      N <- length(phase_diffs)
      sum_of_dot_products <- 0
      
      # compute pairwise dot products
      for (j in 1:(N - 1)) {
        for (k in (j + 1):N) {
          dot_product <- cos(phase_diffs)[j] * cos(phase_diffs)[k] + sin(phase_diffs)[j] * sin(phase_diffs)[k]
          sum_of_dot_products <- sum_of_dot_products + dot_product
        }
      }
      
      result[i] <- (2 / (N * (N - 1))) * sum_of_dot_products
      
    }
  }
  
  return(result)
  
}

# try PLV calc. for a sample trial for two different channels
freqband_PLV_sample_LFPs <- calculate_PLV_PLI_PPC(
  LFP_trial_1_freqband_filtered = sample_LFP_1_freqband_filtered, 
  LFP_trial_2_freqband_filtered = sample_LFP_2_freqband_filtered,
  method = "PLV"
)
freqband_PLV_sample_LFPs

# try PLI calc. for a sample trial for two different channels
freqband_PLI_sample_LFPs <- calculate_PLV_PLI_PPC(
  LFP_trial_1_freqband_filtered = sample_LFP_1_freqband_filtered, 
  LFP_trial_2_freqband_filtered = sample_LFP_2_freqband_filtered,
  method = "PLI"
)
freqband_PLV_sample_LFPs

# try PPC calc. for a sample trial for two different channels
freqband_PPC_sample_LFPs <- calculate_PLV_PLI_PPC(
  LFP_trial_1_freqband_filtered = sample_LFP_1_freqband_filtered, 
  LFP_trial_2_freqband_filtered = sample_LFP_2_freqband_filtered,
  method = "PPC"
)
freqband_PPC_sample_LFPs

# compute PLV/PLI/PPC between the channels in each frequency band for all primed and unprimed trials
# calculate_PLV_PLI_PPC_hat <- function(LFP, un_primed_ind, frequency_bands, sampling_rate, filter_order, method) {
#   
#   # filter either primed or unprimed trials
#   LFP_trials <- LFP[, , un_primed_ind]
#   
#   # prepare matrices to store results
#   results <- array(dim = c(dim(LFP)[1], dim(LFP)[1], length(frequency_bands)))
#   sum_results <- array(0 ,dim = c(dim(LFP)[1], dim(LFP)[1], length(frequency_bands)))
#   
#   # loop over all trials
#   for (k in 1:dim(LFP_trials)[3]) {
#     
#     # PLV/PLI/PPC for trial i between each channel
#     for (i in 1:dim(LFP)[1]) {
#       for (j in 1:dim(LFP)[1]) {
#         trial_freqband_filtered_1 <- freqband_filtering(
#           LFP_trial = LFP_trials[i, , k], 
#           frequency_bands = frequency_bands, 
#           sampling_rate = sampling_rate,
#           filter_order = filter_order
#         )
#         trial_freqband_filtered_2 <- freqband_filtering(
#           LFP_trial = LFP_trials[j, , k], 
#           frequency_bands = frequency_bands, 
#           sampling_rate = sampling_rate,
#           filter_order = filter_order
#         )
#         results[i, j, ] <- calculate_PLV_PLI_PPC(
#           LFP_trial_1_freqband_filtered = trial_freqband_filtered_1, 
#           LFP_trial_2_freqband_filtered = trial_freqband_filtered_2,
#           method = method
#         )
#       }
#     }
#     
#     # calculate mean PLV/PLI/PPC over all trials 
#     # see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3674231/#:~:text=PLV%20can%20therefore%20be%20viewed,each%20scaling%20of%20the%20wavelet.
#     sum_results <- sum_results + results
#     avg_results <- sum_results / dim(LFP_trials)[3]
#     
#   }
#   
#   return(avg_results)
#   
# }
# 
# # test for a subset
# temp <- LFP_stationary[, , 1:10]
# temp2 <- unprimed_ind[unprimed_ind <= 10]
# PLV_hat_test <- calculate_PLV_PLI_PPC_hat(
#   LFP = temp,
#   un_primed_ind = temp2,
#   frequency_bands = frequency_bands,
#   sampling_rate = sampling_rate, 
#   filter_order = 2,
#   method = "PLV"
# )
# dim(PLV_hat_test) # dim 1 and 2: channel x channel combi, dim 3: frequency bands
# PLV_hat_test[1:3, 1:3, 1:3]

# compute PLV/PLI/PPC between the channels in each frequency band for all primed and unprimed trials
calculate_PLV_PLI_PPC_hat <- function(LFP, un_primed_ind, frequency_bands, sampling_rate, filter_order, method) {
  
  num_channels <- unique(lengths(LFP))
  num_trials <- length(LFP)
  
  # filter either primed or unprimed trials
  LFP_trials <- LFP[un_primed_ind]
  
  # prepare matrices to store results
  results <- array(dim = c(num_channels, num_channels, length(frequency_bands)))
  sum_results <- array(0 ,dim = c(num_channels, num_channels, length(frequency_bands)))
  
  # loop over all trials
  for (k in 1:num_trials) {
    
    # PLV/PLI/PPC for trial i between each channel
    for (i in 1:num_channels) {
      for (j in 1:num_channels) {
        trial_freqband_filtered_1 <- freqband_filtering(
          LFP_trial = as.numeric(LFP[[k]][[i]]), 
          frequency_bands = frequency_bands, 
          sampling_rate = sampling_rate,
          filter_order = filter_order
        )
        trial_freqband_filtered_2 <- freqband_filtering(
          LFP_trial = as.numeric(LFP[[k]][[j]]), 
          frequency_bands = frequency_bands, 
          sampling_rate = sampling_rate,
          filter_order = filter_order
        )
        results[i, j, ] <- calculate_PLV_PLI_PPC(
          LFP_trial_1_freqband_filtered = trial_freqband_filtered_1, 
          LFP_trial_2_freqband_filtered = trial_freqband_filtered_2,
          method = method
        )
      }
    }
    
    # calculate mean PLV/PLI/PPC over all trials 
    # see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3674231/#:~:text=PLV%20can%20therefore%20be%20viewed,each%20scaling%20of%20the%20wavelet.
    sum_results <- sum_results + results
    avg_results <- sum_results / length(LFP_trials)
    
    print(paste("Trial", k, "of", num_trials, "completed."))
    
  }
  
  return(avg_results)
  
}

# test for a subset
temp <- LFP_stationary[1:10]
temp2 <- unprimed_ind[unprimed_ind <= 10]
PLV_hat_test <- calculate_PLV_PLI_PPC_hat(
  LFP = temp,
  un_primed_ind = temp2,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  method = "PLV"
)
dim(PLV_hat_test) # dim 1 and 2: channel x channel combi, dim 3: frequency bands
PLV_hat_test[1:3, 1:3, 1:3]

plot_heatmap(result_array = PLV_hat_test, frequency_bands = frequency_bands)

# calculate PLV between the channels for all frequency bands for unprimed trials
PLV_hat_unprimed <- calculate_PLV_PLI_PPC_hat(
  LFP = LFP_stationary,
  un_primed_ind = unprimed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  method = "PLV"
)

# calculate PLV between the channels for all frequency bands for primed trials
PLV_hat_primed <- calculate_PLV_PLI_PPC_hat(
  LFP = LFP_stationary,
  un_primed_ind = primed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  method = "PLV"
)

# calculate PLI between the channels for all frequency bands for unprimed trials
PLI_hat_unprimed <- calculate_PLV_PLI_PPC_hat(
  LFP = LFP_stationary,
  un_primed_ind = unprimed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  method = "PLI"
)

# calculate PLI between the channels for all frequency bands for primed trials
PLI_hat_primed <- calculate_PLV_PLI_PPC_hat(
  LFP = LFP_stationary,
  un_primed_ind = primed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  method = "PLI"
)

# calculate PPC between the channels for all frequency bands for unprimed trials
PPC_hat_unprimed <- calculate_PLV_PLI_PPC_hat(
  LFP = LFP_stationary,
  un_primed_ind = unprimed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  method = "PPC"
)

# calculate PPC between the channels for all frequency bands for primed trials
PPC_hat_primed <- calculate_PLV_PLI_PPC_hat(
  LFP = LFP_stationary,
  un_primed_ind = primed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  method = "PPC"
)










