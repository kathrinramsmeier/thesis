
# Calculation of the Phase Coherence, Phase Locking Value (PLV), Phase Lag Index (PLI) and Pairwise Phase Consistency (PPC).
# Coherence and PLV are equivalent in the Gaussian case (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3674231/#:~:text=PLV%20can%20therefore%20be%20viewed,each%20scaling%20of%20the%20wavelet.).
# but here: individual LFP trials are not normally distributed.



# look at the 6th and 14th channel of trial 333 as an example first
trial <- 333
channel <- 6
sample_LFP_1 <- LFP_stationary[channel, , trial]
channel <- 14
sample_LFP_2 <- LFP_stationary[channel, , trial]
plot(sample_LFP_1, type = "l")
lines(sample_LFP_2, type = "l", col = "red")



# Filter the LFP in the different frequency bands -------------------------

freqband_filtering <- function(LFP_trial, frequency_bands, sampling_rate, filter_order) {
  
  LFP_freqband_filtered <- matrix(nrow = length(LFP_trial), ncol = length(frequency_bands))
  colnames(LFP_freqband_filtered) <- names(frequency_bands)
  
  nyquist_freq <- sampling_rate / 2
  
  # avoid numerical problem
  frequency_bands[[1]][1] <- 0.01
  
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



# Coherence ---------------------------------------------------------------

# calculate coherence between two signals
calculate_coherence <- function(LFP_trial_1_freqband_filtered, LFP_trial_2_freqband_filtered, meaned) {
  
  coherences <- matrix(nrow = 9, ncol = ncol(LFP_trial_1_freqband_filtered))
  
  # loop over all frequency bands
  for (i in 1:ncol(LFP_trial_1_freqband_filtered)) {
    
    # cross-power spectral magnitudes and density autospectral densities
    # formulas from:
    # https://en.wikipedia.org/wiki/Coherence_(signal_processing)
    # https://www.psychologie.hhu.de/fileadmin/redaktion/Fakultaeten/Mathematisch-Naturwissenschaftliche_Fakultaet/Psychologie/CompPsy/Papers/vinck2010_wingerden.pdf
    csd <- gsignal::cpsd(cbind(LFP_trial_1_freqband_filtered[, i], LFP_trial_2_freqband_filtered[, i]))$cross
    csd_magnitudes <- abs(csd)
    asd_1 <- gsignal::cpsd(cbind(LFP_trial_1_freqband_filtered[, i], LFP_trial_1_freqband_filtered[, i]))$cross
    asd_2 <- gsignal::cpsd(cbind(LFP_trial_2_freqband_filtered[, i], LFP_trial_2_freqband_filtered[, i]))$cross
    
    # coherence
    coherences[, i] <- csd_magnitudes^2 / (asd_1 * asd_2)
    
  }
  
  # mean coherence in each frequency band !!!! need to check if that is plausible to do !!!!
  if(meaned) {
    coherences <- colMeans(coherences)
    names(coherences) <- names(frequency_bands)
    
  } else {
    colnames(coherences) <- names(frequency_bands)
  }
  
  return(coherences)
  
}

# try coherence calc. for a sample trial for two different channels
freqband_coherences_sample_LFPs <- calculate_coherence(
  LFP_trial_1_freqband_filtered = sample_LFP_1_freqband_filtered, 
  LFP_trial_2_freqband_filtered = sample_LFP_2_freqband_filtered,
  meaned = FALSE
) 
freqband_coherences_sample_LFPs # each coherence value corresponds to the coherence between the two LFP signals at a specific frequency within the frequency band
calculate_coherence(
  LFP_trial_1_freqband_filtered = sample_LFP_1_freqband_filtered, 
  LFP_trial_2_freqband_filtered = sample_LFP_2_freqband_filtered,
  meaned = TRUE
)

# compute coherence between the channels in each frequency band in a sample trial
calculate_trial_coherence <- function(LFP, trial, frequency_bands, sampling_rate, filter_order) {
  
  # prepare matrix to store results
  coherences_trial <- array(dim = c(dim(LFP)[1], dim(LFP)[1], length(frequency_bands)))
  
  # coherence for trial i between each channel
  for (i in 1:dim(LFP)[1]) {
    for (j in 1:dim(LFP)[1]) {
      trial_freqband_filtered_1 <- freqband_filtering(
        LFP_trial = LFP[i, , trial], 
        frequency_bands = frequency_bands, 
        sampling_rate = sampling_rate,
        filter_order = filter_order
      )
      trial_freqband_filtered_2 <- freqband_filtering(
        LFP_trial = LFP[j, , trial], 
        frequency_bands = frequency_bands, 
        sampling_rate = sampling_rate,
        filter_order = filter_order
      )
      coherences_trial[i, j, ] <- calculate_coherence(
        LFP_trial_1_freqband_filtered = trial_freqband_filtered_1, 
        LFP_trial_2_freqband_filtered = trial_freqband_filtered_2,
        meaned = TRUE
      )
    }
  }
  
  return(coherences_trial)
    
}

coherence_sample_trial <- calculate_trial_coherence(
  LFP = LFP_stationary,
  trial = 333, 
  frequency_bands = frequency_bands, 
  sampling_rate = sampling_rate, 
  filter_order = 2
)

dim(coherence_sample_trial) # dim 1 and 2: channel x channel combi, dim 3: frequency bands
coherence_sample_trial[, , 1] # coherence of trial 333 between the different channels in delta band 
coherence_sample_trial[1, 2, 1] # coherence of trial 333 between the channel 1 and 2 in delta band 

# compute coherence between the channels in each frequency band for all primed and unprimed trials
calculate_avg_coherence <- function(LFP, un_primed_ind, frequency_bands, sampling_rate, filter_order) {
  
  # filter either primed or unprimed trials
  LFP_trials <- LFP[, , un_primed_ind]

  # prepare matrices to store results
  coherences <- array(dim = c(dim(LFP)[1], dim(LFP)[1], length(frequency_bands)))
  sum_coherences <- array(0 ,dim = c(dim(LFP)[1], dim(LFP)[1], length(frequency_bands)))
  
  # loop over all trials
  for (k in 1:dim(LFP_trials)[3]) {
    
    # coherence for trial i between each channel
    for (i in 1:dim(LFP)[1]) {
      for (j in 1:dim(LFP)[1]) {
        trial_freqband_filtered_1 <- freqband_filtering(
          LFP_trial = LFP_trials[i, , k], 
          frequency_bands = frequency_bands, 
          sampling_rate = sampling_rate,
          filter_order = filter_order
        )
        trial_freqband_filtered_2 <- freqband_filtering(
          LFP_trial = LFP_trials[j, , k], 
          frequency_bands = frequency_bands, 
          sampling_rate = sampling_rate,
          filter_order = filter_order
        )
        coherences[i, j, ] <- calculate_coherence(
          LFP_trial_1_freqband_filtered = trial_freqband_filtered_1, 
          LFP_trial_2_freqband_filtered = trial_freqband_filtered_2,
          meaned = TRUE
        )
      }
    }
    
    # calculate mean coherence over all trials (need to think again if that is plausible to do!!!)
    sum_coherences <- sum_coherences + coherences
    avg_coherences <- sum_coherences / dim(LFP_trials)[3]
    
  }
  
  return(avg_coherences)
  
}

# test for a subset
temp <- LFP_stationary[, , 1:10]
temp2 <- unprimed_ind[unprimed_ind <= 10]
avg_coherence_test <- calculate_avg_coherence(
  LFP = temp,
  un_primed_ind = temp2,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2
)
dim(avg_coherence_test) # dim 1 and 2: channel x channel combi, dim 3: frequency bands
avg_coherence_test[1:3, 1:3, 1:3]

# calculate coherence_hat between the channels for all frequency bands for unprimed trials
coherence_hat_unprimed <- calculate_avg_coherence(
  LFP = LFP_stationary,
  un_primed_ind = unprimed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2
)

# calculate coherence_hat between the channels for all frequency bands for primed trials
coherence_hat_primed <- calculate_avg_coherence(
  LFP = LFP_stationary,
  un_primed_ind = primed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2
)

setwd("C:/Users/ramsm/Desktop/Master/Thesis/R/data_results/coherence")
saveRDS(coherence_hat_unprimed, paste0(session, "avg_coherence_unprimed.Rds"))
saveRDS(coherence_hat_primed, paste0(session, "avg_coherence_primed.Rds"))

# plot the average coherence between the channels for the different frequency bands as heatmaps
plot_heatmap <- function(result_array, frequency_bands) {
  
  # preparation of the data
  coherence_df <- as.data.frame(as.table(result_array))
  colnames(coherence_df) <- c("channel1", "channel2", "frequency_band", "Coherence")
  coherence_df$channel1 <- as.numeric(as.factor(coherence_df$channel1))
  coherence_df$channel2 <- as.numeric(as.factor(coherence_df$channel2))
  coherence_df$frequency_band <- factor(coherence_df$frequency_band)
  levels(coherence_df$frequency_band) <- names(frequency_bands)
  
  # heatmap plot
  ggplot(coherence_df, aes(x = channel1, y = channel2, fill = Coherence)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") + 
    facet_wrap( ~ frequency_band, scales = "free") +  
    labs(title = "Avg Coherence Between the Different Channels as a Function of Frequency Bands", x = "", y = "", fill = "Coherence") +
    theme_minimal()
  
}

plot_heatmap(result_array = avg_coherence_test, frequency_bands = frequency_bands)
plot_heatmap(result_array = coherence_hat_unprimed, frequency_bands = frequency_bands)
plot_heatmap(result_array = coherence_hat_primed, frequency_bands = frequency_bands)



# Phase Locking Value (PLV) and Phase Lag Index (PLI) ---------------------

# compute PLV and PLI across the frequency bands
# formulas from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3674231/#:~:text=PLV%20can%20therefore%20be%20viewed,each%20scaling%20of%20the%20wavelet.
calculate_PLV_PLI <- function(LFP_trial_1_freqband_filtered, LFP_trial_2_freqband_filtered, PLV_PLI) {
  
  stopifnot(PLV_PLI %in% c("PLV", "PLI"))
  
  result <- rep(NA, dim(LFP_trial_1_freqband_filtered)[2])
  names(result) <- colnames(LFP_trial_1_freqband_filtered)
  
  # loop over all frequency bands
  for (i in 1:dim(LFP_trial_1_freqband_filtered)[2]) {
    
    # Hilbert transform to get the analytical signal
    LFP_1_freqband_ht <- gsignal::hilbert(LFP_trial_1_freqband_filtered[, i])
    LFP_2_freqband_ht <- gsignal::hilbert(LFP_trial_2_freqband_filtered[, i])
    
    # compute phase difference
    phase_diffs <- pracma::angle(LFP_1_freqband_ht * Conj(LFP_2_freqband_ht))  # arg((LFP_1_band_ht * Conj(LFP_2_band_ht)) / (abs(LFP_1_band_ht) * abs(LFP_2_band_ht)))
    
    if (PLV_PLI == "PLV") {
      # compute the PLV for each frequency band
      result[i] <- abs(mean(exp(1i * phase_diffs)))
    } else {
      # compute the PLI for each frequency band
      result[i] <- abs(mean(sign(phase_diffs)))
    }
  }
  
  return(result)
  
}

# try PLV calc. for a sample trial for two different channels
freqband_PLV_sample_LFPs <- calculate_PLV_PLI(
  LFP_trial_1_freqband_filtered = sample_LFP_1_freqband_filtered, 
  LFP_trial_2_freqband_filtered = sample_LFP_2_freqband_filtered,
  PLV_PLI = "PLV"
)
freqband_PLV_sample_LFPs

# try PLI calc. for a sample trial for two different channels
freqband_PLI_sample_LFPs <- calculate_PLV_PLI(
  LFP_trial_1_freqband_filtered = sample_LFP_1_freqband_filtered, 
  LFP_trial_2_freqband_filtered = sample_LFP_2_freqband_filtered,
  PLV_PLI = "PLI"
)
freqband_PLV_sample_LFPs

# compute PPC: PLV/PLI between the channels in each frequency band for all primed and unprimed trials
calculate_PPC <- function(LFP, un_primed_ind, frequency_bands, sampling_rate, filter_order, PLV_PLI) {
  
  # filter either primed or unprimed trials
  LFP_trials <- LFP[, , un_primed_ind]
  
  # prepare matrices to store results
  results <- array(dim = c(dim(LFP)[1], dim(LFP)[1], length(frequency_bands)))
  sum_results <- array(0 ,dim = c(dim(LFP)[1], dim(LFP)[1], length(frequency_bands)))
  
  # loop over all trials
  for (k in 1:dim(LFP_trials)[3]) {
    
    # PLV/PLI for trial i between each channel
    for (i in 1:dim(LFP)[1]) {
      for (j in 1:dim(LFP)[1]) {
        trial_freqband_filtered_1 <- freqband_filtering(
          LFP_trial = LFP_trials[i, , k], 
          frequency_bands = frequency_bands, 
          sampling_rate = sampling_rate,
          filter_order = filter_order
        )
        trial_freqband_filtered_2 <- freqband_filtering(
          LFP_trial = LFP_trials[j, , k], 
          frequency_bands = frequency_bands, 
          sampling_rate = sampling_rate,
          filter_order = filter_order
        )
        results[i, j, ] <- calculate_PLV_PLI(
          LFP_trial_1_freqband_filtered = trial_freqband_filtered_1, 
          LFP_trial_2_freqband_filtered = trial_freqband_filtered_2,
          PLV_PLI = PLV_PLI
        )
      }
    }
    
    # calculate mean PLV/PLI over all trials 
    # see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3674231/#:~:text=PLV%20can%20therefore%20be%20viewed,each%20scaling%20of%20the%20wavelet.
    sum_results <- sum_results + results
    avg_results <- sum_results / dim(LFP_trials)[3]
    
  }
  
  return(avg_results)
  
}

# test for a subset
temp <- LFP_stationary[, , 1:10]
temp2 <- unprimed_ind[unprimed_ind <= 10]
PPC_test <- calculate_PPC(
  LFP = temp,
  un_primed_ind = temp2,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  PLV_PLI = "PLV"
)
dim(PPC_test) # dim 1 and 2: channel x channel combi, dim 3: frequency bands
PPC_test[1:3, 1:3, 1:3]

plot_heatmap(result_array = PPC_test, frequency_bands = frequency_bands)

# calculate PPC - PLV between the channels for all frequency bands for unprimed trials
PPC_PLV_unprimed <- calculate_PPC(
  LFP = LFP_stationary,
  un_primed_ind = unprimed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  PLV_PLI = "PLV"
)

# calculate PPC - PLV between the channels for all frequency bands for primed trials
PPC_PLV_primed <- calculate_PPC(
  LFP = LFP_stationary,
  un_primed_ind = primed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  PLV_PLI = "PLV"
)

# calculate PPC - PLI between the channels for all frequency bands for unprimed trials
PPC_PLI_unprimed <- calculate_PLV_PLI_hat(
  LFP = LFP_stationary,
  un_primed_ind = unprimed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  PLV_PLI = "PLI"
)

# calculate PPC - PLI between the channels for all frequency bands for primed trials
PPC_PLI_primed <- calculate_PLV_PLI_hat(
  LFP = LFP_stationary,
  un_primed_ind = primed_ind,
  frequency_bands = frequency_bands,
  sampling_rate = sampling_rate, 
  filter_order = 2,
  PLV_PLI = "PLI"
)










