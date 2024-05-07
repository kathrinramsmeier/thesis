
# Phase coherence for frequency bands (one example trial) -------------------

# look at the 6th channel and trial 333 and 999 as an example first
trial <- 333
sample_LFP_1 <- LFP_stationary[channel, , trial]
trial <- 999
sample_LFP_2 <- LFP_stationary[channel, , trial]

# extract phase angles in each frequency band
phase_angles_frequency_bands <- function(LFP_trial) {
  
  # Hilbert transform to get the analytical signal
  LFP_hilbert_transformed <- gsignal::hilbert(LFP_trial)
  
  # extract phase angles
  phase_angles <- pracma::angle(LFP_hilbert_transformed)
  
  # extract phase angles in each frequency band
  phase_angles_fb <- list()
  for (i in 1:length(frequency_bands)) {
    phase_angles_fb[[i]] <- phase_angles[frequency_bands[[i]][1] <= freq_axis & freq_axis < frequency_bands[[i]][2]]
  }
  names(phase_angles_fb) <- names(frequency_bands)
  
  return(phase_angles_fb)
  
}

phase_angles_fb_1 <- phase_angles_frequency_bands(sample_LFP_1)
phase_angles_fb_2 <- phase_angles_frequency_bands(sample_LFP_2)

# compute phase coherence using mean resultant length (MRL) in each frequency band
# (phase coherence = phenomenon where a constant phase difference exists between any two signals or waves of the same frequency)
phase_coherence <- function(phase_angles_fb_1, phase_angles_fb_2) {
  mrl <- rep(NA, length(frequency_bands))
  for (i in 1:length(frequency_bands)) {
    # Compute mean resultant length (MRL) between trials
    complex_numbers_1 <- exp(1i * phase_angles_fb_1[[i]])
    complex_numbers_2 <- exp(1i * phase_angles_fb_2[[i]])
    complex_conjugate <- Re(complex_numbers_2) - Im(complex_numbers_2)*1i
    mrl[i] <- abs(mean(complex_numbers_1 * complex_conjugate))
    # https://en.wikipedia.org/wiki/Circular_mean
  }
  return(mrl)
}

phase_coherence(phase_angles_fb_1, phase_angles_fb_2) # higher values: greater coherence



# Phase coherence of different cortical channels --------------------------

# look at trial 333 as an example first
trial <- 333
sample_LFP <- LFP_stationary[, , trial]

# extract phase angles in each frequency band for each channel
phase_angles_channels <- list()
for (i in 1:nrow(sample_LFP)) {
  phase_angles_channels[[i]] <- phase_angles_frequency_bands(sample_LFP[i, ])
}

phase_coherence(phase_angles_channels[[1]], phase_angles_channels[[10]])

# compute pairwise phase coherence using mean resultant length (MRL) in each frequency band
# ...




# next: 
# * analyse how the phase of the signal changes across 
#   * different priming conditions
#   * different channels
# * investigate phase synchrony within specific frequency bands
# * phase coherence plot
# * computation of the PTA (phase triggered average)