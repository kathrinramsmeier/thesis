
# Power calculation for frequency bands (one example trial) ---------------

trial1_filtered <- LFP_stationary[1, , 100]

# ensure the length of trial1 is even by padding with zeros if needed
if (length(trial1_filtered) %% 2 != 0) {
  trial1_filtered <- c(trial1_filtered, 0)
}

# windowing using Hanning window
windowed <- trial1_filtered * gsignal::hanning(length(trial1_filtered))

# Fourier transform
fftransformed <- fft(windowed)

# fftransformed <- fft(trial1_filtered)

magnitudes <- abs(fftransformed[1 : (length(fftransformed) / 2)])
power <- magnitudes^2

par(mfrow = c(2, 1))
plot(
  time,
  trial1_filtered,
  type = "l",
  xlab = "Time",
  ylab = "Amplitude",
  main = "Preprocessed LFP of one Trial"
)
plot(
  power,
  type = "l",
  main = "FFT Power",
  xlab = "Frequency (Hz)",
  ylab = "Power"
) # there are dominant frequencies in the lower indices
par(mfrow = c(1, 1))

# frequency axis
freq_axis <- seq(0, sampling_rate / 2, length.out = length(fftransformed) / 2)

# calculate the sum of magnitudes within each frequency band
bands_magnitudes <- rep(0, 6)
names(bands_magnitudes) <- names(frequency_bands)
for (i in 1:length(frequency_bands)){
  ind <- which(freq_axis >= frequency_bands[[i]][1] & freq_axis < frequency_bands[[i]][2])
  bands_magnitudes[i] <- sum(magnitudes[ind])
}

# normalise power in each frequency band
bands_power <- bands_magnitudes^2
normalised_bands_power <- bands_power / sum(bands_power)
normalised_bands_power



# Power of LFP for bands, channels and priming status ---------------------

bands_power_calc <- function(
    LFP,
    electrode_channel,
    un_primed_ind = c("primed_ind", "unprimed_ind"),
    normalised = c(TRUE, FALSE),
    frequency_bands = frequency_bands
  ) {
  
  # filter either primed or unprimed trials and the electrode channel
  trials <- LFP[electrode_channel, , un_primed_ind]

  # prepare matrix for results
  if(normalised) {
    normalised_bands_power <- matrix(nrow = dim(trials)[2], ncol = length(frequency_bands))
    colnames(normalised_bands_power) <- names(frequency_bands)
    rownames(normalised_bands_power) <- paste0("trial_", 1:dim(trials)[2])
  } else {
    bands_power <- matrix(nrow = dim(LFP)[3], ncol = length(frequency_bands))
    colnames(bands_power) <- names(frequency_bands)
    rownames(bands_power) <- paste0("trial_", 1:dim(LFP)[3])
  }
  
  # iterate over all trials
  for(i in 1:dim(trials)[2]) {
    
    trial <- trials[, i]
    
    # ensure the length of trial is even by padding with zeros if needed
    if (length(trial) %% 2 != 0) {
      trial <- c(trial, 0)
    }
    
    # windowing using Hanning window
    windowed <- trial * gsignal::hanning(length(trial))
    
    # Fourier transform
    fftransformed <- fft(windowed)
    
    # calculate the power
    magnitudes <- abs(fftransformed[1 : (length(fftransformed) / 2)])
    power <- magnitudes^2
    
    # frequency axis
    freq_axis <- seq(0, sampling_rate / 2, length.out = length(fftransformed) / 2)
    
    # calculate the sum of power within each frequency band
    bands_power_j <- rep(0, length(frequency_bands))
    names(bands_power_j) <- names(frequency_bands)
    for (j in 1:length(frequency_bands)){
      ind <- which(freq_axis >= frequency_bands[[j]][1] & freq_axis < frequency_bands[[j]][2])
      bands_power_j[j] <- sum(power[ind])
    }
    
    if(normalised) {
      # relative power in each frequency band
      normalised_bands_power[i, ] <- bands_power_j / sum(bands_power_j)
    } else {
      bands_power[i, ] <- bands_power_j
    }
    
  }
  
  if(normalised) {
    return(normalised_bands_power)
  } else {
    return(bands_power) 
  }
  
}

# testing on channel 1, unprimed trials
rel_bands_power_unprimed_ch1 <- bands_power_calc(
  LFP = LFP_stationary,
  frequency_bands = frequency_bands,
  un_primed_ind = unprimed_ind,
  electrode_channel = 1,
  normalised = TRUE
)
head(rel_bands_power_unprimed_ch1)

# heatmap 
row_order <- rownames(rel_bands_power_unprimed_ch1[1:50, ])
col_order <- colnames(rel_bands_power_unprimed_ch1)
pheatmap::pheatmap(
  rel_bands_power_unprimed_ch1[1:50, ],
  main = "Relative LFP Power in Unprimed Conditions for the First 50 Trials",
  cluster_rows = FALSE,  
  cluster_cols = FALSE,  
  row_order = row_order,
  col_order = col_order
)

# loop over all channels and take the mean over all trials (unprimed trials)
rel_bands_power_unprimed <- matrix(nrow = dim(LFP_stationary)[1], ncol = length(frequency_bands))
colnames(rel_bands_power_unprimed) <- names(frequency_bands)
rownames(rel_bands_power_unprimed) <- 1:dim(LFP_stationary)[1]
for(i in 1:dim(LFP_stationary)[1]){
  rel_bands_power_unprimed_ch <- bands_power_calc(
    LFP = LFP_stationary,
    frequency_bands = frequency_bands, 
    un_primed_ind = unprimed_ind,
    electrode_channel = i,
    normalised = TRUE
  )
  rel_bands_power_unprimed[i, ] <- apply(rel_bands_power_unprimed_ch, 2, mean)
}
rel_bands_power_unprimed

# # loop over all channels (unprimed) - bipolar LFP
# normalised_bands_power_unprimed_bip <- matrix(nrow = 15, ncol = 6)
# colnames(normalised_bands_power_unprimed_bip) <- names(frequency_bands)
# for(i in 1:15){
#   normalised_bands_power_unprimed_ch <- normalised_bands_power_calc(
#     LFP = LFP,
#     frequency_bands = frequency_bands, 
#     un_primed_ind = unprimed_trials,
#     electrode_channel = i
#   )
#   normalised_bands_power_unprimed_bip[i, ] <- apply(normalised_bands_power_unprimed_ch, 2, mean)
# }

# loop over all channels and take the mean over all trials (primed trials)
rel_bands_power_primed <- matrix(nrow = dim(LFP_stationary)[1], ncol = length(frequency_bands))
colnames(rel_bands_power_primed) <- names(frequency_bands)
rownames(rel_bands_power_primed) <- 1:dim(LFP_stationary)[1]
for(i in 1:dim(LFP_stationary)[1]){
  rel_bands_power_primed_ch <- bands_power_calc(
    LFP = LFP_stationary,
    frequency_bands = frequency_bands, 
    un_primed_ind = primed_ind,
    electrode_channel = i,
    normalised = TRUE
  )
  rel_bands_power_primed[i, ] <- apply(rel_bands_power_primed_ch, 2, mean)
}
rel_bands_power_primed

# # loop over all channels (primed) - bipolar LFP
# normalised_bands_power_primed_bip <- matrix(nrow = 14, ncol = 6)
# colnames(normalised_bands_power_primed_bip) <- names(frequency_bands)
# for(i in 1:14){
#   normalised_bands_power_primed_ch <- normalised_bands_power_calc(
#     LFP = LFP,
#     frequency_bands = frequency_bands, 
#     un_primed_ind = primed_trials,
#     electrode_channel = i
#   )
#   normalised_bands_power_primed[i, ] <- apply(normalised_bands_power_primed_ch, 2, mean)
# }
# normalised_bands_power_primed_bip

# heatmaps
row_order <- rownames(rel_bands_power_unprimed)
col_order <- colnames(rel_bands_power_unprimed)
pheatmap::pheatmap(
  rel_bands_power_unprimed,
  main = "Relative LFP Power in Unprimed Conditions",
  cluster_rows = FALSE,  
  cluster_cols = FALSE,  
  row_order = row_order,
  col_order = col_order
)
pheatmap::pheatmap(
  rel_bands_power_primed,
  main = "Relative LFP Power Heatmap in Primed Conditions",
  cluster_rows = FALSE,  
  cluster_cols = FALSE,  
  row_order = row_order,
  col_order = col_order
)



# Multitaper spectral analysis --------------------------------------------

# demodulates <- demod.dpss(trial1_clipped_demeaned, centreFreq = 1, nw = 3, blockLen = 10)
# 
# # amplitudes at each time point
# amplitudes <- demodulates$amplitude
# plot(amplitudes, type = "l")






# Phase change of LFP for bands, channels and priming status --------------













