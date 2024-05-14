
# Analysis of layer dependent induced spectral power averaged across trials.



# Power Calculation for Frequency Bands (One Example Trial) ---------------

# look at the 6th channel of trial 333 as an example first
channel <- 6
trial <- 333
sample_LFP <- LFP_stationary[channel, , trial]

# Hanning windowing
LFP_hanning_windowed <- sample_LFP * gsignal::hanning(n = length(sample_LFP), method = "periodic")

par(mfrow = c(2, 1))
plot(time[trial, ], sample_LFP, type = "l", xlab = "Time", ylab = "", main = "LFP of one Sample Trial")
plot(time[trial, ], LFP_hanning_windowed, type = "l", xlab = "Time", ylab = "", main = "Hanning Windowed LFP of one Sample Trial")
par(mfrow = c(1, 1))

# Fourier transform to get power
LFP_fftransformed <- fft(LFP_hanning_windowed) # frequencies ranging from 0 Hz to the sampling rate

# spacing between adjacent frequency bins in the resulting power spectrum (ratio of the sampling rate to the length of the data or the filter window used in the Fourier transform operation)
freq_resolution <- sampling_rate / length(LFP_hanning_windowed)

# frequency bins (frequencies at which the power spectrum will be computed)
# these frequencies are evenly spaced across the frequency range from 0 Hz to the Nyquist frequency and they represent the bins in which the power spectrum is organised
# Nyquist-Shannon Sampling Theorem: signal must be sampled at least twice as fast as the bandwidth of the signal to accurately reconstruct the waveform
# -> highest frequency that can be accurately represented in the frequency domain is half of the sampling rate 
freq_bins <- seq(0, nyquist_freq, by = freq_resolution)

# calculate the power spectrum (for the first half of the spectrum)
power <- abs(LFP_fftransformed[1:length(freq_bins)])^2

par(mfrow = c(2, 1))
plot(time[trial, ], sample_LFP, type = "l", xlab = "Time", ylab = "", main = "LFP of one Sample Trial")
plot(power, type = "l", main = "Corresponding FFT Power", xlab = "Frequency (Hz)", ylab = "Power") # there are dominant frequencies in the lower indices
par(mfrow = c(1, 1))

# calculate the sum of power within each frequency band
bands_power <- rep(0, length(frequency_bands))
names(bands_power) <- names(frequency_bands)
for (i in 1:length(frequency_bands)){
  ind <- which(freq_bins >= frequency_bands[[i]][1] & freq_bins < frequency_bands[[i]][2])
  bands_power[i] <- sum(power[ind])
}
bands_power # maybe more observations are needed such that freq_resolution is bigger 



# Approach 1 --------------------------------------------------------------

calculate_freqband_trial_power <- function(LFP_trial, frequency_bands, sampling_rate) {
  
  # Hanning windowing
  LFP_hanning_windowed <- LFP_trial * gsignal::hanning(n = 30, method = "periodic")
  
  # Fourier transform to get power
  LFP_fftransformed <- fft(LFP_hanning_windowed) 
  
  # frequency bins
  freq_resolution <- sampling_rate / length(LFP_hanning_windowed)
  freq_bins <- seq(0, nyquist_freq, by = freq_resolution)
  
  # calculate the power spectrum for the first half of the spectrum
  power <- abs(LFP_fftransformed[1:length(freq_bins)])^2
  # power <- (abs(LFP_fftransformed)^2)[1:round(length(LFP_fftransformed)/2)]
  
  # plot(freq_bins, power, type = "l", xlab = "Hz", ylab = "power")
  
  # calculate the sum of power within each frequency band
  bands_power <- rep(0, length(frequency_bands))
  names(bands_power) <- names(frequency_bands)
  for (i in 1:length(frequency_bands)){
    ind <- which(freq_bins >= frequency_bands[[i]][1] & freq_bins < frequency_bands[[i]][2])
    bands_power[i] <- sum(power[ind])
  }
  
  # Return the power within each band
  return(bands_power)
  
}

calculate_freqband_trial_power(
  LFP_trial = sample_LFP,
  frequency_bands = frequency_bands, 
  sampling_rate = sampling_rate
)



# Approach 2 --------------------------------------------------------------

calculate_freqband_trial_power <- function(LFP_trial, frequency_bands, sampling_rate, filter_order) {

  freqband_power <- rep(0, length(frequency_bands))
  names(freqband_power) <- names(frequency_bands)
  
  nyquist_freq <- sampling_rate / 2
  
  # avoid numerical problem
  frequency_bands[[1]][1] <- 0.01
  
  # loop over all frequency bands
  for (i in 1:length(frequency_bands)) {
    
    freq_band <- c(frequency_bands[[i]][1], frequency_bands[[i]][2])
    
    # butterworth filtering
    butterworth <- gsignal::butter(n = filter_order, w = freq_band / nyquist_freq, type = "pass")
    
    # forward and reverse filtering (correction for phase distortion)
    LFP_freqband_filtered <- gsignal::filtfilt(filt = butterworth$b, a = butterworth$a, x = LFP_trial)
    
    # calculate power for each frequency band
    freqband_power[i] <- sum(abs(LFP_freqband_filtered)^2)
  }
  
  return(freqband_power)
  
}

calculate_freqband_trial_power(
  LFP_trial = sample_LFP,
  frequency_bands = frequency_bands, 
  sampling_rate = sampling_rate, 
  filter_order = 2 
) 



# Power of LFP for Bands, Channels and Priming Status (Approach 1) --------

calculate_bands_power <- function(LFP, electrode_channel, un_primed_ind, hanning_windowed, frequency_bands) {
  
  # filter either primed or unprimed trials and the electrode channel
  LFP_trials <- LFP[electrode_channel, , un_primed_ind]
  
  # prepare matrix for results
  bands_power <- matrix(nrow = ncol(LFP_trials), ncol = length(frequency_bands))
  colnames(bands_power) <- names(frequency_bands)
  rownames(bands_power) <- paste0("trial_", 1:ncol(LFP_trials))
  
  # iterate over all trials
  for(i in 1:ncol(LFP_trials)) {
    
    trial <- LFP_trials[, i]
    
    if(hanning_windowed) {
      # Hanning windowing
      trial <- trial * gsignal::hanning(n = 30, method = "periodic")
    }
    
    # Fourier transform to get power
    LFP_fftransformed <- fft(trial)
    
    # frequency bins
    freq_resolution <- sampling_rate / length(LFP_hanning_windowed)
    freq_bins <- seq(0, nyquist_freq, by = freq_resolution)
    
    # calculate the power spectrum for the first half of the spectrum
    power <- abs(LFP_fftransformed[1:length(freq_bins)])^2
    
    # calculate the sum of power within each frequency band
    bands_power_trial_i <- rep(0, length(frequency_bands))
    for (j in 1:length(frequency_bands)){
      ind <- which(freq_bins >= frequency_bands[[j]][1] & freq_bins < frequency_bands[[j]][2])
      bands_power_trial_i[j] <- sum(power[ind])
    }
    
    bands_power[i, ] <- bands_power_trial_i
    
  }
  
  return(bands_power) 
  
}

# testing on channel 6, unprimed trials
rel_bands_power_unprimed_ch6 <- calculate_bands_power(
  LFP = LFP_stationary,
  frequency_bands = frequency_bands,
  un_primed_ind = unprimed_ind,
  hanning_windowed = FALSE,
  electrode_channel = channel
)
head(rel_bands_power_unprimed_ch6)
apply(rel_bands_power_unprimed_ch6, 2, mean)

# mean bands power in every channel
calculate_mean_bands_power <- function(LFP, un_primed_ind, hanning_windowed, frequency_bands) {
  
  # prepare matrix to store the results
  rel_bands_power_un_primed <- matrix(nrow = dim(LFP)[1], ncol = length(frequency_bands))
  colnames(rel_bands_power_un_primed) <- names(frequency_bands)
  rownames(rel_bands_power_un_primed) <- 1:dim(LFP)[1]
  
  # take the mean over all trials and loop over all channels (mean bands power in every channel)
  for(i in 1:dim(LFP)[1]){
    rel_bands_power_un_primed_ch <- calculate_bands_power(
      LFP = LFP,
      frequency_bands = frequency_bands, 
      un_primed_ind = un_primed_ind,
      hanning_windowed = hanning_windowed,
      electrode_channel = i
    )
    rel_bands_power_un_primed[i, ] <- apply(rel_bands_power_un_primed_ch, 2, mean)
  }
  
  return(rel_bands_power_un_primed)
  
}

# mean band power in every channel for unprimed trials
mean_bands_power_unprimed <- calculate_mean_bands_power(
  LFP = LFP_stationary, 
  un_primed_ind = unprimed_ind,
  hanning_windowed = TRUE,
  frequency_bands = frequency_bands
)

# mean band power in every channel for primed trials
mean_bands_power_primed <- calculate_mean_bands_power(
  LFP = LFP_stationary, 
  un_primed_ind = primed_ind,
  hanning_windowed = TRUE,
  frequency_bands = frequency_bands
)

# heatmaps
plot_power_heatmaps <- function(mean_bands_power_unprimed, mean_bands_power_primed, un_primed) {
  
  min_value <- min(c(mean_bands_power_unprimed, mean_bands_power_primed))
  max_value <- max(c(mean_bands_power_unprimed, mean_bands_power_primed))
  breaks <- seq(min_value, max_value, length.out = 100)
  
  if (un_primed == "unprimed") {
    
    row_order <- rownames(mean_bands_power_unprimed)
    col_order <- colnames(mean_bands_power_unprimed)
    
    plot <- pheatmap::pheatmap(
      mean_bands_power_unprimed,
      main = paste0("Relative LFP Power Heatmap in Unprimed Conditions (Session ", session, ")"),
      cluster_rows = FALSE,  
      cluster_cols = FALSE,  
      row_order = row_order,
      col_order = col_order,
      breaks = breaks
    )
    
  } else if (un_primed == "primed") {
    
    row_order <- rownames(mean_bands_power_primed)
    col_order <- colnames(mean_bands_power_primed)
    
    plot <- pheatmap::pheatmap(
      mean_bands_power_primed,
      main = paste0("Relative LFP Power Heatmap in Primed Conditions (Session ", session, ")"),
      cluster_rows = FALSE,  
      cluster_cols = FALSE,  
      row_order = row_order,
      col_order = col_order,
      breaks = breaks
    )
    
  }
  
  print(plot)
  
}

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

# plot power in each channel for primed and unprimed trials for each frequency band
plot_power <- function(mean_bands_power_unprimed, mean_bands_power_primed) {
  
  min_val <- min(min(mean_bands_power_unprimed), min(mean_bands_power_primed))
  max_val <- max(max(mean_bands_power_unprimed), max(mean_bands_power_primed))
  
  library(gridExtra)
  
  plots <- list()
  
  for (i in 1:length(frequency_bands)) {
    plots[[i]] <- ggplot(aes(x = !!as.name(names(frequency_bands)[i]), y = 1:15), data = as.data.frame(mean_bands_power_primed)) +
      geom_point(size = 3, colour = "blue") +
      geom_point(aes(x = !!as.name(names(frequency_bands)[i]), y = 1:15), data = as.data.frame(mean_bands_power_unprimed), size = 3, colour = "lightblue") +
      scale_y_continuous(breaks = seq(1, 15, 1)) +
      xlim(min_val, max_val) +
      labs(x = "Power", y = "") +
      theme_minimal() +
      theme(axis.text.x = element_blank()) 
  }
  
  plot <- grid.arrange(grobs = plots, ncol = 2)
  print(plot)
  
}

plot_power(
  mean_bands_power_unprimed = mean_bands_power_unprimed, 
  mean_bands_power_primed = mean_bands_power_primed
) # primed is dark blue, unprimed light blue



# Power of LFP for Bands, Channels and Priming Status (Approach 2) --------

calculate_freqband_power <- function(LFP, electrode_channel, un_primed_ind, frequency_bands = frequency_bands) {
  
  # filter either primed or unprimed trials and the electrode channel
  LFP_trials <- LFP[electrode_channel, , un_primed_ind]
  
  # prepare matrix for results
  bands_power <- matrix(nrow = ncol(LFP_trials), ncol = length(frequency_bands))
  colnames(bands_power) <- names(frequency_bands)
  rownames(bands_power) <- paste0("trial_", 1:ncol(LFP_trials))
  
  # iterate over all trials
  for(i in 1:ncol(LFP_trials)) {
    # calculate bands power for trial i and fixed channel
    bands_power[i, ] <- calculate_freqband_trial_power(
      LFP_trial = LFP_trials[, i],
      frequency_bands = frequency_bands, 
      sampling_rate = sampling_rate, 
      filter_order = 2
    ) 
  }
  
  return(bands_power) 
  
}

# testing on channel 6, unprimed trials
rel_bands_power_unprimed_ch6 <- calculate_freqband_power(
  LFP = LFP_stationary,
  frequency_bands = frequency_bands,
  un_primed_ind = unprimed_ind,
  electrode_channel = channel
)
head(rel_bands_power_unprimed_ch6)
apply(rel_bands_power_unprimed_ch6, 2, mean)

# mean bands power in every channel
calculate_mean_freqbands_power <- function(LFP, un_primed_ind, frequency_bands) {
  
  # prepare matrix to store the results
  rel_bands_power_un_primed <- matrix(nrow = dim(LFP)[1], ncol = length(frequency_bands))
  colnames(rel_bands_power_un_primed) <- names(frequency_bands)
  rownames(rel_bands_power_un_primed) <- 1:dim(LFP)[1]
  
  # take the mean over all trials and loop over all channels (mean bands power in every channel)
  for(i in 1:dim(LFP)[1]){
    rel_bands_power_un_primed_ch <- calculate_freqband_power(
      LFP = LFP,
      frequency_bands = frequency_bands, 
      un_primed_ind = un_primed_ind,
      electrode_channel = i
    )
    rel_bands_power_un_primed[i, ] <- apply(rel_bands_power_un_primed_ch, 2, mean)
  }
  
  return(rel_bands_power_un_primed)
  
}

# mean band power in every channel for unprimed trials
mean_bands_power_unprimed <- calculate_mean_freqbands_power(
  LFP = LFP_stationary, 
  un_primed_ind = unprimed_ind,
  frequency_bands = frequency_bands
)

# mean band power in every channel for primed trials
mean_bands_power_primed <- calculate_mean_freqbands_power(
  LFP = LFP_stationary, 
  un_primed_ind = primed_ind,
  frequency_bands = frequency_bands
)

# heatmaps
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

# plot power in each channel for primed and unprimed trials for each frequency band
plot_power(
  mean_bands_power_unprimed = mean_bands_power_unprimed, 
  mean_bands_power_primed = mean_bands_power_primed
) # primed is dark blue, unprimed light blue
