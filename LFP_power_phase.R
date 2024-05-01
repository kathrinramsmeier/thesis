
trial1 <- LFP[1, , 100]

# clipping the LFP - taking away data that is close to the event: exclude 20 milliseconds around the stimulus
time_clipped_ind <- which((time < -10) | (time > 10))
trial1_clipped <- trial1[time_clipped_ind]

# demeaning
trial1_clipped_demeaned <- trial1_clipped - mean(trial1_clipped)

par(mfrow = c(3, 1))
plot(time, trial1, type = "l", xlab = "Time", ylab = "Amplitude", main = "LFP of one Trial")
plot(time[time_clipped_ind], trial1_clipped, type = "l", xlab = "Time", ylab = "Amplitude", main = "Clipped LFP of one Trial")
plot(time[time_clipped_ind], trial1_clipped_demeaned, type = "l", xlab = "Time", ylab = "Amplitude", main = "Demeaned and Clipped LFP of one Trial")
par(mfrow = c(1, 1))

# only investigate the ultimate time after the stimulus
# trial1_clipped_demeaned <- trial1_clipped_demeaned[which(time == 0) : which(round(time) == 500)]

# # kernel smoothing
# dt <- index(trial1)
# y <- coredata(trial1)
# plot(dt, y, type = "l", xlab = "", ylab = "", main = "Kernel Smoothing")
# lines(ksmooth(dt, y, "normal", bandwidth = 20), col = "red", type = "l")

# check for non-stationarity
adf.test(trial1_clipped_demeaned)$p.value # if < 0.05 -> stationary

# frequency bands
# https://ncbi.nlm.nih.gov/pmc/articles/PMC3122299/#:~:text=Neural%20oscillations%20are%20electrical%20activities,gamma%20(%3E80%20Hz).
frequency_bands <- list(
  delta = c(0, 4),
  theta = c(4, 8),
  alpha = c(8, 12),
  beta = c(12, 30),
  low_gamma = c(30, 80),
  high_gamma = c(80, 200)
)

# sampling rate
sampling_rate <- recordinginfo$data_sampling_rate_hz

# filter all primed and unprimed trials
primed_trials <- task$trial_number_count[task$block_trial_count %in% c(1, 2)]
unprimed_trials <- task$trial_number_count[!(task$block_trial_count %in% c(1, 2))]



# Fast Fourier Transform --------------------------------------------------

# ensure the length of trial1 is even by padding with zeros if needed
if (length(trial1_clipped_demeaned) %% 2 != 0) {
  trial1_clipped_demeaned <- c(trial1_clipped_demeaned, 0)
}

# windowing using Hanning window
windowed <- trial1_clipped * gsignal::hanning(length(trial1_clipped))

# Fourier transform
fftransformed <- fft(windowed)

# fftransformed <- fft(trial1_clipped_demeaned)

magnitudes <- abs(fftransformed[1 : (length(fftransformed) / 2)])
power <- magnitudes^2

par(mfrow = c(2, 1))
plot(
  time[time_clipped_ind],
  trial1_clipped_demeaned,
  type = "l",
  xlab = "Time",
  ylab = "Amplitude",
  main = "Demeaned Clipped LFP of one Trial"
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


# Fast Fourier Transformation using segmented LFP -------------------------

# length of each segment
segment_length <- 1000

# overlap percentage between the segments
overlap <- 0.25

# segment the LFP
num_segments <- floor((length(trial1_clipped_demeaned) - segment_length) / (segment_length * (1 - overlap)))
segments <- matrix(0, nrow = num_segments, ncol = segment_length)
for (i in 1:num_segments) {
  start_index <- floor((i - 1) * segment_length * (1 - overlap)) + 1
  end_index <- start_index + segment_length - 1
  segments[i, ] <- trial1_clipped_demeaned[start_index:end_index]
}

# hanning windowing each segment
windowed_segments <- apply(segments, 1, function(seg) seg * gsignal::hanning(length(seg)))

# Fourier transform for each windowed segment
fftransformed <- apply(windowed_segments, 1, function(seg) fft(seg))

magnitudes <- abs(fftransformed[1 : (length(fftransformed) / 2)])
power <- magnitudes^2

par(mfrow = c(2, 1))
plot(
  time[time_clipped_ind],
  trial1_clipped_demeaned,
  type = "l",
  xlab = "Time",
  ylab = "Amplitude",
  main = "Demeaned Clipped LFP of one Trial"
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

normalised_bands_power_calc <- function(
    LFP,
    frequency_bands,
    un_primed_ind,
    electrode_channel,
    time_clipped_ind,
    trials_ind
  ) {
  
  # filter either primed or unprimed trials and the electrode channel
  un_primed_index <- which(trials_ind %in% un_primed_ind)
  trials <- LFP[electrode_channel, , un_primed_index]

  # prepare matrix for results
  normalised_bands_power <- matrix(nrow = dim(trials)[2], ncol = length(frequency_bands))
  colnames(normalised_bands_power) <- names(frequency_bands)
  rownames(normalised_bands_power) <- paste0("trial_", 1:dim(trials)[2])
  
  # iterate over all trials
  for(i in 1:dim(trials)[2]) {
    
    trial <- trials[, i]
    
    # clipping the LFP - taking away data that is close to the event
    trial <- trial[time_clipped_ind]
    
    # demean
    trial <- trial - mean(trial)
    
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
    bands_power <- rep(0, length(frequency_bands))
    names(bands_power) <- names(frequency_bands)
    for (j in 1:length(frequency_bands)){
      ind <- which(freq_axis >= frequency_bands[[j]][1] & freq_axis < frequency_bands[[j]][2])
      bands_power[j] <- sum(power[ind])
    }
    
    # relative power in each frequency band
    normalised_bands_power[i, ] <- bands_power / sum(bands_power)
    
  }
  
  return(normalised_bands_power)
  
}

# test
normalised_bands_power_primed_ch1 <- normalised_bands_power_calc(
  LFP = LFP,
  frequency_bands = frequency_bands,
  un_primed_ind = primed_trials,
  electrode_channel = 1,
  time_clipped_ind = time_clipped_ind,
  trials_ind = trials_ind
)
head(normalised_bands_power_primed_ch1)

# loop over all channels (primed)
normalised_bands_power_primed <- matrix(nrow = dim(LFP)[1], ncol = 6)
colnames(normalised_bands_power_primed) <- names(frequency_bands)
rownames(normalised_bands_power_primed) <- 1:dim(LFP)[1]
for(i in 1:dim(LFP)[1]){
  normalised_bands_power_primed_ch <- normalised_bands_power_calc(
    LFP = LFP,
    frequency_bands = frequency_bands, 
    un_primed_ind = primed_trials,
    electrode_channel = i,
    time_clipped_ind = time_clipped_ind,
    trials_ind = trials_ind
  )
  normalised_bands_power_primed[i, ] <- apply(normalised_bands_power_primed_ch, 2, mean)
}
normalised_bands_power_primed

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

# loop over all channels (unprimed)
normalised_bands_power_unprimed <- matrix(nrow = dim(LFP)[1], ncol = 6)
colnames(normalised_bands_power_unprimed) <- names(frequency_bands)
rownames(normalised_bands_power_unprimed) <- 1:dim(LFP)[1]
for(i in 1:dim(LFP)[1]){
  normalised_bands_power_unprimed_ch <- normalised_bands_power_calc(
    LFP = LFP,
    frequency_bands = frequency_bands, 
    un_primed_ind = unprimed_trials,
    electrode_channel = i,
    time_clipped_ind = time_clipped_ind,
    trials_ind = trials_ind
  )
  normalised_bands_power_unprimed[i, ] <- apply(normalised_bands_power_unprimed_ch, 2, mean)
}
normalised_bands_power_unprimed

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

row_order <- rownames(normalised_bands_power_unprimed)
col_order <- colnames(normalised_bands_power_unprimed)
pheatmap::pheatmap(
  normalised_bands_power_unprimed,
  main = "Relative LFP Power in Unprimed Conditions",
  cluster_rows = FALSE,  
  cluster_cols = FALSE,  
  row_order = row_order,
  col_order = col_order
)
pheatmap::pheatmap(
  normalised_bands_power_primed,
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








# # generate a spectrogram - chopping signal into slices, windowing and Fourier transformation
# spectrogram <- signal::specgram(
#   x = trial1,
#   # size of the fourier transform window
#   # n = 500,
#   # sampling rate - number of samples per unit of time during the LFP recording
#   Fs = sampling_rate,
#   # window size 
#   # window = 100,
#   # overlap
#   # overlap = 100
# )
# 
# # phase information
# p <- abs(spectrogram$S)
# 
# # normalise
# p <- p/max(p)
# 
# # frequency indices
# f <- spectrogram$f
# 
# # calculate total power within each frequency band
# power_spectra <- rep(0, 6)
# names(power_spectra) <- names(frequency_bands)
# for (i in 1:6) {
#   freq_ind <- which(f >= frequency_bands[[i]][1] & f <= frequency_bands[[i]][2])
#   power_spectra[i] <- sum(p[freq_ind])
# }
# power_spectra



# # Fast Fourier Transform --------------------------------------------------
# 
# fftransformed <- fft(trial1)
# 
# par(mfrow = c(2, 1))
# plot(trial1, type = "l", main = "Original Signal")
# plot(
#   Mod(fftransformed),
#   type = "l",
#   main = "FFT Magnitudes",
#   xlab = "Frequency (Hz)",
#   ylab = "Magnitudes"
# ) # there are dominant frequencies in the lower inidces
# par(mfrow = c(1, 1))



# # Band-pass filter --------------------------------------------------------
# 
# # apply band-pass filters for each frequency band
# filtered_signals <- list()
# for(i in 2:6) {
#   # filter frequency
#   W <- c(1 / frequency_bands[i][[1]][1], 1 / frequency_bands[i][[1]][2])  
#   # bandpass filter
#   filtered_data <- dplR::pass.filt(fftransformed, W = W, type = "pass")
#   # Inverse Fourier transform - get back which signal produced this at its Fourier transform 
#   filtered_signals[[i]] <- Re(fft(filtered_data, inverse = TRUE))
# }
# W_delta <- c(0.001, 1 / frequency_bands[1][[1]][2])
# filtered_data_delta <- dplR::pass.filt(fftransformed, W = W_delta, type = "pass")
# filtered_signals[[1]] <- Re(fft(filtered_data_delta, inverse = TRUE))
# 
# ymin <- min(unlist(lapply(filtered_signals, min)))
# ymax <- max(unlist(lapply(filtered_signals, max)))
# par(mfrow = c(3, 2))
# for (i in 1:length(filtered_signals)) {
#   plot(
#     filtered_signals[[i]],
#     type = "l",
#     ylim = c(ymin, ymax),
#     xlab = "Frequency (Hz)",
#     ylab = "",
#     main = names(frequency_bands[i])
#   )
# }
# par(mfrow = c(1, 1))



# # Band-pass filter --------------------------------------------------------
# 
# # Apply band-pass filters for each frequency band
# filtered_signals <- list()
# for (i in 2:6) {
#   # Define the band-pass filter frequency range
#   W <- c((2 * pi) * frequency_bands[[i]][1] / sampling_rate, (2 * pi) * frequency_bands[[i]][2] / sampling_rate)
#   # Band-pass filter
#   filtered_data <- dplR::pass.filt(fftransformed, W = W, type = "pass")
#   # Inverse Fourier transform to get back the filtered signal in the time domain
#   filtered_signals[[i]] <- Re(fft(filtered_data, inverse = TRUE))
# }
# 
# # Define the band-pass filter for the delta band (special case)
# W_delta <- c((2 * pi) * frequency_bands[[1]][1] / sampling_rate, (2 * pi) * frequency_bands[[1]][2] / sampling_rate)
# filtered_data_delta <- dplR::pass.filt(fftransformed, W = W_delta, type = "pass")
# 
# filtered_signals[[1]] <- Re(fft(filtered_data_delta, inverse = TRUE))
# 
# # Find min and max values across all filtered signals for y-axis scaling
# ymin <- min(unlist(lapply(filtered_signals, min)))
# ymax <- max(unlist(lapply(filtered_signals, max)))
# 
# par(mfrow = c(3, 2))
# # Plot filtered signals with consistent y-axis scaling
# for (i in 1:length(filtered_signals)) {
#   plot(
#     filtered_signals[[i]],
#     type = "l",
#     ylim = c(ymin, ymax),
#     xlab = "Frequency (Hz)",
#     ylab = "",
#     main = names(frequency_bands[i])
#   )
# }
# par(mfrow = c(1, 1))




