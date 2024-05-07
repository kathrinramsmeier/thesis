
# Analysis of layer dependent induced spectral power averaged across trials



# Power calculation for frequency bands (one example trial) ---------------

# look at the 6th channel and 333 trial as an example first
channel <- 6
trial <- 333
sample_LFP <- LFP_stationary[channel, , trial]

# # ensure the length of the trial is even by padding with zeros if needed
# if (length(LFP[channel, , trial]) %% 2 != 0) {
#   LFP[channel, , trial] <- c(LFP[channel, , trial], 0)
# }

# Hanning windowing
LFP_hanning_windowed <- sample_LFP * gsignal::hanning(length(sample_LFP))

par(mfrow = c(2, 1))
plot(time, sample_LFP, type = "l", xlab = "Time", ylab = "", main = "LFP of one Sample Trial")
plot(LFP_hanning_windowed, type = "l", xlab = "Time", ylab = "", main = "Hanning Windowed LFP of one Sample Trial")
par(mfrow = c(1, 1))

# Fourier transform
LFP_fftransformed <- fft(LFP_hanning_windowed)

magnitudes <- abs(LFP_fftransformed[1 : (length(LFP_fftransformed) / 2)])
power <- magnitudes^2

par(mfrow = c(2, 1))
plot(time, sample_LFP, type = "l", xlab = "Time", ylab = "", main = "LFP of one Sample Trial")
plot(power, type = "l", main = "Corresponding FFT Power", xlab = "Frequency (Hz)", ylab = "Power") # there are dominant frequencies in the lower indices
par(mfrow = c(1, 1))

# frequency axis
freq_axis <- seq(0, sampling_rate, length.out = length(LFP_fftransformed) / 2)

# calculate the sum of power within each frequency band
bands_power <- rep(0, length(frequency_bands))
names(bands_power) <- names(frequency_bands)
for (i in 1:length(frequency_bands)){
  ind <- which(freq_axis >= frequency_bands[[i]][1] & freq_axis < frequency_bands[[i]][2])
  bands_power[i] <- sum(power[ind])
}

# normalise power in each frequency band
normalised_bands_power <- bands_power / sum(bands_power)
normalised_bands_power






# width <- 1000
# step <- 5
# 
# # windowing (each row corresponds to a windowed section)
# LFP_windowed <- rollapply(
#   sample_LFP,
#   width = width,
#   FUN = identity,
#   by = step,
#   align = "left"
# )
# 
# # Hanning windowing each windowed section
# LFP_hanning_windowed <- apply(LFP_windowed, 1, function(x) x * gsignal::hanning(length(x)))
# 
# # FFT each Hanning windowed section
# LFP_fftransformed <- apply(LFP_hanning_windowed, 2, function(x) fft(x))
# 
# # extract power and magnitudes
# magnitudes <- abs(LFP_fftransformed[1:(nrow(LFP_fftransformed) / 2), ])
# power <- magnitudes^2
# 
# freq_axis <- seq(0, sampling_rate / 2, length.out = nrow(magnitudes))
# 
# mean_power_bands <- matrix(nrow = 500, ncol = length(frequency_bands))
# colnames(mean_power_bands) <- names(frequency_bands)
# for (i in names(frequency_bands)) {
#   band <- frequency_bands[[i]]
#   band_indices <- freq_axis >= band[1] & freq_axis < band[2]
#   band_power <- power[, which(band_indices)]
#   mean_power_bands[, i] <- rowSums(band_power) / sum(band_indices)
# }
# # each value = average power within the corresponding frequency band for one window
# 
# # iterate over each frequency band
# for (i in seq_along(frequency_bands)) {
#   average_power_bands[i] <- mean(mean_power_bands[, i])
# }
# names(average_power_bands) <- names(frequency_bands)
# average_power_bands



# Power of LFP for bands, channels and priming status ---------------------

bands_power_calc <- function(
    LFP,
    electrode_channel,
    un_primed_ind = c("primed_ind", "unprimed_ind"),
    normalised = c(TRUE, FALSE),
    frequency_bands = frequency_bands
) {
  
  # filter either primed or unprimed trials and the electrode channel
  LFP_trials <- LFP[electrode_channel, , un_primed_ind]
  
  # prepare matrix for results
  bands_power <- matrix(nrow = ncol(LFP_trials), ncol = length(frequency_bands))
  colnames(bands_power) <- names(frequency_bands)
  rownames(bands_power) <- paste0("trial_", 1:ncol(LFP_trials))
  
  # iterate over all trials
  for(i in 1:ncol(LFP_trials)) {
    
    trial <- LFP_trials[, i]
    
    # # ensure the length of trial is even by padding with zeros if needed
    # if (length(trial) %% 2 != 0) {
    #   trial <- c(trial, 0)
    # }
    
    # windowing using Hanning window
    LFP_hanning_windowed <- trial * gsignal::hanning(length(trial))
    
    # Fourier transform
    LFP_fftransformed <- fft(LFP_hanning_windowed)
    
    # calculate the power
    magnitudes <- abs(LFP_fftransformed[1 : (length(LFP_fftransformed) / 2)])
    power <- magnitudes^2
    
    # frequency axis
    freq_axis <- seq(0, sampling_rate, length.out = length(LFP_fftransformed) / 2)
    
    # calculate the sum of power within each frequency band
    bands_power_j <- rep(0, length(frequency_bands))
    names(bands_power_j) <- names(frequency_bands)
    for (j in 1:length(frequency_bands)){
      ind <- which(freq_axis >= frequency_bands[[j]][1] & freq_axis < frequency_bands[[j]][2])
      bands_power_j[j] <- sum(power[ind])
    }
    
    if(normalised) {
      # relative power in each frequency band
      bands_power[i, ] <- bands_power_j / sum(bands_power_j)
    } else {
      bands_power[i, ] <- bands_power_j
    }
    
  }
  
  return(bands_power) 
  
}

# # testing on channel 6, unprimed trials
# rel_bands_power_unprimed_ch1 <- bands_power_calc(
#   LFP = LFP_stationary,
#   frequency_bands = frequency_bands,
#   un_primed_ind = unprimed_ind,
#   electrode_channel = channel,
#   normalised = FALSE
# )
# head(rel_bands_power_unprimed_ch1)
# apply(rel_bands_power_unprimed_ch1, 2, mean)

# take the mean over all trials and loop over all channels (unprimed trials)
# (mean bands power in every channel)
rel_bands_power_unprimed <- matrix(nrow = dim(LFP_stationary)[1], ncol = length(frequency_bands))
colnames(rel_bands_power_unprimed) <- names(frequency_bands)
rownames(rel_bands_power_unprimed) <- 1:dim(LFP_stationary)[1]
for(i in 1:dim(LFP_stationary)[1]){
  rel_bands_power_unprimed_ch <- bands_power_calc(
    LFP = LFP_stationary,
    frequency_bands = frequency_bands, 
    un_primed_ind = unprimed_ind,
    electrode_channel = i,
    normalised = FALSE
  )
  rel_bands_power_unprimed[i, ] <- apply(rel_bands_power_unprimed_ch, 2, mean)
}
rel_bands_power_unprimed

# take the mean over all trials and loop over all channels (primed trials)
# (mean bands power in every channel)
rel_bands_power_primed <- matrix(nrow = dim(LFP_stationary)[1], ncol = length(frequency_bands))
colnames(rel_bands_power_primed) <- names(frequency_bands)
rownames(rel_bands_power_primed) <- 1:dim(LFP_stationary)[1]
for(i in 1:dim(LFP_stationary)[1]){
  rel_bands_power_primed_ch <- bands_power_calc(
    LFP = LFP_stationary,
    frequency_bands = frequency_bands, 
    un_primed_ind = primed_ind,
    electrode_channel = i,
    normalised = FALSE
  )
  rel_bands_power_primed[i, ] <- apply(rel_bands_power_primed_ch, 2, mean)
}
rel_bands_power_primed

# heatmaps
row_order <- rownames(rel_bands_power_unprimed)
col_order <- colnames(rel_bands_power_unprimed)
min_value <- min(c(rel_bands_power_unprimed, rel_bands_power_primed))
max_value <- max(c(rel_bands_power_unprimed, rel_bands_power_primed))
breaks <- seq(min_value, max_value, length.out = 100)
pheatmap::pheatmap(
  rel_bands_power_unprimed,
  main = "Relative LFP Power Heatmap in Unprimed Conditions (Session C190127)",
  cluster_rows = FALSE,  
  cluster_cols = FALSE,  
  row_order = row_order,
  col_order = col_order,
  breaks = breaks
)
pheatmap::pheatmap(
  rel_bands_power_primed,
  main = "Relative LFP Power Heatmap in Primed Conditions (Session C190127)",
  cluster_rows = FALSE,  
  cluster_cols = FALSE,  
  row_order = row_order,
  col_order = col_order,
  breaks = breaks
)

# plot power in each channel for primed and unprimed trials for each frequency band
min_val <- min(min(rel_bands_power_unprimed), min(rel_bands_power_primed))
max_val <- max(max(rel_bands_power_unprimed), max(rel_bands_power_primed))
library(gridExtra)
plots <- list()
for (i in 1:length(frequency_bands)) {
  plots[[i]] <- ggplot(aes(x = !!as.name(names(frequency_bands)[i]), y = 1:15), data = as.data.frame(rel_bands_power_primed)) +
    geom_point(size = 3, colour = "blue") +
    geom_point(aes(x = !!as.name(names(frequency_bands)[i]), y = 1:15), data = as.data.frame(rel_bands_power_unprimed), size = 3, colour = "lightblue") +
    scale_y_continuous(breaks = seq(1, 15, 1)) +
    xlim(min_val, max_val) +
    labs(x = "Power", y = "") +
    theme_minimal() +
    theme(axis.text.x = element_blank()) 
}
grid.arrange(grobs = plots, ncol = 2)











