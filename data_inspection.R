library(reticulate) # load npy data
np <- import("numpy")
library(tidyverse)  # data manipulation
library(vars)       # VAR model
library(tseries)    # ADF test
library(bruceR)     # Granger causality
library(multitaper) # Multitaper spectral analysis



# Session C190127 ---------------------------------------------------------

# LFP
LFP <- np$load("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\lfp_array_uv.npy")

# electrode channel
probe <- read.csv("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\probe.csv", header = FALSE)
probe_header <- read.csv("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\probe_header.csv", header = FALSE)
colnames(probe) <- probe_header
rm(list = "probe_header")

# time
time <- np$load("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\time_array_ms.npy")
time <- time[1, ]

# behaviour
behaviour <- np$load("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\behavior.npy")
behaviour_header <- read.csv("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\behavior_header.csv", header = FALSE)
colnames(behaviour) <- behaviour_header
behaviour <- as.data.frame(behaviour)
rm(list = "behaviour_header")

# task
task <- read.csv("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\task.csv", header = FALSE)
task_header <- read.csv("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\task_header.csv", header = FALSE)
colnames(task) <- task_header
rm(list = "task_header")

# recordings info
recordinginfo <- read.csv("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\recordinginfo.csv", header = FALSE)
recordinginfo_header <- read.csv("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\recordinginfo_header.csv", header = FALSE)
colnames(recordinginfo) <- recordinginfo_header
rm(list = "recordinginfo_header")



# Data inspection ---------------------------------------------------------

LFP[1:5, 1:5, 1:5] # electrode channel (1:15), time, trial
dim(LFP)

summary(time) # 0 is stimulus, before ore-stimulus, after past-stimulus
length(time) # 2 values less than LFP
dim(behaviour)



# Data preprocessing ------------------------------------------------------

# add time value at the start and at the end
time <- c(time[1] - diff(time[c(1, 2)]), time, rev(time)[1] + diff(rev(time)[c(1, 2)]))

# exclude incorrect trials and catch trials (where the colours were the same) from LFP and task
correct_trials_ind <- which(behaviour$accuracy_logical == 1)
non_catch_trails_ind <- task$trial_number_count[(task$catch_trial_logical == 0)]
trials_ind <- intersect(correct_trials_ind, non_catch_trails_ind)
LFP <- LFP[, , trials_ind]
task <- task[trials_ind, ] 

dim(LFP)
dim(task)
length(time)



# Bipolar LFP -------------------------------------------------------------

# test: for first trial
LFP_bipolar[ , ,1] <- apply(LFP[, , 1], 2, function(x) {
  # difference between each channel and its neighbouring channel
  diff(x) 
})

# Initialize an empty matrix with the same dimensions as LFP
LFP_bipolar <- array(NA, dim = c(dim(LFP)[1] - 1, dim(LFP)[2], dim(LFP)[3]))

# Loop over the third dimension of LFP
for (k in 1:dim(LFP)[3]) {
  # Apply the function to compute differences between each channel and its neighboring channel for each trial
  LFP_bipolar[, , k] <- apply(LFP[, , k], 2, function(x) diff(x))
}



# Dealing with non-stationarity -------------------------------------------

# # subtraction of the ensemble means and dividing by the SDs for each trial
# standardize_y_slice <- function(slice) {
#   row_means <- apply(slice, MARGIN = 1:2, mean)
#   row_stds <- apply(slice, MARGIN = 1:2, sd)
#   standardized_slice <- sweep(slice, MARGIN = c(1, 2), STATS = row_means, FUN = "-") / row_stds
#   return(standardized_slice)
# }
# LFP_standardised <- apply(LFP, MARGIN = 3, FUN = standardize_y_slice)



# Identifying the priming condition ---------------------------------------

# i.e., when the macaques was shown red target and green distractor many times consecutively
# in each block the target and distractor colour remained the same, at the end of the block it changed
# so maybe priming is identified in all trials in a block but the last? -> everywhere where block_trial_count is 1 is unprimed

# filter all trials with unprimed and primed condition
unprimed_trials <- task |> 
  dplyr::filter(block_trial_count == 1) |>
  select(trial_number_count) |>
  unlist()
primed_trials <- task |>
  dplyr::filter(block_trial_count != 1) |>
  select(trial_number_count) |>
  unlist()



# Check for non-stationarity ----------------------------------------------

# Augmented-Dickey-Fuller unit root test on each time series
p_values <- rep(0, 15)
for(i in 1:15) {
  p_values[i] <- adf.test(timeseries[, i])$p.value
}
p_values # all > 0.05 -> all seem to be sufficiently stationary






# VAR model incorporating the priming condition ---------------------------

# two separate VAR models? or priming as a variable in the VAR model?

LFP_primed <- LFP_demeaned[, , primed_trials]
LFP_unprimed <- LFP_demeaned[, , unprimed_trials]

# ensemble means for primed and unprimed (mean over trials)
ensemble_means_primed <- apply(LFP_primed, c(1, 2), mean) 
ensemble_means_unprimed <- apply(LFP_unprimed, c(1, 2), mean) 

# create time series for each electrode channel for primed
plot(
  NULL,
  xlim = c(1, ncol(ensemble_means_primed)),
  ylim = c(-200, 200),
  xlab = "Time",
  ylab = "LFP ensemble means (primed)",
  main = "LFP ensemble means for each electrode channel (primed)"
)
timeseries_primed <- matrix(nrow = 15, ncol = dim(ensemble_means_primed)[2])
rownames(timeseries_primed) <- 1:15
for (i in 1:15) {
  ts_obj <- ts(ensemble_means_primed[i, ])
  timeseries_primed[i, ] <- ts_obj
  lines(ts_obj, col = i)
}
legend("topright", legend = 1:15, col = 1:15, lty = 1, title = "Electrode channel")

# create time series for each electrode channel for unprimed
plot(
  NULL,
  xlim = c(1, ncol(ensemble_means_unprimed)),
  ylim = c(-200, 200),
  xlab = "Time",
  ylab = "LFP ensemble means (unprimed)",
  main = "LFP ensemble means for each electrode channel (unprimed)"
)
timeseries_unprimed <- matrix(nrow = 15, ncol = dim(ensemble_means_unprimed)[2])
rownames(timeseries_unprimed) <- 1:15
for (i in 1:15) {
  ts_obj <- ts(ensemble_means_unprimed[i, ])
  timeseries_unprimed[i, ] <- ts_obj
  lines(ts_obj, col = i)
}
legend("topright", legend = 1:15, col = 1:15, lty = 1, title = "Electrode channel")

# select lag for VAR model using AIC for primed
timeseries_primed_t <- t(timeseries_primed)
select_lag <- VARselect(timeseries_primed_t, lag.max = 30)
select_lag$selection
plot(select_lag$criteria[1, ]) # 13 seems to be sufficient already

# select lag for VAR model using AIC for unprimed
timeseries_unprimed_t <- t(timeseries_unprimed)
select_lag <- VARselect(timeseries_unprimed_t, lag.max = 30)
select_lag$selection
plot(select_lag$criteria[1, ]) # 13 seems to be sufficient already

# create the VAR models
VAR_model_primed <- VAR(timeseries_primed, p = 13)
VAR_model_unprimed <- VAR(timeseries_unprimed, p = 13)
summary(VAR_model_primed)
summary(VAR_model_unprimed)



# # ensemble means for both primed and unprimed 
# ensemble_means[primed_trials] <- apply(LFP_primed, c(1, 2), mean)
# ensemble_means[unprimed_trials] <- apply(LFP_unprimed, c(1, 2), mean)
# 
# # create time series for each electrode channel
# timeseries <- matrix(nrow = 15, ncol = dim(ensemble_means)[2])
# rownames(timeseries) <- 1:15
# plot(
#   NULL,
#   xlim = c(1, ncol(ensemble_means)),
#   ylim = range(ensemble_means),
#   xlab = "Time",
#   ylab = "LFP ensemble means",
#   main = "LFPeEnsemble means for each channel"
# )
# for (i in 1:15) {
#   ts_obj <- ts(ensemble_means[i, ])
#   timeseries[i, ] <- ts_obj
#   lines(ts_obj, col = i)
# }
# legend("topright", legend = 1:15, col = 1:15, lty = 1, title = "Channel")
# 
# # select lag for VAR model using AIC
# timeseries <- t(timeseries)
# select_lag <- VARselect(timeseries, lag.max = 30)
# select_lag$selection
# plot(select_lag$criteria[1, ]) # 13 seems to be sufficient already
# 
# # # add priming indicator
# # timeseries$priming_indicator[primed_trials] <- 1
# # timeseries$priming_indicator[unprimed_trials] <- 0
# 
# # create the VAR model
# VAR_model <- VAR(timeseries, p = 13)
# summary(VAR_model)



# Spectrogram, power, phase changes ---------------------------------------

# # Check for non-stationarity: Augmented-Dickey-Fuller unit root test on each time series
# p_values <- matrix(nrow = 15, ncol = dim(LFP_primed)[3])
# for(i in 1:15) {
#   for(j in 1:dim(LFP_primed)[3]) {
#     p_values[i, j] <- adf.test(LFP_primed[i, , j])$p.value
#   }
# }
# p_values # if > 0.05 -> stationary

plot(LFP[1, , 1], type = "l")
# demeaning
trial1 <- LFP[1, , 1] - mean(LFP[1, , 1])
plot(trial1, type = "l", xlab = "time", ylab = "Amplitude")

# check for non-stationarity
adf.test(trial1)$p.value # > 0.05 -> stationary

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


### 1. attempt

periodogramm <- spec.pgram(
  trial1, 
  # sampling rate - number of samples per unit of time during the LFP recording
  Fs = 1/(time[2] - time[1]),
  taper = 0, 
  pad = 0
)

# calculate total power within each frequency band for primed and unprimed condition
power_spectra <- matrix(nrow = 6, ncol = 2)
colnames(power_spectra) <- c("primed", "unprimed")
rownames(power_spectra) <- names(frequency_bands)
for (i in 1:6){
  power_spectra[i, 1] <- sum(periodogramm$spec[
    periodogramm$freq >= frequency_bands[i][[1]][1] &
      periodogramm$freq <= frequency_bands[i][[1]][2]
  ])
}
power_spectra



### 2. attempt

# periodogram of fourier transformed ts
fourier <- fft(trial1)
periodogram <- abs(fourier)^2 / length(timeseries_primed)
plot(periodogram, type = "l", xlab = "Frequency", ylab = "Power", main = "Periodogram")

# extract frequencies from the periodogram
freq <- seq(0, (length(periodogram) - 1) / length(periodogram), 1 / length(periodogram))

power_spectra <- matrix(nrow = 6, ncol = 2)
colnames(power_spectra) <- c("primed", "unprimed")
rownames(power_spectra) <- names(frequency_bands)
for (i in 1:6){
  ind <- which(freq >= frequency_bands[i][[1]][1] & freq < frequency_bands[i][[1]][2])
  power_spectra[i, 1] <- sum(periodogram[ind])
}
power_spectra



### 3. attempt

# generate a spectrogram - chopping signal into slices, windowing and Fourier transformation
spectrogram <- signal::specgram(
  x = trial1,
  # size of the fourier transform window
  # n = 500,
  # sampling rate - number of samples per unit of time during the LFP recording
  Fs = 1000,
  # window size 
  # window = 100,
  # overlap
  # overlap = 100
)

# phase information
P <- abs(spectrogram$S)

# normalise
P = P/max(P)

# convert to dB
P = 10*log10(P)

# config time axis
t = spectrogram$t

# plot spectrogram
oce::imagep(x = t,
       y = spectrogram$f,
       z = t(P),
       col = oce::oce.colorsViridis,
       ylab = 'Frequency',
       xlab = 'Time',
       drawPalette = T,
       decimate = F
)

# extract spectrogram data
spec_freq_ind <- spectrogram$f  
spec_time_ind <- spectrogram$t 
spec_power <- spectrogram$S     

power_spectra <- matrix(nrow = 6, ncol = 2)
colnames(power_spectra) <- c("primed", "unprimed")
rownames(power_spectra) <- names(frequency_bands)
for (i in 1:6){
  ind <- which(spec_freq_ind >= frequency_bands[i][[1]][1] & spec_freq_ind < frequency_bands[i][[1]][2])
  power_spectra[i, 1] <- sum(power[ind])
}
power_spectra




power_spectra <- matrix(nrow = length(frequency_bands), ncol = length(time))
rownames(power_spectra) <- names(frequency_bands)
colnames(power_spectra) <- time

# Calculate power spectra for each frequency band for each time point
for (i in 1:length(frequency_bands)) {
  ind <- freq >= frequency_bands[[i]][1] & freq < frequency_bands[[i]][2]
  power_spectra[i, ] <- apply(power[ind, , drop = FALSE], 2, sum)
}

power_spectra

# Plot power spectra for each frequency band
matplot(
  time,
  power_spectra,
  type = "l",
  xlab = "Time",
  ylab = "Power",
  main = "Power Spectra for Frequency Bands"
)
legend(
  "topright",
  legend = rownames(power_spectra),
  col = 1:length(power_spectra),
  lty = 1
)





# Compute Fourier Transform to get power spectrum
power_spectrum <- abs(fft(trial1))

# Compute Hilbert Transform to get instantaneous phase
phase <- Hilbert(lfp_data$signal)$phase

# Define frequency bands (adjust these bands according to your needs)
bands <- list(delta = c(1, 4),
              theta = c(4, 8),
              alpha = c(8, 12),
              beta = c(12, 30),
              gamma = c(30, 100))

# Define sampling rate and length of your LFP data
sampling_rate <- time[2] - time[1]  # Example sampling rate in Hz
data_length <- length(lfp_data$signal)

# Create frequency vector
frequency_vector <- (0:(data_length - 1)) * sampling_rate / data_length

# Compute power and phase within each frequency band
band_power <- sapply(bands, function(band) {
  freq_indices <- which(frequency_vector >= band[1] & frequency_vector <= band[2])
  sum(power_spectrum[freq_indices])
})

band_phase <- sapply(bands, function(band) {
  phase_indices <- which(frequency_vector >= band[1] & frequency_vector <= band[2])
  mean(phase[phase_indices])
})


### concatenating the trials

y_dim <- dim(LFP)[2] * dim(LFP)[3]
LFP_conc <- matrix(nrow = 15, ncol = y_dim)
for(i in 1:15) {
  LFP_conc[i, ] <- c(t(LFP[i, , ]))
}

trial1 <- LFP_conc[1, ]
plot(trial1, type = "l")

# Check for non-stationarity: Augmented-Dickey-Fuller unit root test on each time series
p_values <- rep(0, 15)
for(i in 1:15) {
  p_values[i] <- adf.test(LFP_conc[i, ])$p.value
}
p_values # if > 0.05 -> stationary

# generate a spectrogram - chopping signal into slices, windowing and Fourier transformation
spectrogram <- signal::specgram(
  x = trial1,
  # size of the fourier transform window
  # n = ,
  # sampling rate - number of samples per unit of time during the LFP recording
  Fs = time[2] - time[1]
)

# generate a spectrogram - chopping signal into slices, windowing and Fourier transformation
spectrogram <- signal::specgram(
  x = trial1,
  # size of the fourier transform window
  # n = ,
  # sampling rate - number of samples per unit of time during the LFP recording
  Fs = time[2] - time[1]
)

# extract spectrogram data
spec_freq <- spectrogram$f  
spec_time <- spectrogram$t 
spec_power <- spectrogram$S     

power_spectra <- matrix(nrow = 6, ncol = 2)
colnames(power_spectra) <- c("primed", "unprimed")
rownames(power_spectra) <- names(frequency_bands)
for (i in 1:6){
  ind <- which(spec_freq >= frequency_bands[i][[1]][1] & spec_freq < frequency_bands[i][[1]][2])
  power_spectra[i, 1] <- sum(spec_power[ind])
}
power_spectra















# NEXT STEPS:

# investigate whether the FT looks that weird
# maybe use multitaper (package: multitaper)

# investigate phase changes

# Granger Causality








