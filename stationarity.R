
# sampling rate
sampling_rate <- recordinginfo$data_sampling_rate_hz

# extract the LFP from the first channel
channel1 <- t(LFP[1, , ])
dim(channel1) # trial, time

# clipping the LFP
time_clipped_ind <- which((time < -10) | (time > 10))
channel1_clipped <- channel1[time_clipped_ind, ]



# Check for stationarity --------------------------------------------------

# Augmented Dickey-Fuller test - H0: non-stationay
adf_stationarity_test <- function(LFP_data) {
  p_values <- rep(NA, nrow(LFP_data))
  for(i in 1:nrow(LFP_data)) {
    p_values[i] <- adf.test(LFP_data[i, ])$p.value
  }
  percentage_p_values <- sum(p_values < 0.05) / nrow(LFP_data)
  return(percentage_p_values)
}
adf_result_baseline <- adf_stationarity_test(channel1_clipped)
adf_result_baseline # most (but not all) are stationary

# Kwiatkowski-Phillips-Schmidt-Shin test - H0: stationary
# kpss.test()

# examining the autocorrelation function of some sample trials
# non-CS variables typically will have an autocorrelation that declines slowly with increasing lags and most CS variables will have a sharply declining autocorrelation function 
examine_autocorr <- function(LFP_data) {
  set.seed(42)
  sample_trials <- sample(1:dim(LFP_data)[1], size = 6)
  par(mfrow = c(3, 2))
  for(i in 1:6){
    acf(LFP_data[i, ])
  }
  par(mfrow = c(1, 1))
}
examine_autocorr(channel1_clipped)


# Dealing with non-stationarity - removing linear trends ------------------

channel1_detrended <- pracma::detrend(t(channel1_clipped), tt = "linear")
channel1_detrended <- t(channel1_detrended)

par(mfrow = c(2, 1))
plot(time, channel1_clipped[1, ], type = "l", xlab = "Time", ylab = "Amplitude", main = "Clipped LFP of one Trial")
plot(time, channel1_detrended[1, ], type = "l", xlab = "Time", ylab = "Amplitude", main = "Detrended and Clipped LFP of one Trial")
par(mfrow = c(1, 1))

# ADF test for stationarity
adf_result_detrended <- adf_stationarity_test(channel1_detrended) 
adf_result_detrended # no difference to baseline

# examining the autocorrelation function of some sample trials
examine_autocorr(channel1_detrended) # no difference to baseline



# Dealing with non-stationarity - removal of temporal mean and div --------
# ision of temporal sd-----------------------------------------------------

# subtract the mean of each LFP trial
channel1_demeaned <- channel1_detrended - rowMeans(channel1_detrended)

# divide by the temporal sd
channel1_normalised <- channel1_demeaned / sd(channel1_demeaned)

par(mfrow = c(2, 1))
plot(time, channel1_clipped[1, ], type = "l", xlab = "Time", ylab = "Amplitude", main = "Clipped LFP of one Trial")
plot(time, channel1_normalised[1, ], type = "l", xlab = "Time", ylab = "Amplitude", main = "Normalised, Detrended and Clipped LFP of one Trial")
par(mfrow = c(1, 1))

# ADF test for stationarity
adf_result_detrended <- adf_stationarity_test(channel1_normalised) 
adf_result_detrended # no difference to baseline

# examining the autocorrelation function of some sample trials
examine_autocorr(channel1_normalised) # no difference to baseline



# Dealing with non-stationarity - ensemble mean removal, division  --------
# of ensemble sd ----------------------------------------------------------

# subtract the ensemble mean
channel1_ensemble_mean_removed <- t(t(channel1_normalised) - colMeans(channel1_normalised))

# divide by the ensemble sd
channel1_ensemble_adjusted <- t(t(channel1_ensemble_mean_removed) / sd(t(channel1_ensemble_mean_removed)))

par(mfrow = c(2, 1))
plot(time, channel1_clipped[1, ], type = "l", xlab = "Time", ylab = "Amplitude", main = "Clipped LFP of one Trial")
plot(time, channel1_ensemble_adjusted[1, ], type = "l", xlab = "Time", ylab = "Amplitude", main = "Normalised, Detrended and Clipped LFP of one Trial, Ensemble Mean Removed, Divided by the Ensemble SD")
par(mfrow = c(1, 1))

# ADF test for stationarity
adf_result_detrended <- adf_stationarity_test(channel1_ensemble_adjusted) 
adf_result_detrended # only slightly better

# examining the autocorrelation function of some sample trials
examine_autocorr(channel1_ensemble_adjusted) # almost no difference to baseline



# Dealing with non-stationarity - differencing ----------------------------

# differencing one time
channel1_differenced <- diff(
  x = t(channel1_ensemble_adjusted),
  lag = 1,
  # order of differencing - remove linear trends
  differences = 1
)

channel1_differenced <- t(channel1_differenced)

par(mfrow = c(2, 1))
plot(time, channel1_clipped_demeaned[1, ], type = "l", xlab = "Time", ylab = "Amplitude", main = "Demeaned and Clipped LFP of one Trial")
plot(time[1 : (length(time) - 1)], channel1_differenced[1, ], type = "l", xlab = "Time", ylab = "Amplitude", main = "Differenced Demeaned and Clipped LFP of one Trial")
par(mfrow = c(1, 1))

# ADF test for stationarity
adf_result_differenced <- adf_stationarity_test(channel1_differenced) 
adf_result_differenced # all are now stationary

# examining the autocorrelation function of some sample trials
examine_autocorr(channel1_differenced) # better

# problem with differencing: change of interpretation - causal interactions are now among changes in the LFPs rather then in the LFPs per se (2) 



# Dealing with non-stationarity - notch filtering -------------------------

# goal: removing line noise using a notch filter

# remove frequencies around 60 Hz (= electrical line noise)
fc <- 60

# Nyquist frequency: maximum frequency that can be represented in the digital signal
nyquist_freq <- sampling_rate / 2
fc_norm <- fc / nyquist_freq

notch_filter <- signal::fir1(
  # filter order - higher -> filters with sharper frequency responses, better attenuation, more computation; lower -> smoother frequency responses, may not achieve as much attenuation
  n = 10, 
  # band edges
  w = fc_norm,
  # notch filter
  type = "stop"
)

# notch filtering each LFP trial
channel1_filtered <- array(dim = dim(channel1_ensemble_adjusted))
for(i in 1:nrow(channel1_filtered)) {
  channel1_filtered[i, ] <- signal::filter(filt = notch_filter, x = channel1_clipped_demeaned[i, ])
}

par(mfrow = c(2, 1))
plot(time, channel1_clipped_demeaned[1, ], type = "l", xlab = "Time", ylab = "Amplitude", main = "Demeaned and Clipped LFP of one Trial")
plot(time, channel1_filtered[1, ], type = "l", xlab = "Time", ylab = "Amplitude", main = "Notch-Filtered, Demeaned and Clipped LFP of one Trial")
par(mfrow = c(1, 1))

# ADF test for stationarity
adf_result_notch <- adf_stationarity_test(channel1_filtered) 
adf_result_notch # almost all LFPs are now stationary

# examining the autocorrelation function of some sample trials
examine_autocorr(channel1_filtered) # not much difference compared to baseline

# # find the best filter order
# adf_result_notch <- rep(NA, 500)
# for(j in 10:500){
#   notch_filter <- signal::fir1(
#     # filter order - higher -> filters with sharper frequency responses, better attenuation, more computation; lower -> smoother frequency responses, may not achieve as much attenuation
#     n = j, 
#     # band edges
#     w = fc_norm,
#     # notch filter
#     type = "stop"
#   )
#   # notch filtering each LFP trial
#   channel1_filtered <- array(dim = dim(channel1_ensemble_adjusted))
#   for(i in 1:nrow(channel1_filtered)) {
#     channel1_filtered[i, ] <- signal::filter(filt = notch_filter, x = channel1_ensemble_adjusted[i, ])
#   }
#   # ADF test for stationarity
#   adf_result_notch[j] <- adf_stationarity_test(channel1_filtered) # less ts are now stationary compared to baseline
# }
# max(adf_result_notch)



### other ideas: windowing, spectral factorisation, wavelet transformations, adaptive recursive least-squares modelling
