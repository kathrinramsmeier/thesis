
# extract the LFP from the first channel
channel1 <- LFP[1, , ]
dim(channel1) # time, trial



# Check for stationarity --------------------------------------------------

# Augmented Dickey-Fuller test - H0: non-stationary
adf_stationarity_test <- function(LFP_data) {
  p_values <- rep(NA, nrow(LFP_data))
  for(i in 1:nrow(LFP_data)) {
    p_values[i] <- adf.test(LFP_data[, i])$p.value
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
    acf(LFP_data[, i])
  }
  par(mfrow = c(1, 1))
}
examine_autocorr(channel1_clipped)


# Dealing with non-stationarity - removing linear trends ------------------

channel1_detrended <- pracma::detrend(t(channel1_clipped), tt = "linear")
channel1_detrended <- t(channel1_detrended)

par(mfrow = c(2, 1))
plot(time, channel1_clipped[, 1], type = "l", xlab = "Time", ylab = "Amplitude", main = "Clipped LFP of one Trial")
plot(time, channel1_detrended[, 1], type = "l", xlab = "Time", ylab = "Amplitude", main = "Detrended and Clipped LFP of one Trial")
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
plot(time, channel1_clipped[, 1], type = "l", xlab = "Time", ylab = "Amplitude", main = "Clipped LFP of one Trial")
plot(time, channel1_normalised[, 1], type = "l", xlab = "Time", ylab = "Amplitude", main = "Normalised, Detrended and Clipped LFP of one Trial")
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
plot(time, channel1_clipped[, 1], type = "l", xlab = "Time", ylab = "Amplitude", main = "Clipped LFP of one Trial")
plot(time, channel1_ensemble_adjusted[, 1], type = "l", xlab = "Time", ylab = "Amplitude", main = "Normalised, Detrended and Clipped LFP of one Trial, Ensemble Mean Removed, Divided by the Ensemble SD")
par(mfrow = c(1, 1))

# ADF test for stationarity
adf_result_detrended <- adf_stationarity_test(channel1_ensemble_adjusted) 
adf_result_detrended # only slightly better

# examining the autocorrelation function of some sample trials
examine_autocorr(channel1_ensemble_adjusted) # almost no difference to baseline



# Dealing with non-stationarity - differencing ----------------------------

# differencing one time
channel1_differenced <- diff(
  x = channel1_ensemble_adjusted,
  lag = 1,
  # order of differencing - remove linear trends
  differences = 1
)

par(mfrow = c(2, 1))
plot(time, channel1_clipped[, 1], type = "l", xlab = "Time", ylab = "Amplitude", main = "Demeaned and Clipped LFP of one Trial")
plot(time[1 : (length(time) - 1)], channel1_differenced[, 1], type = "l", xlab = "Time", ylab = "Amplitude", main = "Differenced Demeaned and Clipped LFP of one Trial")
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
for(i in 1:ncol(channel1_filtered)) {
  channel1_filtered[, i] <- signal::filter(filt = notch_filter, x = channel1_ensemble_adjusted[, i])
}

par(mfrow = c(2, 1))
plot(time, channel1_clipped[, 1], type = "l", xlab = "Time", ylab = "Amplitude", main = "Demeaned and Clipped LFP of one Trial")
plot(time, channel1_filtered[, 1], type = "l", xlab = "Time", ylab = "Amplitude", main = "Notch-Filtered, Demeaned and Clipped LFP of one Trial")
par(mfrow = c(1, 1))

# ADF test for stationarity
adf_result_notch <- adf_stationarity_test(channel1_filtered) 
adf_result_notch # almost all LFPs are now stationary

# examining the autocorrelation function of some sample trials
examine_autocorr(channel1_filtered) # not much difference compared to baseline

# exclude the remaining few non-stationary trials
p_values <- rep(NA, nrow(channel1_filtered))
for(i in 1:nrow(channel1_filtered)) {
  p_values[i] <- adf.test(channel1_filtered[, i])$p.value
}
channel1_filtered <- channel1_filtered[, which(p_values < 0.05)]



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



# Stationarity transformation on the whole LFP ----------------------------

transform_stationary <- function(LFP, 
                                 # hard_remove_all = c(TRUE, FALSE), 
                                 time_clipped_ind = time_clipped_ind) {
  
  # prepare matrix for the transformed LFPs
  LFP_stationary <- array(dim = c(length(time_clipped_ind), dim(LFP)[3]))
  
  # prepare list for the remaining non-stationary LFPs before removing them
  stationary_trials_ind <- list()
  
  # loop over all channels
  for(i in 1:dim(LFP)[1]) {
    
    # extract the LFP from the channel i
    channel <- LFP[i, , ]
    
    # removing linear trends
    channel_detrended <- pracma::detrend(t(channel), tt = "linear")
    channel_detrended <- t(channel_detrended)
    
    # removal of temporal mean and division of temporal sd
    channel_demeaned <- channel_detrended - rowMeans(channel_detrended)
    channel_normalised <- channel_demeaned / sd(channel_demeaned)
    
    # ensemble mean removal and division of ensemble sd 
    channel_ensemble_mean_removed <- t(t(channel_normalised) - colMeans(channel_normalised))
    channel_ensemble_adjusted <- t(t(channel_ensemble_mean_removed) / sd(t(channel_ensemble_mean_removed)))
    
    # notch filtering each LFP trial
    fc <- 60 # remove frequencies around 60 Hz (= electrical line noise)
    nyquist_freq <- sampling_rate / 2 # Nyquist frequency: maximum frequency that can be represented in the digital signal
    fc_norm <- fc / nyquist_freq
    notch_filter <- signal::fir1(
      # filter order - higher -> filters with sharper frequency responses, better attenuation, more computation; lower -> smoother frequency responses, may not achieve as much attenuation
      n = 10, 
      # band edges
      w = fc_norm,
      # notch filter
      type = "stop"
    )
    channel_filtered <- array(dim = dim(channel_ensemble_adjusted))
    for(j in 1:ncol(channel_filtered)) {
      channel_filtered[, j] <- signal::filter(filt = notch_filter, x = channel_ensemble_adjusted[, j])
    }
    
    # display the remaining non-stationary LFPs
    p_values <- rep(NA, nrow(channel_filtered))
    for (j in 1:nrow(channel_filtered)) {
      p_values[j] <- adf.test(channel_filtered[, j])$p.value
    }
    # channel_name <- paste0("channel", i)
    # stationary_trials_ind[[channel_name]] <- channel_filtered[, which(p_values > 0.05)]
    
    print(paste0("Channel ", i, ": Number of LFP trials still non-stationary: ", sum((p_values > 0.05))))
    
    # # exclude the remaining few non-stationary trials
    # if (hard_remove_all) {
    #   channel_filtered <- channel_filtered[, which(p_values <= 0.05)]
    # }
    
    LFP_stationary <- abind::abind(LFP_stationary, channel_filtered, along = 3)
    
  }
  
  LFP_stationary <- LFP_stationary[, , 2:(dim(LFP)[1] + 1)]
  LFP_stationary <- aperm(LFP_stationary, c(3, 1, 2))
  return(LFP_stationary)
  # return(stationary_channels_ind)
  
}

LFP_stationary <- transform_stationary(LFP = LFP, time_clipped_ind = time_clipped_ind)
dim(LFP_stationary)
LFP_stationary[1:3, 1:3, 1:3]













