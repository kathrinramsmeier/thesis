
# look at the 6th channel and 333 trial as an example first
channel <- 6
trial <- 333



# Check for stationarity --------------------------------------------------

# Augmented Dickey-Fuller test (H0: non-stationary) either for all channels or a specific one
adf_stationarity_test <- function(LFP, channel = "all") {
  
  stopifnot(channel %in% c(1:dim(LFP)[1], "all"))
  stopifnot(length(dim(LFP)) %in% c(2, 3))
  
  if (length(dim(LFP)) == 3) {
    
    # ADF test for all channels
    if(channel == "all") {
      
      p_values <- matrix(NA, nrow = dim(LFP)[3], ncol = dim(LFP)[1])
      colnames(p_values) <- paste0("Channel_", 1:dim(LFP)[1])
      
      # loop over all channels
      for (j in 1:dim(LFP)[1]) {
        
        # loop over all trials
        for (i in 1:dim(LFP)[3]) {
          p_values[i, j] <- adf.test(LFP[j, , i])$p.value
        }
        
      }
      
      # percentage of LFPs that are stationary
      percentage_stationary <- colSums(p_values < 0.05) / nrow(p_values)
      
      # ADF test for one specific channel
    } else if (channel %in% 1:15) {
      
      p_values <- rep(NA, dim(LFP)[3])
      
      for (i in 1:dim(LFP)[3]) {
        p_values[i] <- adf.test(LFP[channel, , i])$p.value
      }
      
      # percentage of LFPs that are stationary
      percentage_stationary <- sum(p_values < 0.05) / dim(LFP)[3]
      
    }
    
  } else if (length(dim(LFP)) == 2) {
    
    p_values <- rep(NA, dim(LFP)[2])
    
    for (i in 1:dim(LFP)[2]) {
      p_values[i] <- adf.test(LFP[, i])$p.value
    }
    
    # percentage of LFPs that are stationary
    percentage_stationary <- sum(p_values < 0.05) / dim(LFP)[2]
    
  }
  
  return(percentage_stationary)
  
}
  
adf_result_baseline <- adf_stationarity_test(LFP = LFP, channel = "all")
adf_result_baseline # most are stationary - but can be better in some channels

# Kwiatkowski-Phillips-Schmidt-Shin test - H0: stationary
# kpss.test()

# examining the autocorrelation function of some sample trials in a given channel
# non-CS variables typically will have an autocorrelation that declines slowly with increasing lags and most CS variables will have a sharply declining autocorrelation function 
examine_autocorr <- function(LFP, channel) {
  set.seed(42)
  if(length(dim(LFP)) == 3) {
    sample_trials <- sample(1:dim(LFP)[3], size = 6)
    par(mfrow = c(3, 2))
    for(i in 1:6){
      acf(LFP[channel, , i])
    }
    par(mfrow = c(1, 1))
  } else if (length(dim(LFP)) == 2) {
    sample_trials <- sample(1:dim(LFP)[2], size = 6)
    par(mfrow = c(3, 2))
    for(i in 1:6){
      acf(LFP[channel, i])
    }
  }
}
examine_autocorr(LFP = LFP, channel = 1)
examine_autocorr(LFP = LFP, channel = 6) # worse than in channel 1



# Dealing with non-stationarity - removing linear trends ------------------

LFP_detrended <- pracma::detrend(LFP[channel, , ], tt = "linear")

par(mfrow = c(2, 1))
plot(time, LFP[channel, , trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Clipped LFP of one Sample Trial")
plot(time, LFP_detrended[, trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Detrended and Clipped LFP of one Sample Trial")
par(mfrow = c(1, 1))
# is this really good or does this maybe delete valuable information?

# ADF test for stationarity
adf_result_detrended <- adf_stationarity_test(LFP = LFP_detrended) 
adf_result_detrended - adf_result_baseline[channel] # no difference to baseline

# examining the autocorrelation function of some sample trials
examine_autocorr(LFP_detrended) # no difference to baseline



# Dealing with non-stationarity - removal of temporal mean and div --------
# ision of temporal sd-----------------------------------------------------

# subtract the mean of each LFP trial
LFP_demeaned <- LFP_detrended - colMeans(LFP_detrended)

# divide by the temporal sd
LFP_normalised <- t(LFP_demeaned) / sd(t(LFP_demeaned))
LFP_normalised <- t(LFP_normalised)

par(mfrow = c(2, 1))
plot(time, LFP[channel, , trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Clipped LFP of one Sample Trial")
plot(time, LFP_normalised[, trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Normalised, Detrended and Clipped LFP of one Sample Trial")
par(mfrow = c(1, 1))

# ADF test for stationarity
adf_result_normalised <- adf_stationarity_test(LFP_normalised) 
adf_result_normalised - adf_result_baseline[channel] # no difference to baseline

# examining the autocorrelation function of some sample trials
examine_autocorr(LFP_normalised) # no difference to baseline



# Dealing with non-stationarity - ensemble mean removal, division  --------
# of ensemble sd ----------------------------------------------------------

# subtract the ensemble mean
LFP_ensemble_mean_removed <- LFP_normalised - rowMeans(LFP_normalised)

# divide by the ensemble sd
LFP_ensemble_adjusted <- LFP_ensemble_mean_removed / sd(LFP_ensemble_mean_removed)

par(mfrow = c(2, 1))
plot(time, LFP[channel, , trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Clipped LFP of one Sample Trial")
plot(time, LFP_ensemble_adjusted[, trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Normalised, Detrended and Clipped LFP of one Sample Trial, Ensemble Mean Removed, Divided by the Ensemble SD")
par(mfrow = c(1, 1))

# ADF test for stationarity
adf_result_detrended <- adf_stationarity_test(LFP_ensemble_adjusted) 
adf_result_detrended - adf_result_baseline[channel] # only slightly better

# examining the autocorrelation function of some sample trials
examine_autocorr(LFP_ensemble_adjusted) # almost no difference to baseline



# Dealing with non-stationarity - differencing ----------------------------

# differencing one time
LFP_differenced <- diff(
  x = LFP_ensemble_adjusted,
  lag = 1,
  # order of differencing - remove linear trends
  differences = 1
)

par(mfrow = c(2, 1))
plot(time, LFP[channel, , trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Demeaned and Clipped LFP of one Sample Trial")
plot(time[1 : (length(time) - 1)], LFP_differenced[, trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Differenced Demeaned and Clipped LFP of one Sample Trial")
par(mfrow = c(1, 1))

# ADF test for stationarity
adf_result_differenced <- adf_stationarity_test(LFP_differenced) 
adf_result_differenced # all are now stationary (by definition)

# examining the autocorrelation function of some sample trials
examine_autocorr(LFP_differenced) # better

# problem with differencing: change of interpretation - causal interactions are now among changes in the LFPs rather then in the LFPs per se (2) 



# Dealing with non-stationarity - notch filtering -------------------------

# goal: removing line noise using a notch filter

# remove frequencies around 60 Hz (= electrical line noise)
line_noise <- 60

# Nyquist frequency: maximum frequency that can be represented in the digital signal
nyquist_freq <- sampling_rate / 2
line_noise_norm <- line_noise / nyquist_freq

notch_filter <- signal::fir1(
  # filter order - higher -> filters with sharper frequency responses, better attenuation, more computation; lower -> smoother frequency responses, may not achieve as much attenuation
  n = 10, 
  # band edges
  w = line_noise_norm,
  # notch filter
  type = "stop"
)

# notch filtering each LFP trial
LFP_notch_filtered <- array(dim = dim(LFP_ensemble_adjusted))
for(i in 1:ncol(LFP_notch_filtered)) {
  LFP_notch_filtered[, i] <- signal::filter(filt = notch_filter, x = LFP_ensemble_adjusted[, i])
}

par(mfrow = c(2, 1))
plot(time, LFP[channel, , trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Demeaned and Clipped LFP of one Sample Trial")
plot(time, LFP_notch_filtered[, trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Notch-Filtered, Demeaned and Clipped LFP of one Sample Trial")
par(mfrow = c(1, 1))
# information may be lost

# ADF test for stationarity
adf_result_notch <- adf_stationarity_test(LFP_notch_filtered) 
adf_result_notch - adf_result_baseline[channel] # more LFPs are now stationary

# examining the autocorrelation function of some sample trials
examine_autocorr(LFP_notch_filtered) # not much difference compared to baseline

# # exclude the remaining few non-stationary trials
# p_values <- rep(NA, nrow(LFP_notch_filtered))
# for(i in 1:nrow(LFP_notch_filtered)) {
#   p_values[i] <- adf.test(LFP_notch_filtered[, i])$p.value
# }
# LFP_notch_filtered <- LFP_notch_filtered[, which(p_values < 0.05)]

# find the best filter order for each channel
adf_result_notch <- rep(NA, dim(LFP)[1])
for(k in 1:dim(LFP)[1]) {
  
  LFP_stationary_channel_k <- LFP_stationary[k, , ]
  adf_result_notch_channel_k <- rep(NA, 30)
  
  # loop over all trials
  for(j in 1:30){
    notch_filter <- signal::fir1(
      # filter order - higher -> filters with sharper frequency responses, better attenuation, more computation; lower -> smoother frequency responses, may not achieve as much attenuation
      n = j,
      # band edges
      w = line_noise_norm,
      # notch filter
      type = "stop"
    )
    # notch filtering each LFP trial
    LFP_filtered <- array(dim = dim(LFP_stationary_channel_k))
    for(i in 1:ncol(LFP_filtered)) {
      LFP_filtered[, i] <- signal::filter(filt = notch_filter, x = LFP_stationary_channel_k[, i])
    }
    # ADF test for stationarity
    adf_result_notch_channel_k[j] <- adf_stationarity_test(LFP_filtered)
  }
  adf_result_notch[k] <- which.max(adf_result_notch_channel_k)
}



### other ideas: windowing, spectral factorisation, wavelet transformations, adaptive recursive least-squares modelling



# Stationarity transformation on the whole LFP ----------------------------

transform_stationary <- function(
    LFP,
    detrending,
    normalising,
    ensemble_adjusting,
    notch_filtering,
    sampling_rate,
    print_remaining_non_stationary
    ) {
  
  # prepare matrix for the transformed LFPs
  LFP_stationary <- array(dim = c(dim(LFP)[2], dim(LFP)[3]))
  
  # # prepare list for the remaining non-stationary LFPs (before removing them)
  # stationary_trials_ind <- list()
  
  # loop over all channels
  for(i in 1:dim(LFP)[1]) {
    
    # extract the LFP from channel i
    LFP_channel_i <- LFP[i, , ]
    
    # removing linear trends
    if (detrending) {
      LFP_channel_i <- pracma::detrend(LFP_channel_i, tt = "linear")
    }
    
    # removal of temporal mean and division of temporal sd
    if (normalising) {
      LFP_channel_i <- LFP_channel_i - colMeans(LFP_channel_i)
      LFP_channel_i <- t(LFP_channel_i) / sd(t(LFP_channel_i))
      LFP_channel_i <- t(LFP_channel_i)
    }
    
    # ensemble mean removal and division of ensemble sd
    if (ensemble_adjusting) {
      LFP_channel_i <- LFP_channel_i - rowMeans(LFP_channel_i)
      LFP_channel_i <- LFP_channel_i / sd(LFP_channel_i)
    }
    
    # notch filtering each LFP trial
    if (notch_filtering) {
      
      # remove frequencies around 60 Hz (= electrical line noise)
      line_noise <- 60
      
      # Nyquist frequency: maximum frequency that can be represented in the digital signal
      nyquist_freq <- sampling_rate / 2
      line_noise_norm <- line_noise / nyquist_freq
      
      notch_filter <- signal::fir1(
        # filter order - higher -> filters with sharper frequency responses, better attenuation, more computation; lower -> smoother frequency responses, may not achieve as much attenuation
        n = 10, 
        # band edges
        w = line_noise_norm,
        # notch filter
        type = "stop"
      )
      
      # notch filtering each LFP trial
      LFP_channel_i_filtered <- array(dim = dim(LFP_channel_i))
      for(j in 1:ncol(LFP_channel_i_filtered)) {
        LFP_channel_i_filtered[, j] <- signal::filter(filt = notch_filter, x = LFP_channel_i[, j])
      }
      LFP_channel_i <- LFP_channel_i_filtered
      
    }
    
    if (print_remaining_non_stationary) {
      # display the number of remaining non-stationary LFPs
      p_values <- rep(NA, ncol(LFP_channel_i))
      for (j in 1:ncol(LFP_channel_i)) {
        p_values[j] <- adf.test(LFP_channel_i[, j])$p.value
      }
      print(paste0("Channel ", i, ": Number of LFP trials still non-stationary: ", sum(p_values > 0.05)))
    }
    
    # exclude the remaining few non-stationary trials
    # ...
    
    LFP_stationary <- abind::abind(LFP_stationary, LFP_channel_i, along = 3)
    
  }
  
  LFP_stationary <- LFP_stationary[, , 2:(dim(LFP)[1] + 1)]
  LFP_stationary <- aperm(LFP_stationary, c(3, 1, 2))
  return(LFP_stationary)
  
}

LFP_stationary <- transform_stationary(
  LFP = LFP, 
  detrending = TRUE, 
  normalising = TRUE, 
  ensemble_adjusting = TRUE, 
  notch_filtering = TRUE,
  sampling_rate = sampling_rate,
  print_remaining_non_stationary = FALSE
)

all(dim(LFP_stationary) == dim(LFP))
LFP_stationary[1:3, 1:3, 1:3]





