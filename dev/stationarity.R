
session_data <- data_preprocessing(
  session_data = session_data, 
  same_lengths = TRUE,
  ms_clipping = 10
)
list2env(session_data, envir = .GlobalEnv)
rm(list = "session_data")

# look at the 6th channel and 333 trial as an example first
channel <- 6
trial <- 333



# Check for Stationarity --------------------------------------------------

# Augmented Dickey-Fuller test (H0: non-stationary) either for all channels or a specific one
adf_stationarity_test <- function(LFP, channel) {
  
  # for the case that LFP was not cut equally
  if (is.list(LFP)) {
    
    num_channels <- unique(lengths(LFP))
    num_trials <- length(LFP)
    
    p_values <- matrix(NA, nrow = num_trials, ncol = num_channels)
    colnames(p_values) <- paste0("Channel_", 1:num_channels)

    for (j in 1:num_channels) {
      for (i in 1:num_trials) {
        p_values[i, j] <- tseries::adf.test(LFP[[i]][[j]])$p.value
      }
    }
    
    # ratio of LFPs that are stationary
    percentage_stationary <- sum(p_values < 0.05) / (num_channels * num_trials)
    
  # for the case that LFP was cut equally
  } else {
    
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
            p_values[i, j] <- tseries::adf.test(LFP[j, , i])$p.value
          }
          
        }
        
        # percentage of LFPs that are stationary
        percentage_stationary <- colSums(p_values < 0.05) / nrow(p_values)
        
        # ADF test for one specific channel
      } else if (channel %in% 1:15) {
        
        p_values <- rep(NA, dim(LFP)[3])
        
        for (i in 1:dim(LFP)[3]) {
          p_values[i] <- tseries::adf.test(LFP[channel, , i])$p.value
        }
        
        # percentage of LFPs that are stationary
        percentage_stationary <- sum(p_values < 0.05) / dim(LFP)[3]
        
      }
      
    } else if (length(dim(LFP)) == 2) {
      
      p_values <- rep(NA, dim(LFP)[2])
      
      for (i in 1:dim(LFP)[2]) {
        p_values[i] <- tseries::adf.test(LFP[, i])$p.value
      }
      
      # percentage of LFPs that are stationary
      percentage_stationary <- sum(p_values < 0.05) / dim(LFP)[2]
      
    }
    
  }
  
  return(percentage_stationary)
  
}

# percentage of LFPs that are stationary
adf_result_baseline <- adf_stationarity_test(LFP = LFP, channel = "all")
adf_result_baseline # most are stationary - but can be better in some channels

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
examine_autocorr(LFP = LFP, channel = 6)



# Dealing With Non-Stationarity - Removing Linear Trends ------------------

LFP_detrended <- pracma::detrend(LFP[channel, , ], tt = "linear")

par(mfrow = c(2, 1))
plot(time[trial, ], LFP[channel, , trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Clipped LFP of one Sample Trial")
plot(time[trial, ], LFP_detrended[, trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Detrended and Clipped LFP of one Sample Trial")
par(mfrow = c(1, 1))

# ADF test for stationarity
adf_result_detrended <- adf_stationarity_test(LFP = LFP_detrended, channel = "all") 
adf_result_detrended - adf_result_baseline[channel] # no difference to baseline

# examining the autocorrelation function of some sample trials
examine_autocorr(LFP_detrended) # no difference to baseline



# Dealing With Non-Stationarity - Removal of Temporal Mean and Div --------
# ision of Temporal SD -----------------------------------------------------

# subtract the mean of each LFP trial
LFP_demeaned <- LFP_detrended - colMeans(LFP_detrended)

# divide by the temporal sd
LFP_normalised <- t(LFP_demeaned) / sd(t(LFP_demeaned))
LFP_normalised <- t(LFP_normalised)

par(mfrow = c(2, 1))
plot(time[trial, ], LFP[channel, , trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Clipped LFP of one Sample Trial")
plot(time[trial, ], LFP_normalised[, trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Normalised, Detrended and Clipped LFP of one Sample Trial")
par(mfrow = c(1, 1))

# ADF test for stationarity
adf_result_normalised <- adf_stationarity_test(LFP_normalised, channel = "all") 
adf_result_normalised - adf_result_baseline[channel] # no difference to baseline

# examining the autocorrelation function of some sample trials
examine_autocorr(LFP_normalised) # no difference to baseline



# Dealing With Non-Stationarity - Ensemble Mean Removal, Division  --------
# of Ensemble SD ----------------------------------------------------------

# subtract the ensemble mean
LFP_ensemble_mean_removed <- LFP_normalised - rowMeans(LFP_normalised)

# divide by the ensemble sd
LFP_ensemble_adjusted <- LFP_ensemble_mean_removed / sd(LFP_ensemble_mean_removed)

par(mfrow = c(2, 1))
plot(time[trial, ], LFP[channel, , trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Clipped LFP of one Sample Trial")
plot(time[trial, ], LFP_ensemble_adjusted[, trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Normalised, Detrended and Clipped LFP of one Sample Trial, Ensemble Mean Removed, Divided by the Ensemble SD")
par(mfrow = c(1, 1))

# ADF test for stationarity
adf_result_detrended <- adf_stationarity_test(LFP_ensemble_adjusted, channel = "all") 
adf_result_detrended - adf_result_baseline[channel] # only slightly better

# examining the autocorrelation function of some sample trials
examine_autocorr(LFP_ensemble_adjusted) # almost no difference to baseline



# Dealing With Non-Stationarity - Notch Filtering -------------------------

# goal: removing line noise using a notch filter

# remove frequencies around 60 Hz (= electrical line noise)
line_noise <- 60
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
plot(time[trial, ], LFP[channel, , trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Demeaned and Clipped LFP of one Sample Trial")
plot(time[trial, ], LFP_notch_filtered[, trial], type = "l", xlab = "Time", ylab = "Amplitude", main = "Notch-Filtered, Demeaned and Clipped LFP of one Sample Trial")
par(mfrow = c(1, 1))

# ADF test for stationarity
adf_result_notch <- adf_stationarity_test(LFP_notch_filtered, channel = "all") 
adf_result_notch - adf_result_baseline[channel] # more LFPs are now stationary

# examining the autocorrelation function of some sample trials
examine_autocorr(LFP_notch_filtered) # not much difference compared to baseline

# find the best filter order for each channel (first run transform_stationary() with notch_filtering = FALSE)
find_best_filter_order <- function(LFP){
  
  adf_result_notch <- rep(NA, dim(LFP)[1])
  for(k in 1:dim(LFP)[1]) {
    
    LFP_stationary_channel_k <- LFP[k, , ]
    adf_result_notch_channel_k <- rep(NA, 30)
    
    # try different values for the filter order
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
      adf_result_notch_channel_k[j] <- adf_stationarity_test(LFP_filtered, channel = "all") # shouldn't channel be set to k???????
    }
    # get the filter order for which the most trials are stationary
    adf_result_notch[k] <- which.max(adf_result_notch_channel_k)
  }
  
  return(adf_result_notch)
  
}

adf_result_notch <- find_best_filter_order(LFP = LFP_stationary)



# Stationarity Transformation on the Whole LFP ----------------------------

transform_stationary <- function(
    LFP,
    detrending,
    normalising,
    ensemble_adjusting,
    notch_filtering,
    notch_filter_order,
    sampling_rate,
    print_remaining_non_stationary
) {
  
  num_channels <- unique(lengths(LFP))
  num_trials <- length(LFP)

  LFP_stationary <- vector("list", length = num_trials)

  for(i in 1:num_channels) {
    
    # extract the LFP from channel i
    LFP_channel_i <- lapply(LFP, function(x) x[[i]]) # trial, time
    
    # removing linear trends
    if (detrending) {
      # detrending each trial
      for (j in 1:num_trials) {
        LFP_channel_i[[j]] <- pracma::detrend(LFP_channel_i[[j]], tt = "linear")
      }
    }
    
    # removal of temporal mean and division of temporal sd
    if (normalising) {
      trial_means <- lapply(LFP_channel_i, mean)
      trial_sds <- lapply(LFP_channel_i, sd)
      for (j in 1:num_trials) {
        LFP_channel_i[[j]] <- LFP_channel_i[[j]] - trial_means[[j]]
        LFP_channel_i[[j]] <- LFP_channel_i[[j]] / trial_sds[[j]]
      }
    }
    
    # ensemble mean removal and division of ensemble sd
    if (ensemble_adjusting) {

      # transform LFP_channel_i into a matrix, fill the rest with NAs
      max_len <- max(sapply(LFP_channel_i, length))
      LFP_channel_i_mat <- matrix(NA, nrow = length(LFP_channel_i), ncol = max_len)
      for (j in seq_along(LFP_channel_i)) {
        LFP_channel_i_mat[j, 1:length(LFP_channel_i[[j]])] <- LFP_channel_i[[j]]
      } # trials, time
      
      ensemble_means <- colMeans(LFP_channel_i_mat, na.rm = TRUE)
      ensenlble_sds <- apply(LFP_channel_i_mat, 2, sd, na.rm = TRUE)
      
      # handle columns with zero standard deviation by setting sd to 1 to avoid division by zero
      ensenlble_sds[ensenlble_sds == 0] <- 1
      
      # substract the ensemble mean and divide by ensemble sd
      LFP_channel_i_mat <- t(t(LFP_channel_i_mat) - ensemble_means)
      LFP_channel_i_mat <- t(t(LFP_channel_i_mat) / ensenlble_sds)
      
      # convert back to a list
      LFP_channel_i <- vector("list", nrow(LFP_channel_i_mat))
      for (j in seq_len(nrow(LFP_channel_i_mat))) {
        LFP_channel_i[[j]] <- LFP_channel_i_mat[
          j, 
          !is.na(LFP_channel_i_mat[j, ]) & !is.nan(LFP_channel_i_mat[j, ]) & is.finite(LFP_channel_i_mat[j, ])
        ]
      }

    }
    
    # notch filtering each LFP trial
    if (notch_filtering) {
      
      # remove frequencies around 60 Hz (= electrical line noise)
      line_noise <- 60
      line_noise_norm <- line_noise / nyquist_freq
      
      # notch filtering each LFP trial
      LFP_channel_i_filtered <- list()
      for (j in 1:num_trials) {
        
        if (notch_filter_order[j] == 0) {
          # no notch filter
          LFP_channel_i_filtered[[j]] <- LFP_channel_i[[j]]
        } else {
          # filter
          notch_filter <- signal::fir1(
            # filter order - higher -> filters with sharper frequency responses, better attenuation, more computation; lower -> smoother frequency responses, may not achieve as much attenuation
            n = notch_filter_order[j], 
            # band edges
            w = line_noise_norm,
            # notch filter
            type = "stop"
          )
          LFP_channel_i_filtered[[j]] <- signal::filter(filt = notch_filter, x = LFP_channel_i[[j]])
        }
        
      }
      LFP_channel_i <- LFP_channel_i_filtered
      
    }
    
    if (print_remaining_non_stationary) {
      # display the number of remaining non-stationary LFPs
      p_values <- rep(NA, length(LFP_channel_i))
      for (j in 1:1363) {
        p_values[j] <- tseries::adf.test(LFP_channel_i[[j]])$p.value
      }
      print(paste0("Channel ", i, ": Number of LFP trials still non-stationary (before removing all): ", sum(p_values > 0.05)))
    }
    
    # add LFP_channel_i to the list
    for (j in 1:num_trials) {
      LFP_stationary[[j]][[i]] <- LFP_channel_i[[j]]
    }
    
  }
  
  return(LFP_stationary)
  
}

session_data <- load_data(session = session)
session_data <- data_preprocessing(
  session_data = session_data, 
  same_lengths = FALSE,
  ms_clipping = 10
)
list2env(session_data, envir = .GlobalEnv)
rm(list = "session_data")

LFP_stationary <- transform_stationary(
  LFP = LFP, 
  detrending = TRUE, 
  normalising = TRUE, 
  ensemble_adjusting = FALSE,
  notch_filtering = FALSE,
  notch_filter_order = NULL,
  sampling_rate = sampling_rate,
  print_remaining_non_stationary = FALSE
)

adf_stationarity_test(LFP = LFP_stationary)
adf_stationarity_test(LFP = LFP)

# find the best filter order for each trial (first run transform_stationary() with notch_filtering = FALSE)
find_best_filter_order <- function(LFP, max_filter_order, nyquist_freq){
  
  num_channels <- unique(lengths(LFP))
  num_trials <- length(LFP)
  
  adf_result_notch <- rep(NA, num_trials)
  
  # loop over all trials
  for(k in 1:num_trials) {
    
    # extract the LFP from trial k
    LFP_trial_k <- LFP[[k]] # channel, time
    
    adf_result_notch_trial_k <- rep(NA, max_filter_order + 1)
    
    # remove frequencies around 60 Hz (= electrical line noise)
    line_noise <- 60
    line_noise_norm <- line_noise / nyquist_freq
    
    # no filtering
    LFP_not_filtered <- LFP_trial_k
    
    # ADF test for stationarity
    p_values <- rep(NA, length(LFP_not_filtered))
    for (i in 1:length(LFP_not_filtered)) {
      p_values[i] <- tseries::adf.test(LFP_not_filtered[[i]])$p.value
    }
    
    adf_result_notch_trial_k[1] <- sum(p_values < 0.05) / length(p_values)
    
    # try different values for the filter order
    for(j in 1:max_filter_order){
      notch_filter <- signal::fir1(
        # filter order - higher -> filters with sharper frequency responses, better attenuation, more computation; lower -> smoother frequency responses, may not achieve as much attenuation
        n = j,
        # band edges
        w = line_noise_norm,
        # notch filter
        type = "stop"
      )
      
      # notch filtering each LFP signal of each channel
      LFP_filtered <- list()
      for (i in 1:num_channels) {
        LFP_filtered[[i]] <- signal::filter(filt = notch_filter, x = LFP_trial_k[[i]])
      }
      
      # ADF test for stationarity
      p_values <- rep(NA, length(LFP_filtered))
      for (i in 1:length(LFP_filtered)) {
        p_values[i] <- tseries::adf.test(LFP_filtered[[i]])$p.value
      }
      adf_result_notch_trial_k[j + 1] <- sum(p_values < 0.05) / length(p_values)
    }
    
    # get the filter order for each trial for which the most LFP channels are stationary
    adf_result_notch[k] <- which.max(adf_result_notch_trial_k) - 1
    
  }
  
  return(adf_result_notch)
  
}

adf_result_notch <- find_best_filter_order(
  LFP = LFP_stationary,
  max_filter_order = 20, 
  nyquist_freq = nyquist_freq
)
head(adf_result_notch)
table(adf_result_notch)

LFP_stationary <- transform_stationary(
  LFP = LFP, 
  detrending = TRUE, 
  normalising = TRUE, 
  ensemble_adjusting = FALSE,
  notch_filtering = TRUE,
  notch_filter_order = adf_result_notch,
  sampling_rate = sampling_rate,
  print_remaining_non_stationary = FALSE
)

adf_stationarity_test(LFP = LFP_stationary)
adf_result_baseline

# exclude the remaining few non-stationary trials (not applied in analysis)
exclude_remaining_non_stationary <- function(LFP) {
  
  num_trials <- length(LFP)
  num_channels <- unique(lengths(LFP))
  
  p_values <- matrix(NA, nrow = num_trials, ncol = num_channels)
  
  for (i in 1:num_trials) {
    for (j in 1:num_channels) {
      # ADF test
      p_values[i, j] <- tseries::adf.test(LFP[[i]][[j]])$p.value
    }
  }
  
  all_channels_stationary <- rep(NA, nrow(p_values))
  
  for (i in 1:nrow(p_values)) {
    
    # check if trials in all channels are stationary
    all_channels_stationary[i] <- all(p_values[i, ] < 0.05)
    stationary_ind <- which(all_channels_stationary)
    
  }
  
  num_removed <- length(LFP) - length(LFP[stationary_ind])
  print(paste(num_removed, "non-stationary trials removed."))
  
  # filter the stationary trials
  LFP_stationary <- LFP[stationary_ind]
  
  return(LFP_stationary)
  
}

LFP_stationary_hra <- exclude_remaining_non_stationary(LFP = LFP_stationary)

adf_stationarity_test(LFP = LFP_stationary_hra)

# adapt unprimed and primed ind to removed trials
unprimed_ind_hra <- unprimed_ind[which(unprimed_ind %in% c(1:dim(LFP_stationary_hra)[3]))]
primed_ind_hra <- primed_ind[which(primed_ind %in% c(1:dim(LFP_stationary_hra)[3]))]
length(unprimed_ind_hra)
length(primed_ind_hra)



# Dealing With Non-Stationarity - Differencing ----------------------------

# problem with differencing: change of interpretation - causal interactions are now among changes in the LFPs rather then in the LFPs per se (2) 
# (not applied in analysis)

trial <- 15

LFP_sample_trial <- LFP_stationary[[trial]]
LFP_sample_trial <- do.call(cbind, lapply(LFP_sample_trial, function(x) as.numeric(x)))
dim(LFP_sample_trial) # trial, channel

# differencing one time
LFP_differenced <- diff(x = LFP_sample_trial, differences = 1)

par(mfrow = c(2, 1))
plot(LFP_sample_trial[, 5], type = "l", xlab = "", ylab = "", main = "Demeaned and Notch Filtered LFP of one Sample Trial")
plot(LFP_differenced[, 5], type = "l", xlab = "", ylab = "", main = "Demeaned, Notch Filtered and Differenced LFP of one Sample Trial")
par(mfrow = c(1, 1))

# ADF test for stationarity
tseries::adf.test(LFP_sample_trial[, 12])$p.value # not stationary
tseries::adf.test(LFP_differenced[, 12])$p.value # now stationary

# examining the autocorrelation function
examine_autocorr(LFP_sample_trial)
examine_autocorr(LFP_differenced) # better

# in a function for all trials
differencing <- function(LFP, order) {
  
  num_trials <- length(LFP)
  LFP_differenced <- vector("list", num_trials)
  
  for (i in 1:num_trials) {
    
    # difference each trial
    LFP_trial_i <- LFP[[i]]
    LFP_trial_i <- do.call(cbind, lapply(LFP_trial_i, function(x) as.numeric(x)))
    LFP_trial_i_differenced_mat <- diff(x = LFP_trial_i, differences = order)
    
    # transform back into a list
    col_lists <- apply(LFP_trial_i_differenced_mat, 2, function(x) as.list(x))
    LFP_differenced_trial_i <- lapply(1:ncol(LFP_trial_i_differenced_mat), function(i) col_lists[[i]])
    LFP_differenced[[i]] <- LFP_differenced_trial_i
    
  }
  
  return(LFP_differenced)
  
}

LFP_differenced <- differencing(LFP = LFP_stationary, order = 1)

# check stationarity
adf_stationarity_test(LFP = LFP_stationary)
adf_stationarity_test(LFP = LFP_differenced)
