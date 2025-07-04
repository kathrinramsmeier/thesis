
library(tidyverse)  # data manipulation
library(vars)       # VAR model
library(bruceR)     # Granger causality
library(ggplot2)    # plots

options(scipen = 1000)



# Load Session Data -------------------------------------------------------

load_data <- function(session) {
  
  # detach("package:reticulate")
  library(reticulate)
  np <- import("numpy")
  
  setwd("~\\data\\")
  
  # LFP
  LFP <- np$load(paste0(session, "-npy\\lfp_array_uv.npy"))
  
  # electrode channel
  probe <- read.csv(paste0(session, "-npy\\probe.csv"), header = FALSE)
  probe_header <- read.csv(paste0(session, "-npy\\probe_header.csv"), header = FALSE)
  colnames(probe) <- probe_header
  rm(list = "probe_header")
  
  # time
  time <- np$load(paste0(session, "-npy\\time_array_ms.npy"))
  time <- time[1, ]
  
  # behaviour
  behaviour <- np$load(paste0(session, "-npy\\behavior.npy"))
  behaviour_header <- read.csv(paste0(session, "-npy\\behavior_header.csv"), header = FALSE)
  colnames(behaviour) <- behaviour_header
  behaviour <- as.data.frame(behaviour)
  rm(list = "behaviour_header")
  
  # task
  task <- read.csv(paste0(session, "-npy\\task.csv"), header = FALSE)
  task_header <- read.csv(paste0(session, "-npy\\task_header.csv"), header = FALSE)
  colnames(task) <- task_header
  rm(list = "task_header")
  
  # recordings info
  recordinginfo <- read.csv(paste0(session, "-npy\\recordinginfo.csv"), header = FALSE)
  recordinginfo_header <- read.csv(paste0(session, "-npy\\recordinginfo_header.csv"), header = FALSE)
  colnames(recordinginfo) <- recordinginfo_header
  rm(list = "recordinginfo_header")
  
  # probe
  probe <- read.csv(paste0(session, "-npy\\probe.csv"), header = FALSE)
  probe_header <- read.csv(paste0(session, "-npy\\probe_header.csv"), header = FALSE)
  colnames(probe) <- probe_header
  rm(list = "probe_header")
  
  return(list(
    LFP = LFP, 
    probe = probe, 
    time = time,
    behaviour = behaviour,
    task = task,
    recordinginfo = recordinginfo,
    probe = probe
  ))
  
}



# Data Preprocessing ------------------------------------------------------

data_preprocessing <- function(session_data, same_lengths, ms_clipping) {
  
  # add time value at the start and at the end
  session_data$time <- c(
    session_data$time[1] - diff(session_data$time[c(1, 2)]), 
    session_data$time, 
    rev(session_data$time)[1] + diff(rev(session_data$time)[c(1, 2)])
  )
  
  # exclude incorrect trials and catch trials (where the colours were the same) from LFP and task
  correct_trials_ind <- which(session_data$behaviour$accuracy_logical == 1)
  non_catch_trails_ind <- session_data$task$trial_number_count[(session_data$task$catch_trial_logical == 0)]
  trials_ind <- intersect(correct_trials_ind, non_catch_trails_ind)
  session_data$LFP <- session_data$LFP[, , trials_ind]
  session_data$task <- session_data$task[trials_ind, ] 
  session_data$task$trial_number_count <- 1:nrow(session_data$task)
  session_data$behaviour <- session_data$behaviour[trials_ind, ]
  
  # exclude ms_clipping ms before the reaction
  clipped_reaction <- session_data$behaviour$reaction_time_ms - ms_clipping
  
  # for clipping: get the ind in the time vector that is closest to the reaction time
  clipped_reaction_time_ind <- rep(NA, dim(session_data$LFP)[3])
  for (i in 1:dim(session_data$LFP)[3]) {
    clipped_reaction_time_ind[i] <- which.min(abs(session_data$time - clipped_reaction[i]))
  }
  
  if (same_lengths) {
    
    # get the times from start to the clipped reaction time and pad LFP signals to the same lengths
    time_ind <- 1:min(clipped_reaction_time_ind)
    time_i <- session_data$time[time_ind]
    
    # only keep LFP after stimulus onset
    time_ind <- which(time_i > 0)
    
    # now filter the LFP
    time <- matrix(nrow = dim(session_data$LFP)[3], ncol = length(time_ind))
    LFP <- array(dim = c(dim(session_data$LFP)[1], length(time_ind), dim(session_data$LFP)[3]))
    for (i in 1:dim(session_data$LFP)[3]) {
      time[i, ] <- time_i[time_ind]
      for (j in 1:dim(LFP)[1]) {
        LFP[j, , i] <- session_data$LFP[j, time_ind, i]
      }
    }
    
  } else {
    
    time <- list()
    LFP <- vector("list", length = dim(session_data$LFP)[3])
    
    for (i in 1:dim(session_data$LFP)[3]) {
      
      # get the times from start to the clipped reaction time for each trial
      time_ind <- 1:clipped_reaction_time_ind[i]
      time_i <- session_data$time[time_ind]
      
      # only keep LFP after stimulus onset
      time_ind <- which(time_i > 0)
      time[[i]] <- time_i[time_ind]
      
      # loop over all channels
      for (j in 1:15) {
        LFP[[i]][[j]] <- session_data$LFP[j, time_ind, i]
      }
    }
    # LFP: trials, channels, time
    
  }
  
  session_data$time <- time
  session_data$LFP <- LFP
  rm(list = c("time", "LFP"))
  
  if (same_lengths) {
    
    # checks
    stopifnot(
      dim(session_data$LFP)[1] == 15 &
        dim(session_data$LFP)[2] == dim(session_data$time)[2] &
        dim(session_data$LFP)[3] == nrow(session_data$task)
    )
    
  }
  
  return(session_data)
  
}



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



# Stationarity Transformation ---------------------------------------------

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

# differencing
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

# exclude the remaining few non-stationary trials
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




# Power of LFP for Bands, Channels and Priming Status ---------------------

calculate_bands_power <- function(LFP, electrode_channel, un_primed_ind, hanning_windowing, frequency_bands) {
  
  # filter either primed or unprimed trials 
  LFP_trials <- LFP[un_primed_ind]
  
  # filter the electrode channel
  LFP_trials <- lapply(LFP_trials, function(x) x[[electrode_channel]]) # trial, time
  
  num_trials <- length(LFP_trials)
  
  # prepare matrix for results
  bands_power <- matrix(nrow = num_trials, ncol = length(frequency_bands))
  colnames(bands_power) <- names(frequency_bands)
  rownames(bands_power) <- paste0("trial_", 1:num_trials)
  
  # iterate over all trials
  for(i in 1:num_trials) {
    
    trial <- LFP_trials[[i]]
    n <- length(trial)
    
    if(hanning_windowing) {
      # Hanning windowing
      trial <- trial * gsignal::hanning(n = 30, method = "periodic")
    }
    
    # Fourier transform to get power
    LFP_fftransformed <- fft(trial)
    
    # frequency bins
    freq_resolution <- sampling_rate / n
    freq_bins <- seq(0, nyquist_freq, by = freq_resolution)
    
    # calculate the noramlised power spectrum (for the first half of the spectrum) (normalising = dividing by n)
    power <- abs(LFP_fftransformed[1:length(freq_bins)])^2 / n
    
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

# mean bands power in every channel
calculate_mean_bands_power <- function(LFP, un_primed_ind, hanning_windowing, frequency_bands) {
  
  num_channels <- unique(lengths(LFP))
  
  # prepare matrix to store the results
  rel_bands_power_un_primed <- matrix(nrow = num_channels, ncol = length(frequency_bands))
  colnames(rel_bands_power_un_primed) <- names(frequency_bands)
  rownames(rel_bands_power_un_primed) <- 1:num_channels
  
  # take the mean over all trials and loop over all channels (mean bands power in every channel)
  for(i in 1:num_channels){
    rel_bands_power_un_primed_ch <- calculate_bands_power(
      LFP = LFP,
      frequency_bands = frequency_bands, 
      un_primed_ind = un_primed_ind,
      hanning_windowing = hanning_windowing,
      electrode_channel = i
    )
    rel_bands_power_un_primed[i, ] <- apply(rel_bands_power_un_primed_ch, 2, mean)
  }
  
  return(rel_bands_power_un_primed)
  
}

# plot power in each channel for primed and unprimed trials for each frequency band
plot_power <- function(mean_bands_power_unprimed, mean_bands_power_primed, frequency_bands) {
  
  min_val <- min(min(mean_bands_power_unprimed), min(mean_bands_power_primed))
  max_val <- max(max(mean_bands_power_unprimed), max(mean_bands_power_primed))
  
  library(gridExtra)
  
  plots <- list()
  
  fb_len <- length(frequency_bands)
  
  plot_names <- rep(NA, fb_len)
  for (i in 1:fb_len) {
    plot_names[i] <- paste0(names(frequency_bands[i]), " (", frequency_bands[[i]][1], " - ", frequency_bands[[i]][2], " Hz)")
  }
  plot_names[fb_len] <- "gamma (> 30 Hz)"
  
  for (i in 1:fb_len) {
    plots[[i]] <- ggplot(aes(x = !!as.name(names(frequency_bands)[i]), y = 1:15), data = as.data.frame(mean_bands_power_primed)) +
      geom_point(size = 3, shape = 16, colour = "darkgoldenrod4") +
      geom_point(aes(x = !!as.name(names(frequency_bands)[i]), y = 1:15), data = as.data.frame(mean_bands_power_unprimed), size = 3, shape = 17, colour = "darkgoldenrod2") +
      scale_y_continuous(breaks = seq(1, 15, 1)) +
      scale_y_reverse(breaks = seq(1, 15, 1)) +
      xlim(min_val, max_val) +
      labs(x = "", y = "", title = plot_names[i]) +
      theme_minimal() 
  }
  
  plot <- grid.arrange(grobs = plots, ncol = 4)
  print(plot)
  
}



# Coherence ---------------------------------------------------------------

# filter the LFP in the different frequency bands
freqband_filtering <- function(LFP_trial, frequency_bands, sampling_rate, filter_order) {
  
  LFP_freqband_filtered <- matrix(nrow = length(LFP_trial), ncol = length(frequency_bands))
  colnames(LFP_freqband_filtered) <- names(frequency_bands)
  
  nyquist_freq <- sampling_rate / 2
  
  # avoid numerical problem (only needed when all frequency bands are examined)
  # frequency_bands[[1]][1] <- 0.01
  
  # loop over all frequency bands
  for (i in 1:length(frequency_bands)) {
    
    freq_band <- c(frequency_bands[[i]][1], frequency_bands[[i]][2])
    
    # butterworth filtering
    butterworth <- gsignal::butter(n = filter_order, w = freq_band / nyquist_freq, type = "pass")
    
    # forward and reverse filtering (correction for phase distortion)
    LFP_freqband_filtered[, i] <- gsignal::filtfilt(filt = butterworth$b, a = butterworth$a, x = LFP_trial)
    
  }
  
  return(LFP_freqband_filtered)
  
}

# calculate coherence between two signals
calculate_coherence <- function(LFP_trial_1_freqband_filtered, LFP_trial_2_freqband_filtered, meaned, frequency_bands) {
  
  coherences <- list()
  
  # loop over all frequency bands
  for (i in 1:ncol(LFP_trial_1_freqband_filtered)) {
    
    # cross-power spectral magnitudes and density autospectral densities
    # formulas from: https://www.psychologie.hhu.de/fileadmin/redaktion/Fakultaeten/Mathematisch-Naturwissenschaftliche_Fakultaet/Psychologie/CompPsy/Papers/vinck2010_wingerden.pdf
    csd <- gsignal::cpsd(cbind(LFP_trial_1_freqband_filtered[, i], LFP_trial_2_freqband_filtered[, i]))$cross
    csd_magnitudes <- abs(csd)
    asd_1 <- gsignal::cpsd(cbind(LFP_trial_1_freqband_filtered[, i], LFP_trial_1_freqband_filtered[, i]))$cross
    asd_2 <- gsignal::cpsd(cbind(LFP_trial_2_freqband_filtered[, i], LFP_trial_2_freqband_filtered[, i]))$cross
    
    # coherence
    coherences[[i]] <- as.numeric(csd_magnitudes^2 / (asd_1 * asd_2))
    
  }
  
  # mean coherence in each frequency band
  if(meaned) {
    coherences <- unlist(lapply(coherences, mean))
    names(coherences) <- names(frequency_bands)
  }
  
  return(coherences)
  
}

# compute coherence between the channels in each frequency band for all primed and unprimed trials
calculate_avg_coherence <- function(LFP, un_primed_ind, frequency_bands, sampling_rate, filter_order) {
  
  # filter either primed or unprimed trials
  LFP_trials <- LFP[un_primed_ind]
  
  num_channels <- unique(lengths(LFP_trials))
  num_trials <- length(LFP_trials)
  
  # prepare matrices to store results
  coherences <- array(dim = c(num_channels, num_channels, length(frequency_bands)))
  sum_coherences <- array(0, dim = c(num_channels, num_channels, length(frequency_bands)))
  
  # loop over all trials
  for (k in 1:num_trials) {
    
    # coherence for trial i between each channel
    for (i in 1:num_channels) {
      for (j in 1:num_channels) {
        trial_freqband_filtered_1 <- freqband_filtering(
          LFP_trial = as.numeric(LFP_trials[[k]][[i]]), 
          frequency_bands = frequency_bands, 
          sampling_rate = sampling_rate,
          filter_order = filter_order
        )
        trial_freqband_filtered_2 <- freqband_filtering(
          LFP_trial = as.numeric(LFP_trials[[k]][[j]]), 
          frequency_bands = frequency_bands, 
          sampling_rate = sampling_rate,
          filter_order = filter_order
        )
        coherences[i, j, ] <- calculate_coherence(
          LFP_trial_1_freqband_filtered = trial_freqband_filtered_1, 
          LFP_trial_2_freqband_filtered = trial_freqband_filtered_2,
          meaned = TRUE,
          frequency_bands = frequency_bands
        )
      }
    }
    
    sum_coherences <- sum_coherences + coherences
    
    print(paste("Trial", k, "of", num_trials, "completed."))
    
  }
  
  # calculate mean coherence over all trials
  avg_coherences <- sum_coherences / num_trials
  return(avg_coherences)
  
}



# PLV, PLI and PPC --------------------------------------------------------

# compute PLV and PLI across the frequency bands
# formulas from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3674231/#:~:text=PLV%20can%20therefore%20be%20viewed,each%20scaling%20of%20the%20wavelet.
calculate_PLV_PLI_PPC <- function(LFP_trial_1_freqband_filtered, LFP_trial_2_freqband_filtered, method) {
  
  stopifnot(method %in% c("PLV", "PLI", "PPC"))
  
  num_freq_bands <- dim(LFP_trial_1_freqband_filtered)[2]
  
  result <- rep(NA, num_freq_bands)
  names(result) <- colnames(LFP_trial_1_freqband_filtered)
  
  # loop over all frequency bands
  for (i in 1:num_freq_bands) {
    
    # Hilbert transform to get the analytical signal
    LFP_1_freqband_ht <- gsignal::hilbert(LFP_trial_1_freqband_filtered[, i])
    LFP_2_freqband_ht <- gsignal::hilbert(LFP_trial_2_freqband_filtered[, i])
    
    # compute relative phase = phase difference
    phase_diffs <- pracma::angle(LFP_1_freqband_ht * Conj(LFP_2_freqband_ht))  # resembles arg((LFP_1_band_ht * Conj(LFP_2_band_ht)) / (abs(LFP_1_band_ht) * abs(LFP_2_band_ht)))
    
    if (method == "PLV") {
      
      # compute the PLV for each frequency band
      result[i] <- abs(mean(exp(1i * phase_diffs)))
      
    } else if (method == "PLI") {
      
      # compute the PLI for each frequency band
      result[i] <- abs(mean(sign(phase_diffs)))
      
    } else if (method == "PPC") {
      
      N <- length(phase_diffs)
      ppc_sum <- 0
      for (j in 1:(N-1)) {
        for (k in (j+1):N) {
          ppc_sum <- ppc_sum + (cos(phase_diffs[j]) * cos(phase_diffs[k]) + sin(phase_diffs[j]) * sin(phase_diffs[k]))
        }
      }
      result[i] <- (2 / (N * (N - 1))) * ppc_sum
      
    }
  }
  
  return(result)
  
}

# compute PLV/PLI/PPC between the channels in each frequency band for all primed and unprimed trials
calculate_PLV_PLI_PPC_hat <- function(LFP, un_primed_ind, frequency_bands, sampling_rate, filter_order, method) {
  
  # filter either primed or unprimed trials
  LFP_trials <- LFP[un_primed_ind]
  
  num_channels <- unique(lengths(LFP_trials))
  num_trials <- length(LFP_trials)
  num_freq_bands <- length(frequency_bands)
  
  # prepare matrices to store results
  sum_results <- array(0 ,dim = c(num_channels, num_channels, num_freq_bands))
  
  # loop over all trials
  for (k in 1:num_trials) {
    
    results <- array(dim = c(num_channels, num_channels, num_freq_bands))
    
    # PLV/PLI/PPC for trial i between each channel
    for (i in 1:num_channels) {
      for (j in 1:num_channels) {
        trial_freqband_filtered_1 <- freqband_filtering(
          LFP_trial = as.numeric(LFP[[k]][[i]]), 
          frequency_bands = frequency_bands, 
          sampling_rate = sampling_rate,
          filter_order = filter_order
        )
        trial_freqband_filtered_2 <- freqband_filtering(
          LFP_trial = as.numeric(LFP[[k]][[j]]), 
          frequency_bands = frequency_bands, 
          sampling_rate = sampling_rate,
          filter_order = filter_order
        )
        results[i, j, ] <- calculate_PLV_PLI_PPC(
          LFP_trial_1_freqband_filtered = trial_freqband_filtered_1, 
          LFP_trial_2_freqband_filtered = trial_freqband_filtered_2,
          method = method
        )
      }
    }
    
    sum_results <- sum_results + results
    
    print(paste("Trial", k, "of", num_trials, "completed."))
    
  }
  
  # calculate mean PLV/PLI/PPC over all trials 
  # see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3674231/#:~:text=PLV%20can%20therefore%20be%20viewed,each%20scaling%20of%20the%20wavelet.
  avg_results <- sum_results / length(LFP_trials)
  
  return(avg_results)
  
}

plot_heatmap <- function(result_array, frequency_bands, label_name) {
  
  len_fb <- length(frequency_bands)
  
  # preparation of the data
  coherence_df <- as.data.frame(as.table(result_array))
  colnames(coherence_df) <- c("channel1", "channel2", "frequency_band", "Coherence")
  coherence_df$channel1 <- as.numeric(as.factor(coherence_df$channel1))
  coherence_df$channel2 <- as.numeric(as.factor(coherence_df$channel2))
  coherence_df$frequency_band <- factor(coherence_df$frequency_band)
  
  plot_names <- rep(NA, len_fb)
  for (i in 1:len_fb) {
    plot_names[i] <- paste0(names(frequency_bands[i]), " (", frequency_bands[[i]][1], " - ", frequency_bands[[i]][2], " Hz)")
  }
  levels(coherence_df$frequency_band) <- plot_names
  
  # heatmap plots for all frequency bands
  p <- ggplot(coherence_df, aes(x = channel1, y = channel2, fill = Coherence)) +
    geom_tile() +
    scale_fill_gradientn(colors = c("white", "#ADD8E6", "#3948A9")) +
    facet_wrap(~ frequency_band, ncol = len_fb) +
    labs(x = "", y = "", fill = label_name) +
    theme_minimal() +
    scale_x_continuous(breaks = 1:15, labels = as.character(1:15), position = "top") + 
    scale_y_reverse(breaks = 1:15, labels = as.character(1:15)) +  
    theme(
      axis.title = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      plot.margin = margin(10, 10, 10, 10),
      panel.spacing = unit(0, "lines")
    )
  
  print(p)
  
}

# plot the difference of the average coherence between the channels for the different frequency bands between unprimed and primed trials as heatmaps
plot_diff_heatmap <- function(result_array_unprimed, results_array_primed, frequency_bands, label_name) {
  
  len_fb <- length(frequency_bands)
  
  # preparation of the data
  coherence_df <- as.data.frame(as.table(results_array_primed - result_array_unprimed))
  colnames(coherence_df) <- c("channel1", "channel2", "frequency_band", "Delta_Coherence")
  coherence_df$channel1 <- as.numeric(as.factor(coherence_df$channel1))
  coherence_df$channel2 <- as.numeric(as.factor(coherence_df$channel2))
  coherence_df$frequency_band <- factor(coherence_df$frequency_band)
  
  plot_names <- rep(NA, len_fb)
  for (i in 1:len_fb) {
    plot_names[i] <- paste0(names(frequency_bands[i]), " (", frequency_bands[[i]][1], " - ", frequency_bands[[i]][2], " Hz)")
  }
  levels(coherence_df$frequency_band) <- plot_names
  
  # heatmap plots for all frequency bands
  p <- ggplot(coherence_df, aes(x = channel1, y = channel2, fill = Delta_Coherence)) +
    geom_tile() +
    scale_fill_gradientn(colors = c("orange", "white", "#3948A9")) +
    facet_wrap(~ frequency_band, ncol = len_fb) +
    labs(x = "", y = "", fill = label_name) +
    theme_minimal() +
    scale_x_continuous(breaks = 1:15, labels = as.character(1:15)) +
    scale_y_reverse(breaks = 1:15, labels = as.character(1:15)) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      plot.margin = margin(10, 10, 10, 10),
      panel.spacing = unit(0, "lines")
    )
  
  print(p)
  
}



# Granger Causality Analysis ----------------------------------------------

# from Seth et. al.

# 1.) note the amount of variance accounted for by the VAR model in terms of the adjusted sum-square-error
calculate_adj_sum_square_err <- function(VAR_model) {
  
  residuals <- residuals(VAR_model)
  lags <- VAR_model$p
  
  # Calculate residual sum of squares for each variable
  RSS <- apply(residuals, 1, function(x) sum(x^2))
  
  n_obs <- dim(VAR_model$y)[1]
  n_var <- dim(VAR_model$y)[2]
  
  df_error <- (n_obs - lags) - (n_var * lags)
  df_total <- (n_obs - lags)
  
  rss_adj <- numeric(n_var)
  
  for (i in 1:n_var) {
    xvec <- VAR_model$y[(lags + 1):n_obs, i]
    rss2 <- sum(xvec^2)
    rss_adj[i] <- 1 - ((RSS[i] / df_error) / (rss2 / df_total))
  }
  
  return(rss_adj)
  
}

# 2.) checking the VAR model's consistency (Code from the matlab GCCA toolbox converted into R)
check_consistency <- function(VAR_model) {
  
  VAR_model_values <- VAR_model$datamat[, 1:15]
  residuals <- residuals(VAR_model)
  
  predictions <- as.matrix(as.matrix(VAR_model_values) - as.matrix(residuals))
  
  VAR_model_values <- as.matrix(VAR_model_values)
  
  # covariance estimates
  R_r <- VAR_model_values %*% t(VAR_model_values)  
  R_s <- predictions %*% t(predictions)   
  
  # compare matrix norms
  cons <- 1 - norm(R_s - R_r, type = "F") / norm(R_r, type = "F") 
  return(cons)
  
}

# 3.) Durbin-Watson statistic
dW_VAR_test <- function(VAR_model) {
  d <- rep(NA, 15)
  for (i in 1:15) {
    string <- paste0("VAR_model$varresult$y", i)
    coeff <- eval(parse(text = string))
    d[i] <- lmtest::dwtest(coeff)$statistic
  }
  return(d)
}

# 4.) Portmanteau- and Breusch-Godfrey test
# H0: no autocorrelation in the residuals
# Breusch, T . S. (1978), Testing for autocorrelation in dynamic linear models, Australian Economic Papers, 17: 334-355.
# Godfrey, L. G. (1978), Testing for higher order serial correlation in regression equations when the regressors include lagged dependent variables, Econometrica, 46: 1303-1313.
P_BG_test <- function(VAR_model) {
  p_values <- rep(NA, 2)
  names(p_values) <- c("Portmanteu", "Breusch-Godfrey")
  p_values[1] <- vars::serial.test(VAR_model_sample_trial, lags.pt = 2, type = "PT.adjusted")$serial$p.value
  p_values[2] <- vars::serial.test(VAR_model_sample_trial, type = "BG")$serial$p.value
  return(p_values)
}

# 5.) investigate residuals
plot_residuals <- function(res) {
  par(mfrow = c(3, 5))  
  for (i in 1:15) {  
    plot(res[, i], type = "p", xlab = "Time", ylab = "Residuals")
  }
  par(mfrow = c(1, 1))
}

GC_analysis <- function(LFP, lag_order) {
  
  num_trials <- length(LFP)
  num_channels <- unique(lengths(LFP))
  
  # prepare matrix for results
  combinations <- matrix(nrow = 15 * 15 - 15, ncol = 2)
  k <- 1
  for (i in 1:15) {
    for (j in 1:15) {
      # exclude identical pairs
      if (i != j) {
        combinations[k, 1] <- i
        combinations[k, 2] <- j
        k <- k + 1
      }
    }
  }
  combinations <- combinations[, c(2, 1)]
  p_values_trials <- matrix(nrow = 15 * 15 - 15, ncol = num_trials)
  GC_p_values <- cbind(combinations, p_values_trials)
  colnames(GC_p_values) <- c("influencing_ch", "influenced_ch", paste0("trial_", 1:num_trials, "_p_value"))
  
  adj_sum_square_err <- matrix(nrow = num_trials, ncol = num_channels)
  consistency <- rep(NA, num_trials)
  # dw_VAR_test_res <- matrix(nrow = num_trials, ncol = num_channels)
  portmanteu_test_res <- rep(NA, num_trials)
  bg_test_res <- rep(NA, num_trials)
  res <- list()
  
  for (i in 1:num_trials) {
    
    # extract trial i and transform into a matrix
    LFP_trial_i <- LFP[[i]]
    LFP_trial_i <- do.call(cbind, lapply(LFP_trial_i, function(x) as.numeric(x)))
    
    # VAR model creation
    VAR_model_trial_i <- VAR(LFP_trial_i, type = "const", p = lag_order)
    
    # checks
    adj_sum_square_err[i, ] <- calculate_adj_sum_square_err(VAR_model_trial_i)
    consistency[i] <- check_consistency(VAR_model_trial_i)
    # dw_VAR_test_res[i, ] <- dW_VAR_test(VAR_model_trial_i)
    portmanteu_test_res[i] <- vars::serial.test(VAR_model_trial_i, type = "PT.adjusted")$serial$p.value
    bg_test_res[i] <- vars::serial.test(VAR_model_trial_i, type = "BG")$serial$p.value
    res[[i]] <- residuals(VAR_model_trial_i)
    
    # Granger causality
    GC <- granger_causality(VAR_model_trial_i)
    GC_trial_i <- GC$result
    
    # extract p-values from F-Test
    p_values_F <- round(GC_trial_i$p.F[-(seq(from = 15, to = 225, by = 15))], 6)
    
    # Bonferroni adjust the p values
    adj_p_values_F <- p.adjust(p_values_F, method = "bonferroni")
    
    GC_p_values[, i + 2] <- adj_p_values_F
    
    print(paste("GC calculation of trial", i, "of", num_trials, "succeeded."))
    
  }
  
  return(list(
    p_values = GC_p_values, 
    adj_sum_square_err = adj_sum_square_err, 
    consistency = consistency, 
    # dw_VAR_test_res = dw_VAR_test_res,
    portmanteu_test_res = portmanteu_test_res,
    bg_test_res = bg_test_res,
    res = res
  ))
  
}

calculate_GC_p_values <- function(LFP, un_primed_ind, lag_order) {
  
  # filter either primed or unprimed trials
  LFP_trials <- LFP[un_primed_ind]
  
  # GC analysis
  GC_p_values <- GC_analysis(LFP = LFP_trials, lag_order = lag_order)
  
  return(GC_p_values)
  
}



# Network Visualisation of GC Influences (Layer Levels) -------------------

transform_layer_levels <- function(GC_p_values) {
  
  # prepare matrix with layer information and p values
  layer_influencing <- rep(NA, nrow(GC_p_values$p_values))
  layer_influencing[GC_p_values$p_values[, 1] %in% upper_channels] <- "upper"
  layer_influencing[GC_p_values$p_values[, 1] %in% middle_channels] <- "middle"
  layer_influencing[GC_p_values$p_values[, 1] %in% deep_channels] <- "deep"
  layer_influenced <- rep(NA, nrow(GC_p_values$p_values))
  layer_influenced[GC_p_values$p_values[, 2] %in% upper_channels] <- "upper"
  layer_influenced[GC_p_values$p_values[, 2] %in% middle_channels] <- "middle"
  layer_influenced[GC_p_values$p_values[, 2] %in% deep_channels] <- "deep"
  GC_p_values_layer_levels <- cbind(layer_influencing, layer_influenced)
  GC_p_values_layer_levels <- as.data.frame(GC_p_values_layer_levels)
  
  # calculate the percentage of significant GC combis
  percentage_sign <- cbind(
    c("upper", "upper", "upper", "middle", "middle", "middle", "deep", "deep", "deep"),
    c("upper", "middle", "deep", "upper", "middle", "deep", "upper", "middle", "deep")
  )
  percentage_sign_layer_level <- rep(NA, nrow(percentage_sign))
  for (i in 1:nrow(percentage_sign)) {
    ind <- GC_p_values_layer_levels[, 1] %in% percentage_sign[i, 1] & GC_p_values_layer_levels[, 2] %in% percentage_sign[i, 2]
    p_values_layer_level <- GC_p_values$p_values[ind, 3:ncol(GC_p_values$p_values)]
    percentage_sign_layer_level[i] <- sum(p_values_layer_level < 0.05) / (nrow(p_values_layer_level) * ncol(p_values_layer_level))
  }
  
  percentage_sign_layer_level <- cbind(percentage_sign, percentage_sign_layer_level)
  
  return(percentage_sign_layer_level)
  
}

# plot the influences (unprimed, primed separate)
plot_layer_influences <- function(percentage_sign, un_primed) {
  
  # assign colours to vertices based on influence
  min_influence <- min(as.numeric(percentage_sign[, 3]))
  max_influence <- max(as.numeric(percentage_sign[, 3]))
  colour_range <- colorRampPalette(c("lightblue", "darkblue"))(10)
  influence_colours <- colour_range[cut(
    as.numeric(percentage_sign[, 3]),
    breaks = seq(min_influence, max_influence, length.out = length(colour_range) - 1),
    include.lowest = TRUE
  )]
  
  # plot the network with coloured nodes based on influence
  plot(
    graph_from_data_frame(percentage_sign[, 1:2], directed = TRUE),
    layout = layout.circle,
    vertex.label.cex = 1.2,
    vertex.size = 50,
    edge.arrow.size = 1.2,
    edge.width = 3,
    edge.color = influence_colours,
    edge.curved = rep(0.1, nrow(percentage_sign)),
    vertex.color = "lightblue",
    vertex.frame.color = "lightblue"
  )
  legend(
    "bottomright",
    legend = c("Smaller Influence", "Bigger Influence"),
    col = c(colour_range[1], colour_range[length(colour_range)]),
    lty = 1,
    lwd = 2,
    cex = 0.8
  )
  
}
