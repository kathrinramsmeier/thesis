
library(tidyverse)  # data manipulation
library(vars)       # VAR model
library(bruceR)     # Granger causality
library(multitaper) # Multitaper spectral analysis
library(igraph)     # graphs for GC
library(ggplot2)    # plots

options(scipen = 1000)



# Load Session Data -------------------------------------------------------

# delete data from previous session to get storage space
rm(list = c("LFP", "probe", "time", "behaviour", "task", "recordinginfo", "probe"))

load_data <- function(session) {
  
  # detach("package:reticulate")
  library(reticulate)
  np <- import("numpy")
  
  setwd("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\")
  
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

data_preprocessing <- function(session_data, same_lengths) {
  
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
  
  # for clipping: get the ind in the time vector that is closest to the reaction time
  reaction_time_ind <- rep(NA, dim(session_data$LFP)[3])
  for (i in 1:dim(session_data$LFP)[3]) {
    reaction_time_ind[i] <- which.min(abs(session_data$time - session_data$behaviour$reaction_time_ms[i]))
  }
  
  if (same_lengths) {
    
    # exclude about 10 ms before the reaction
    clip <- reaction_time_ind - 10
    time_ind <- 1:min(clip)
    time_i <- session_data$time[time_ind]
    
    # only keep LFP after stimulus onset
    time_ind <- which(time_i > 0)
    
    # only keep the LFP after stimulus onset up to the minimum of 10 ms before reaction times
    time <- matrix(nrow = dim(session_data$LFP)[3], ncol = length(time_ind))
    LFP <- array(dim = c(dim(session_data$LFP)[1], length(time_ind), dim(session_data$LFP)[3]))
    for (i in 1:dim(session_data$LFP)[3]) {
      time[i, ] <- time_i[time_ind]
      for (j in 1:dim(LFP)[1]) {
        LFP[j, , i] <- session_data$LFP[j, time_ind, i]
      }
    }
    
  } else {
    
    # exclude about 10 ms before the reaction
    clip <- reaction_time_ind - 10
    time <- list()
    LFP <- vector("list", length = dim(session_data$LFP)[3])
    
    for (i in 1:dim(session_data$LFP)[3]) {
      
      # exclude about 10 ms before the reaction
      time_ind <- 1:clip[i]
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



# Necessary Variables -----------------------------------------------------

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



# Check for Stationarity --------------------------------------------------

# Augmented Dickey-Fuller test (H0: non-stationary) either for all channels or a specific one
adf_stationarity_test <- function(LFP, channel) {
  
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
  
  return(percentage_stationary)
  
}



# Dealing With Non-Stationarity - Notch Filtering -------------------------

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
      adf_result_notch_channel_k[j] <- adf_stationarity_test(LFP_filtered, channel = "all")
    }
    # get the filter order for which the most trials are stationary
    adf_result_notch[k] <- which.max(adf_result_notch_channel_k)
  }
  
  return(adf_result_notch)
  
}



# Stationarity Transformation on the Whole LFP ----------------------------

transform_stationary <- function(
    LFP,
    detrending,
    normalising,
    ensemble_adjusting,
    notch_filtering,
    notch_filter_order,
    sampling_rate,
    print_remaining_non_stationary,
    hard_remove_all
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
        n = notch_filter_order[i], 
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
        p_values[j] <- tseries::adf.test(LFP_channel_i[, j])$p.value
      }
      print(paste0("Channel ", i, ": Number of LFP trials still non-stationary (before removing all): ", sum(p_values > 0.05)))
    }
    
    LFP_stationary <- abind::abind(LFP_stationary, LFP_channel_i, along = 3)
    
  }
  
  LFP_stationary <- LFP_stationary[, , 2:(dim(LFP)[1] + 1)]
  LFP_stationary <- aperm(LFP_stationary, c(3, 1, 2))
  
  # exclude the remaining few non-stationary trials
  if (hard_remove_all) {
    
    p_values <- matrix(NA, nrow = dim(LFP_stationary)[3], ncol = dim(LFP_stationary)[1])
    
    # loop over all channels
    for (j in 1:dim(LFP_stationary)[1]) {
      
      # loop over all trials
      for (i in 1:dim(LFP_stationary)[3]) {
        p_values[i, j] <- tseries::adf.test(LFP_stationary[j, , i])$p.value
      }
      
    }
    
    all_channels_stationary <- rep(NA, nrow(p_values))
    for (j in 1:nrow(p_values)) {
      
      # check if trials in all channels are stationary
      all_channels_stationary[j] <- all(p_values[j, ] < 0.05)
      stationary_ind <- which(all_channels_stationary)
      
    }
    
    # filter the stationary trials
    LFP_stationary <- LFP_stationary[, , stationary_ind]
    
  }
  
  return(LFP_stationary)
  
}

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
    freq_resolution <- sampling_rate / length(LFP_fftransformed)
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



# Filter the LFP in the different frequency bands -------------------------

freqband_filtering <- function(LFP_trial, frequency_bands, sampling_rate, filter_order) {
  
  LFP_freqband_filtered <- matrix(nrow = length(LFP_trial), ncol = length(frequency_bands))
  colnames(LFP_freqband_filtered) <- names(frequency_bands)
  
  nyquist_freq <- sampling_rate / 2
  
  # avoid numerical problem
  frequency_bands[[1]][1] <- 0.01
  
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

# Coherence ---------------------------------------------------------------

# calculate coherence between two signals
calculate_coherence <- function(LFP_trial_1_freqband_filtered, LFP_trial_2_freqband_filtered, meaned) {
  
  coherences <- matrix(nrow = 9, ncol = ncol(LFP_trial_1_freqband_filtered))
  
  # loop over all frequency bands
  for (i in 1:ncol(LFP_trial_1_freqband_filtered)) {
    
    # cross-power spectral magnitudes and density autospectral densities
    # formulas from:
    # https://en.wikipedia.org/wiki/Coherence_(signal_processing)
    # https://www.psychologie.hhu.de/fileadmin/redaktion/Fakultaeten/Mathematisch-Naturwissenschaftliche_Fakultaet/Psychologie/CompPsy/Papers/vinck2010_wingerden.pdf
    csd <- gsignal::cpsd(cbind(LFP_trial_1_freqband_filtered[, i], LFP_trial_2_freqband_filtered[, i]))$cross
    csd_magnitudes <- abs(csd)
    asd_1 <- gsignal::cpsd(cbind(LFP_trial_1_freqband_filtered[, i], LFP_trial_1_freqband_filtered[, i]))$cross
    asd_2 <- gsignal::cpsd(cbind(LFP_trial_2_freqband_filtered[, i], LFP_trial_2_freqband_filtered[, i]))$cross
    
    # coherence
    coherences[, i] <- csd_magnitudes^2 / (asd_1 * asd_2)
    
  }
  
  # mean coherence in each frequency band !!!! need to check if that is plausible to do !!!!
  if(meaned) {
    coherences <- colMeans(coherences)
    names(coherences) <- names(frequency_bands)
    
  } else {
    colnames(coherences) <- names(frequency_bands)
  }
  
  return(coherences)
  
}

# compute coherence between the channels in each frequency band for all primed and unprimed trials
calculate_avg_coherence <- function(LFP, un_primed_ind, frequency_bands, sampling_rate, filter_order) {
  
  # filter either primed or unprimed trials
  LFP_trials <- LFP[, , un_primed_ind]
  
  # prepare matrices to store results
  coherences <- array(dim = c(dim(LFP)[1], dim(LFP)[1], length(frequency_bands)))
  sum_coherences <- array(0 ,dim = c(dim(LFP)[1], dim(LFP)[1], length(frequency_bands)))
  
  # loop over all trials
  for (k in 1:dim(LFP_trials)[3]) {
    
    # coherence for trial i between each channel
    for (i in 1:dim(LFP)[1]) {
      for (j in 1:dim(LFP)[1]) {
        trial_freqband_filtered_1 <- freqband_filtering(
          LFP_trial = LFP_trials[i, , k], 
          frequency_bands = frequency_bands, 
          sampling_rate = sampling_rate,
          filter_order = filter_order
        )
        trial_freqband_filtered_2 <- freqband_filtering(
          LFP_trial = LFP_trials[j, , k], 
          frequency_bands = frequency_bands, 
          sampling_rate = sampling_rate,
          filter_order = filter_order
        )
        coherences[i, j, ] <- calculate_coherence(
          LFP_trial_1_freqband_filtered = trial_freqband_filtered_1, 
          LFP_trial_2_freqband_filtered = trial_freqband_filtered_2,
          meaned = TRUE
        )
      }
    }
    
    # calculate mean coherence over all trials (need to think again if that is plausible to do!!!)
    sum_coherences <- sum_coherences + coherences
    avg_coherences <- sum_coherences / dim(LFP_trials)[3]
    
  }
  
  return(avg_coherences)
  
}

# plot the average coherence between the channels for the different frequency bands as heatmaps
plot_heatmap <- function(result_array, frequency_bands) {
  
  # preparation of the data
  coherence_df <- as.data.frame(as.table(result_array))
  colnames(coherence_df) <- c("channel1", "channel2", "frequency_band", "Coherence")
  coherence_df$channel1 <- as.numeric(as.factor(coherence_df$channel1))
  coherence_df$channel2 <- as.numeric(as.factor(coherence_df$channel2))
  coherence_df$frequency_band <- factor(coherence_df$frequency_band)
  levels(coherence_df$frequency_band) <- names(frequency_bands)
  
  # heatmap plot
  ggplot(coherence_df, aes(x = channel1, y = channel2, fill = Coherence)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") + 
    facet_wrap( ~ frequency_band, scales = "free") +  
    labs(title = "Avg Coherence Between the Different Channels as a Function of Frequency Bands", x = "", y = "", fill = "Coherence") +
    theme_minimal()
  
}



# Phase Locking Value (PLV) and Phase Lag Index (PLI) ---------------------

# compute PLV and PLI across the frequency bands
# formulas from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3674231/#:~:text=PLV%20can%20therefore%20be%20viewed,each%20scaling%20of%20the%20wavelet.
calculate_PLV_PLI <- function(LFP_trial_1_freqband_filtered, LFP_trial_2_freqband_filtered, PLV_PLI) {
  
  stopifnot(PLV_PLI %in% c("PLV", "PLI"))
  
  result <- rep(NA, dim(LFP_trial_1_freqband_filtered)[2])
  names(result) <- colnames(LFP_trial_1_freqband_filtered)
  
  # loop over all frequency bands
  for (i in 1:dim(LFP_trial_1_freqband_filtered)[2]) {
    
    # Hilbert transform to get the analytical signal
    LFP_1_freqband_ht <- gsignal::hilbert(LFP_trial_1_freqband_filtered[, i])
    LFP_2_freqband_ht <- gsignal::hilbert(LFP_trial_2_freqband_filtered[, i])
    
    # compute phase difference
    phase_diffs <- pracma::angle(LFP_1_freqband_ht * Conj(LFP_2_freqband_ht))  # arg((LFP_1_band_ht * Conj(LFP_2_band_ht)) / (abs(LFP_1_band_ht) * abs(LFP_2_band_ht)))
    
    if (PLV_PLI == "PLV") {
      # compute the PLV for each frequency band
      result[i] <- abs(mean(exp(1i * phase_diffs)))
    } else {
      # compute the PLI for each frequency band
      result[i] <- abs(mean(sign(phase_diffs)))
    }
  }
  
  return(result)
  
}

# compute PPC: PLV/PLI between the channels in each frequency band for all primed and unprimed trials
calculate_PPC <- function(LFP, un_primed_ind, frequency_bands, sampling_rate, filter_order, PLV_PLI) {
  
  # filter either primed or unprimed trials
  LFP_trials <- LFP[, , un_primed_ind]
  
  # prepare matrices to store results
  results <- array(dim = c(dim(LFP)[1], dim(LFP)[1], length(frequency_bands)))
  sum_results <- array(0 ,dim = c(dim(LFP)[1], dim(LFP)[1], length(frequency_bands)))
  
  # loop over all trials
  for (k in 1:dim(LFP_trials)[3]) {
    
    # PLV/PLI for trial i between each channel
    for (i in 1:dim(LFP)[1]) {
      for (j in 1:dim(LFP)[1]) {
        trial_freqband_filtered_1 <- freqband_filtering(
          LFP_trial = LFP_trials[i, , k], 
          frequency_bands = frequency_bands, 
          sampling_rate = sampling_rate,
          filter_order = filter_order
        )
        trial_freqband_filtered_2 <- freqband_filtering(
          LFP_trial = LFP_trials[j, , k], 
          frequency_bands = frequency_bands, 
          sampling_rate = sampling_rate,
          filter_order = filter_order
        )
        results[i, j, ] <- calculate_PLV_PLI(
          LFP_trial_1_freqband_filtered = trial_freqband_filtered_1, 
          LFP_trial_2_freqband_filtered = trial_freqband_filtered_2,
          PLV_PLI = PLV_PLI
        )
      }
    }
    
    # calculate mean PLV/PLI over all trials 
    # see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3674231/#:~:text=PLV%20can%20therefore%20be%20viewed,each%20scaling%20of%20the%20wavelet.
    sum_results <- sum_results + results
    avg_results <- sum_results / dim(LFP_trials)[3]
    
  }
  
  return(avg_results)
  
}



# GC Analysis for All Trials ----------------------------------------------

GC_analysis <- function(LFP, lag_order) {
  
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
  p_values_trials <- matrix(nrow = 15 * 15 - 15, ncol = dim(LFP)[3])
  GC_p_values <- cbind(combinations, p_values_trials)
  colnames(GC_p_values) <- c("influencing_ch", "influenced_ch", paste0("trial_", 1:dim(LFP)[3], "_p_value"))
  
  # loop over all trials
  for (i in 1: dim(LFP)[3]) {
    
    LFP_trial_i <- LFP[, , i]
    LFP_trial_i <- t(LFP_trial_i) 
    
    # VAR model creation (using AIC to determine the lag order)
    VAR_model_trial_i <- VAR(LFP_trial_i, type = "const", p = lag_order)
    
    # Granger causality
    GC <- granger_causality(VAR_model_trial_i)
    GC_trial_i <- GC$result
    
    # extract p-values from F-Test
    p_values_F <- round(GC_trial_i$p.F[-(seq(from = 15, to = 225, by = 15))], 6)
    
    # Bonferroni adjust the p values
    adj_p_values_F <- p.adjust(p_values_F, method = "bonferroni")
    
    GC_p_values[, i + 2] <- adj_p_values_F
    
    print(paste("GC calculation of trial", i, "of", dim(LFP)[3], "succeeded."))
    
  }
  
  return(GC_p_values)
  
}

calculate_GC_p_values <- function(LFP, un_primed_ind, lag_order) {
  
  # filter either primed or unprimed trials
  LFP_trials <- LFP[, , un_primed_ind]
  
  # GC analysis
  GC_p_values <- GC_analysis(LFP = LFP_trials, lag_order = lag_order)
  
  return(GC_p_values)
  
}



# Network Visualisation of GC Influences (Layer Levels) -------------------

transform_layer_levels <- function(GC_p_values) {
  
  # prepare matrix with layer information and p values
  layer_influencing <- rep(NA, nrow(GC_p_values))
  layer_influencing[GC_p_values[, 1] %in% upper_channels] <- "upper"
  layer_influencing[GC_p_values[, 1] %in% middle_channels] <- "middle"
  layer_influencing[GC_p_values[, 1] %in% deep_channels] <- "deep"
  layer_influenced <- rep(NA, nrow(GC_p_values))
  layer_influenced[GC_p_values[, 2] %in% upper_channels] <- "upper"
  layer_influenced[GC_p_values[, 2] %in% middle_channels] <- "middle"
  layer_influenced[GC_p_values[, 2] %in% deep_channels] <- "deep"
  GC_p_values_layer_levels <- cbind(layer_influencing, layer_influenced, GC_p_values)
  GC_p_values_layer_levels <- as.data.frame(GC_causality_p_values_layer_levels)
  
  # calculate the percentage of significant GC combis
  percentage_sign <- cbind(
    c("upper", "upper", "upper", "middle", "middle", "middle", "deep", "deep", "deep"),
    c("upper", "middle", "deep", "upper", "middle", "deep", "upper", "middle", "deep")
  )
  percentage_sign_layer_level <- rep(NA, nrow(percentage_sign))
  for (i in 1:nrow(percentage_sign)) {
    ind <- GC_p_values_layer_levels[, 1] %in% percentage_sign[i, 1] & GC_p_values_layer_levels[, 2] %in% percentage_sign[i, 2]
    p_values_layer_level <- GC_p_values_layer_levels[ind, ]$adj_p_values_F
    percentage_sign_layer_level[i] <- sum(p_values_layer_level < 0.05) / length(p_values_layer_level)
  }
  
  return(percentage_sign_layer_level)
  
}

# plot the influences (unprimed, primed separate)
plot_layer_influences <- function(percentage_sign, un_primed) {
  
  # assign colours to vertices based on influence
  min_influence <- min(as.numeric(percentage_sign[, 3]))
  max_influence <- max(as.numeric(percentage_sign[, 3]))
  colour_range <- cm.colors(10)
  influence_colours <- colour_range[cut(
    as.numeric(percentage_sign[, 3]),
    breaks = seq(min_influence, max_influence, length.out = length(colour_range) - 1),
    include.lowest = TRUE
  )]
  
  # plot the network with coloured nodes based on influence
  plot(
    graph_from_data_frame(percentage_sign[, 1:2], directed = TRUE),
    layout = layout.circle,
    vertex.label.cex = 1.5,
    vertex.size = 30,
    edge.arrow.size = 0.5,
    edge.width = 3,
    edge.color = influence_colours,
    edge.curved = rep(0.1, nrow(percentage_sign)),
    vertex.color = "white",
    main = paste0(
      "GC Relationships Between Cortical Layers for a Subset of ", un_primed, " Trials"
    )
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