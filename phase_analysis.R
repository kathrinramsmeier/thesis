
# Calculation of the phase coherence using Mean Resultant Length (MRL), Phase Locking Value (PLV) and Phase Lag Index (PLI).
# Coherence and PLV are equivalent in the Gaussian case (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3674231/#:~:text=PLV%20can%20therefore%20be%20viewed,each%20scaling%20of%20the%20wavelet.).
# but here: individual LFP trials are not normally distributed



# Mean Resultant Length (MRL) ---------------------------------------------

# look at the 6th and 7th channel of trial 333 and 6th channel of trial 999 as an example first
trial <- 333
channel <- 6
sample_LFP_1 <- LFP_stationary[channel, , trial]
channel <- 7
sample_LFP_2 <- LFP_stationary[channel, , trial]
trial <- 999
sample_LFP_3 <- LFP_stationary[channel, , trial]
plot(sample_LFP_2, type = "l")
lines(sample_LFP_3, type = "l", col = "red")

# calculate mean resultant length for each frequency band
calculate_MRL <- function(LFP_1, LFP_2, frequency_bands, freq_axis) {
  
  # FFT to obtain phase information
  LFP_fft_1 <- fft(LFP_1)
  LFP_fft_2 <- fft(LFP_2)
  
  # Hilbert transform to extract instantaneous phase
  phase_1 <- atan2(Im(LFP_fft_1), Re(LFP_fft_1)) # Re = magnitude, Im = phase 
  phase_2 <- atan2(Im(LFP_fft_2), Re(LFP_fft_2))
  
  # MRL
  MRL <- rep(NA, length(frequency_bands))
  for (i in 1:length(frequency_bands)) {
    phase_x_band <- phase_1[freq_axis >= frequency_bands[[i]][1] & freq_axis < frequency_bands[[i]][2]]
    phase_y_band <- phase_2[freq_axis >= frequency_bands[[i]][1] & freq_axis < frequency_bands[[i]][2]]
    MRL[[i]] <- abs(mean(exp(1i * (phase_x_band - phase_y_band))))
  }
  
  return(MRL)
}

# try MRL calc. for a sample trial for two different channels
MRL_ch1_ch2 <- calculate_MRL(
  LFP_1 = sample_LFP_1, 
  LFP_2 = sample_LFP_2,
  frequency_bands = frequency_bands,
  freq_axis = freq_axis
)
MRL_ch1_ch2

# try MRL calc. for two different sample trials
MRL_tr1_tr2 <- calculate_MRL(
  LFP_1 = sample_LFP_2, 
  LFP_2 = sample_LFP_3,
  frequency_bands = frequency_bands,
  freq_axis = freq_axis
)
MRL_tr1_tr2

# compute pairwise MRL between the trials in each frequency band in channel 1
MRL_ch_1 <- array(dim = c(dim(LFP_stationary)[3], dim(LFP_stationary)[3], length(frequency_bands)))
for (i in 1:dim(LFP_stationary)[3]) {
  for (j in 1:dim(LFP_stationary)[3]) {
    MRL_ch_1[i, j, ] <- calculate_MRL(
      LFP_1 = LFP_stationary[1, , i], 
      LFP_2 = LFP_stationary[1, , j],
      freq = freq_axis,
      frequency_bands = frequency_bands
    )
  }
}
dim(MRL_ch_1) # dim 1 and 2: trial x trial combi, dim 3: frequency bands
head(MRL_ch_1)

# averaging MRL over trials for channel 1 (TO DO: research and think if this is a plausible step!)
MRL_hat <- apply(MRL_ch_1, 2, mean)
MRL_hat

# compute MRL_hat between the trials in each frequency band in channel i
MRL_hat_calculation <- function(LFP, un_primed_ind, freq, frequency_bands) {
  
  # filter either primed or unprimed trials
  LFP_trials <- LFP[, , un_primed_ind]
  
  # prepare matrix to store results
  MRL_hat_results <- matrix(nrow = dim(LFP_trials)[1], ncol = length(frequency_bands))
  rownames(MRL_hat_results) <- 1:dim(LFP_trials)[1]
  colnames(MRL_hat_results) <- names(frequency_bands)
  
  # compute pairwise MRL between the trials in each frequency band in each channel
  for (i in 1:dim(LFP_trials)[1]) {
    MRL_channel_i <- matrix(ncol = length(frequency_bands), nrow = dim(LFP_trials)[3])
    for (j in 1:dim(LFP_trials)[3]) {
      for (k in 1:dim(LFP_trials)[3]) {
        MRL_channel_i[j, ] <- calculate_MRL(
          LFP_1 = LFP_trials[i, , j], 
          LFP_2 = LFP_trials[i, , k],
          freq = freq,
          frequency_bands = frequency_bands
        )
      }
    }
    
    # average across trials
    MRL_hat_results[i, ] <- apply(MRL_channel_i, 2, mean)
  }
  
  return(MRL_hat_results)
  
}

# # test for a subset
# temp <- LFP_stationary[, , 1:10]
# temp2 <- unprimed_ind[unprimed_ind <= 10]
# MRL_hat_test <- MRL_hat_calculation(
#   LFP = temp,
#   un_primed_ind = temp2,
#   freq = freq_axis,
#   frequency_bands = frequency_bands
# )
# MRL_hat_test

# calculate MRL_hat for all frequency bands for unprimed trials
MRL_hat_unprimed <- MRL_hat_calculation(
  LFP = LFP_stationary,
  un_primed_ind = unprimed_ind,
  freq = freq_axis,
  frequency_bands = frequency_bands
)

# calculate PLI_hat for all frequency bands for primed trials
MRL_hat_primed <- MRL_hat_calculation(
  LFP = LFP_stationary,
  un_primed_ind = primed_ind,
  freq = freq_axis,
  frequency_bands = frequency_bands
)

# compute pairwise MRL between the channels in each frequency band in a sample trial
trial <- 333
MRL_sample_trial <- array(dim = c(dim(LFP_stationary)[1], dim(LFP_stationary)[1], length(frequency_bands)))
for (i in 1:dim(LFP_stationary)[1]) {
  for (j in 1:dim(LFP_stationary)[1]) {
    MRL_sample_trial[i, j, ] <- calculate_MRL(
      LFP_1 = LFP_stationary[i, , trial], 
      LFP_2 = LFP_stationary[j, , trial],
      freq = freq_axis,
      frequency_bands = frequency_bands
    )
  }
}
dim(MRL_sample_trial) # dim 1 and 2: channel x channel combi, dim 3: frequency bands
MRL_sample_trial

# plot phase coherence of channel i with all other channels over the different frequency bands for one specific trial
plot_phase_coherence <- function(pc_trial_matrix) {
  
  plot_list <- list()
  for (i in 1:15) {
    
    pc_ch_i <- as.data.frame(pc_trial_matrix[i, , ])
    colnames(pc_ch_i) <- names(frequency_bands)
    rownames(pc_ch_i) <- 1:dim(LFP)[1]
    pc_ch_i_long <- tidyr::gather(as.data.frame(pc_ch_i), key = "Frequency_band", value = "Phase_coherence", factor_key = TRUE)
    
    plot_i <- ggplot(pc_ch_i_long, aes(x = as.numeric(rownames(pc_ch_i_long)), y = Phase_coherence, color = Frequency_band)) +
      geom_point(size = 2) +
      labs(x = "", y = "Phase Coherence", title = paste("Phase Coherence of Channel", i, "and All Other Channels for a Sample Trial")) +
      scale_color_discrete(name = "Frequency Band") +
      theme_minimal() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 10)
      )
    
    plot_list[[i]] <- plot_i
    
  }
  grid.arrange(grobs = plot_list, ncol = 3) 
  
}
plot_phase_coherence(pc_trial_matrix = MRL_sample_trial)

# # calculate phase coherence for all trials
# calculate_phase_coherence <- function(LFP, un_primed_ind, freq, frequency_bands) {
#   
#   # filter either primed or unprimed trials and the electrode channel
#   LFP_trials <- LFP[, , un_primed_ind]
#   
#   # prepare matrix for results
#   pc_trials <- array(dim = c(dim(LFP_trials)[3], dim(LFP)[1], dim(LFP)[1], length(frequency_bands)))
#   
#   for (i in 1:dim(pc_trials)[1]) {
#     
#     # compute pairwise MRL between the channels in each frequency band in a trial i
#     MRL_trial_i <- array(dim = c(dim(LFP)[1], dim(LFP)[1], length(frequency_bands)))
#     for (j in 1:dim(LFP)[1]) {
#       for (k in 1:dim(LFP)[1]) {
#         MRL_trial_i[j, k, ] <- calculate_MRL(
#           LFP_1 = LFP[j, , i], 
#           LFP_2 = LFP[k, , i],
#           freq = freq_axis,
#           frequency_bands = frequency_bands
#         )
#       }
#     }
#     
#     pc_trials[i, , , ] <- MRL_trial_i
#     
#   }
#   
#   return(pc_trials)
#   
# }
# 
# # test for a subset
# temp <- LFP_stationary[, , 1:10]
# temp2 <- unprimed_ind[unprimed_ind <= 10]
# phase_coherence_test <- calculate_phase_coherence(
#   LFP = temp,
#   un_primed_ind = temp2,
#   freq = freq_axis,
#   frequency_bands = frequency_bands
# )
# dim(phase_coherence_test)
# 
# # phase coherence of channel i with all other channels over the different frequency bands for all unprimed trials
# phase_coherence_unprimed <- calculate_phase_coherence(
#   LFP = LFP_stationary, 
#   un_primed_ind = unprimed_ind,
#   freq = freq_axis,
#   frequency_bands = frequency_bands
# )
# dim(phase_coherence_unprimed) # unprimed trials, channels, channels, frequency bands
# 
# # phase coherence of channel i with all other channels over the different frequency bands for all primed trials
# phase_coherence_primed <- calculate_phase_coherence(
#   LFP = LFP_stationary, 
#   un_primed_ind = primed_ind,
#   freq = freq_axis,
#   frequency_bands = frequency_bands
# )
# dim(phase_coherence_primed) # primed trials, channels, channels, frequency bands
# 
# # plot phase coherence of channel i with all other channels over the different frequency bands for a sample unprimed trial
# set.seed(42)
# trial_ind <- sample(1:length(unprimed_ind), 1)
# plot_phase_coherence(pc_trial_matrix = phase_coherence_unprimed[trial_ind, , , ])
# 
# # plot phase coherence of channel i with all other channels over the different frequency bands for a sample primed trial
# set.seed(42)
# trial_ind <- sample(1:length(primed_ind), 1)
# plot_phase_coherence(pc_trial_matrix = phase_coherence_primed[trial_ind, , , ])



# Phase Locking Value (PLV) and Phase Lag Index (PLI) ---------------------

# compute PLV and PLI across the frequency bands
# formulas from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3674231/#:~:text=PLV%20can%20therefore%20be%20viewed,each%20scaling%20of%20the%20wavelet.
calculate_PLV_PLI <- function(LFP_1, LFP_2, PLV_PLI, freq, frequency_bands) {
  
  stopifnot(PLV_PLI %in% c("PLV", "PLI"))
  
  result <- rep(NA, length(frequency_bands))
  names(result) <- names(frequency_bands)
  
  for (i in 1:length(frequency_bands)) {
    
    # filter the LFPs within the frequency band
    band_ind <- which(freq >= frequency_bands[[i]][1] & freq < frequency_bands[[i]][2])
    LFP_1_band <- LFP_1[band_ind]
    LFP_2_band <- LFP_2[band_ind]
    
    # Hilbert transform to get the analytical signal
    LFP_1_band_ht <- gsignal::hilbert(LFP_1_band)
    LFP_2_band_ht <- gsignal::hilbert(LFP_2_band)
    
    # relative phase
    rel_phase <- pracma::angle(LFP_1_band_ht * Conj(LFP_2_band_ht))  # arg((LFP_1_band_ht * Conj(LFP_2_band_ht)) / (abs(LFP_1_band_ht) * abs(LFP_2_band_ht)))
    
    if(PLV_PLI == "PLV") {
      # compute the PLV for each frequency band
      result[i] <- abs(mean(exp(1i * rel_phase)))
    } else {
      # compute the PLI for each frequency band
      result[i] <- abs(mean(sign(rel_phase)))
    }
    
  }
  
  return(result)
  
}

# try PLV calc. for a sample trial for two different channels
PLV_t1_t2 <- calculate_PLV_PLI(
  LFP_1 = sample_LFP_1, 
  LFP_2 = sample_LFP_2,
  PLV_PLI = "PLV",
  freq = freq_axis,
  frequency_bands = frequency_bands
)
PLV_t1_t2

# try PLI calc. for a sample trial for two different channels
PLI_t1_t2 <- calculate_PLV_PLI(
  LFP_1 = sample_LFP_1, 
  LFP_2 = sample_LFP_2,
  PLV_PLI = "PLI",
  freq = freq_axis,
  frequency_bands = frequency_bands
)
PLI_t1_t2

# try PLV calc. for two different sample trials
PLV_t2_t3 <- calculate_PLV_PLI(
  LFP_1 = sample_LFP_2, 
  LFP_2 = sample_LFP_3,
  PLV_PLI = "PLV",
  freq = freq_axis,
  frequency_bands = frequency_bands
)
PLV_t2_t3

# try PLI calc. for two different sample trials
PLI_t2_t3 <- calculate_PLV_PLI(
  LFP_1 = sample_LFP_2, 
  LFP_2 = sample_LFP_3,
  PLV_PLI = "PLI",
  freq = freq_axis,
  frequency_bands = frequency_bands
)
PLI_t2_t3

# compute pairwise PLV between the trials in each frequency band in channel 1
PLV_channel_1 <- matrix(ncol = length(frequency_bands), nrow = dim(LFP_stationary)[3])
colnames(PLV_channel_1) <- names(frequency_bands)
for (i in 1:dim(LFP_stationary)[3]) {
  for (j in 1:dim(LFP_stationary)[3]) {
    PLV_channel_1[i, ] <- calculate_PLV_PLI(
      LFP_1 = LFP_stationary[1, , i], 
      LFP_2 = LFP_stationary[1, , j],
      PLV_PLI = "PLV",
      freq = freq_axis,
      frequency_bands = frequency_bands
    )
  }
}
dim(PLV_channel_1)
head(PLV_channel_1)

# averaging PLVs over trials for channel 1
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3674231/#:~:text=PLV%20can%20therefore%20be%20viewed,each%20scaling%20of%20the%20wavelet.

PLV_hat <- apply(PLV_channel_1, 2, mean)
PLV_hat

# compute PLV_hat or PLI_hat between the trials in each frequency band in channel i
calculate_PLV_PLI_hat <- function(LFP, PLV_PLI, un_primed_ind, freq, frequency_bands) {
  
  # filter either primed or unprimed trials
  LFP_trials <- LFP[, , un_primed_ind]
  
  # prepare matrix to store results
  PLV_PLI_hat_results <- matrix(nrow = dim(LFP_trials)[1], ncol = length(frequency_bands))
  rownames(PLV_PLI_hat_results) <- 1:dim(LFP_trials)[1]
  colnames(PLV_PLI_hat_results) <- names(frequency_bands)
  
  # compute pairwise PLV/PLI between the trials in each frequency band in each channel
  for (i in 1:dim(LFP_trials)[1]) {
    PLV_channel_i <- matrix(ncol = length(frequency_bands), nrow = dim(LFP_trials)[3])
    for (j in 1:dim(LFP_trials)[3]) {
      for (k in 1:dim(LFP_trials)[3]) {
        PLV_channel_i[j, ] <- calculate_PLV_PLI(
          LFP_1 = LFP_trials[i, , j], 
          LFP_2 = LFP_trials[i, , k],
          PLV_PLI = PLV_PLI,
          freq = freq,
          frequency_bands = frequency_bands
        )
      }
    }
    
    # average across trials
    PLV_PLI_hat_results[i, ] <- apply(PLV_channel_i, 2, mean)
  }
  
  return(PLV_PLI_hat_results)
  
}

# # test for a subset
# temp <- LFP_stationary[, , 1:10]
# temp2 <- unprimed_ind[unprimed_ind <= 10]
# PLV_hat_test <- calculate_PLV_PLI_hat(
#   LFP = temp,
#   PLV_PLI = "PLV",
#   un_primed_ind = temp2,
#   freq = freq_axis,
#   frequency_bands = frequency_bands
# )
# PLV_hat_test

# calculate PLI_hat for all frequency bands for unprimed trials
PLI_hat_unprimed <- calculate_PLV_PLI_hat(
  LFP = LFP_stationary, 
  PLV_PLI = "PLI",
  un_primed_ind = unprimed_ind,
  freq = freq_axis,
  frequency_bands = frequency_bands
)

# calculate PLI_hat for all frequency bands for primed trials
PLI_hat_primed <- calculate_PLV_PLI_hat(
  LFP = LFP_stationary,
  PLV_PLI = "PLI",
  un_primed_ind = primed_ind,
  freq = freq_axis,
  frequency_bands = frequency_bands
)






