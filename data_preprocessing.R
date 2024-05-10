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

session <- "C190127"
session_data <- load_data(session = session)



# Data Inspection ---------------------------------------------------------

session_data$LFP[1:5, 1:5, 1:5] # electrode channel (1:15), time, trial
dim(session_data$LFP)

summary(session_data$time) # before stimulus, 0 is stimulus, past-stimulus
length(session_data$time) # 2 values less than LFP
dim(session_data$behaviour)

session_data$probe # layer 1:5 = upper, 6:10 = middle, 11:15 = deep



# Data Preprocessing ------------------------------------------------------

data_preprocessing <- function(session_data) {
  
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
  
  # only keep the LFP after stimulus onset up to the minimum of 10 ms before reaction times
  time <- matrix(nrow = dim(session_data$LFP)[3], ncol = 131)
  LFP <- array(dim = c(dim(session_data$LFP)[1], 131, dim(session_data$LFP)[3]))
  for (i in 1:dim(session_data$LFP)[3]) {
    
    # exclude about 10 ms before the reaction
    clip_lower <- reaction_time_ind - 10
    time_ind <- 1:min(clip_lower)
    time_i <- session_data$time[time_ind]
    
    # only keep LFP after stimulus onset
    time_ind <- which(time_i > 0)
    time[i, ] <- time_i[time_ind]
    for (j in 1:dim(LFP)[1]) {
      LFP[j, , i] <- session_data$LFP[j, time_ind, i]
    }
  }

  session_data$time <- time
  session_data$LFP <- LFP
  rm(list = c("time", "LFP"))
  
  # checks
  stopifnot(
    dim(session_data$LFP)[1] == 15 & 
    dim(session_data$LFP)[2] == dim(session_data$time)[2] & 
    dim(session_data$LFP)[3] == nrow(session_data$task)
  )
  
  return(session_data)
  
}

session_data <- data_preprocessing(session_data = session_data)

# get the data sets out of the list into the global environment
list2env(session_data, envir = .GlobalEnv)
rm(list = "session_data")



# Necessary Variables -----------------------------------------------------

# filter all primed and unprimed trials (indices)
primed_ind <- task$trial_number_count[task$block_trial_count %in% c(1, 2)]
unprimed_ind <- task$trial_number_count[!(task$block_trial_count %in% c(1, 2))]

# sampling rate
sampling_rate <- recordinginfo$data_sampling_rate_hz

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

# get channel information
upper_channels <- probe$contact[probe$layer_category == "upper"]
middle_channels <- probe$contact[probe$layer_category == "middle"]
middle_channels[c(1, 3)] <- c(6, 8)
deep_channels <- probe$contact[probe$layer_category == "deep"]



# Data exploration --------------------------------------------------------

# LFP of a sample trial over time for all cortical channels over time
par(mfrow = c(4, 4))
for (i in 1:dim(LFP)[1]) {
  plot(time[100, ], LFP[i, , 100], type = "l", main = paste("LFP of a Sample Trial in Channel", i))
}
par(mfrow = c(1, 1))

# LFP of a sample trial over time for all cortical channels over time in one plot
colours <- RColorBrewer::brewer.pal(15, "Greens")
plot(time[100, ], LFP[1, , 100], type = "l", col = colours[1], main = "LFP of a Sample Trial in Different Channels")
for (i in 2:dim(LFP)[1]) {
  lines(time[100, ], LFP[i, , 100], col = colours[i])
}
legend("topright", legend = paste("Channel", 1:dim(LFP)[1]), col = colours, lty = 1, cex = 0.4)

# LFP of a sample trial over time for 3 cortical channels over time in one plot
colours = c("green", "blue", "violet")
plot(time[100, ], LFP[1, , 100], type = "l", lwd = 2, col = colours[1], main = "LFP of a Sample Trial in 3 Different Channels (Lower, Middle and Upper")
lines(time[100, ], LFP[8, , 100], lwd = 2, col = colours[2])
lines(time[100, ], LFP[15, , 100], lwd = 2, col = colours[3])
legend("bottomright", legend = paste("Channel", c(1, 8, 15)), col = colours, lty = 1, cex = 0.8)

# # stationary LFP of a sample trial over time for 3 cortical channels over time in one plot
# colours = c("green", "blue", "violet")
# plot(time[100, ], LFP_stationary[1, , 100], type = "l", lwd = 2, col = colours[1], main = "LFP of a Sample Trial in 3 Different Channels (Lower, Middle and Upper")
# lines(time[100, ], LFP_stationary[8, , 100], lwd = 2, col = colours[2])
# lines(time[100, ], LFP_stationary[15, , 100], lwd = 2, col = colours[3])
# legend("bottomright", legend = paste("Channel", c(1, 8, 15)), col = colours, lty = 1, cex = 0.8)
