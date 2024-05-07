library(reticulate) # load npy data
np <- import("numpy")
library(tidyverse)  # data manipulation
library(vars)       # VAR model
library(tseries)    # ADF test
library(bruceR)     # Granger causality
library(multitaper) # Multitaper spectral analysis
library(igraph)     # graphs for GC

options(scipen = 1000)



# Session C190127 ---------------------------------------------------------

# LFP
LFP <- np$load("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\lfp_array_uv.npy")

# electrode channel
probe <- read.csv("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\probe.csv", header = FALSE)
probe_header <- read.csv("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\probe_header.csv", header = FALSE)
colnames(probe) <- probe_header
# rm(list = "probe_header")

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

# probe
probe <- read.csv("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\probe.csv", header = FALSE)
probe_header <- read.csv("C:\\Users\\ramsm\\Desktop\\Master\\Thesis\\data\\C190127-npy\\probe_header.csv", header = FALSE)
colnames(probe) <- probe_header
rm(list = "probe_header")



# Data inspection ---------------------------------------------------------

LFP[1:5, 1:5, 1:5] # electrode channel (1:15), time, trial
dim(LFP)

summary(time) # 0 is stimulus, before ore-stimulus, after past-stimulus
length(time) # 2 values less than LFP
dim(behaviour)

probe # layer 1:5 = upper, 6:10 = middle, 11:15 = deep



# Data preprocessing ------------------------------------------------------

# add time value at the start and at the end
time <- c(time[1] - diff(time[c(1, 2)]), time, rev(time)[1] + diff(rev(time)[c(1, 2)]))

# exclude incorrect trials and catch trials (where the colours were the same) from LFP and task
correct_trials_ind <- which(behaviour$accuracy_logical == 1)
non_catch_trails_ind <- task$trial_number_count[(task$catch_trial_logical == 0)]
trials_ind <- intersect(correct_trials_ind, non_catch_trails_ind)
LFP <- LFP[, , trials_ind]
task <- task[trials_ind, ] 
task$trial_number_count <- 1:nrow(task)

# # only investigate the LFP after stimulus onset and up to 600 ms after stimulus onset; clipping 10 ms after stimulus onset 
# time_ind <- which(time > 10)
# time <- time[time_ind]
# LFP <- LFP[, time_ind, ]

# clipping 10 ms before and after stimulus onset
time_ind <- which(time > 10 | time < -10)
time <- time[time_ind]
LFP <- LFP[, time_ind, ]



# Checks ------------------------------------------------------------------

dim(LFP)[1] == 15
dim(LFP)[2] == length(time)
dim(LFP)[3] == nrow(task)



# Necessary variables -----------------------------------------------------

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
  plot(time, LFP[i, , 100], type = "l", main = paste("LFP of a Sample Trial in Channel", i))
}
par(mfrow = c(1, 1))

# LFP of a sample trial over time for all cortical channels over time in one plot
bluePalette <- RColorBrewer::colorRampPalette(c("lightblue", "darkblue"))
colours <- bluePalette(15)
plot(time, LFP[1, , 100], type = "l", col = colours[1], main = "LFP of a Sample Trial in Different Channels")
for (i in 2:dim(LFP)[1]) {
  lines(time, LFP[i, , 100], col = colours[i])
}
legend("topright", legend = paste("Channel", 1:dim(LFP)[1]), col = colours, lty = 1, cex = 0.4)

# LFP of a sample trial over time for 3 cortical channels over time in one plot
colours = c("green", "blue", "violet")
plot(time, LFP[1, , 100], type = "l", lwd = 2, col = colours[1], main = "LFP of a Sample Trial in 3 Different Channels (Lower, Middle and Upper")
lines(time, LFP[8, , 100], lwd = 2, col = colours[2])
lines(time, LFP[15, , 100], lwd = 2, col = colours[3])
legend("bottomright", legend = paste("Channel", c(1, 8, 15)), col = colours, lty = 1, cex = 0.8)
