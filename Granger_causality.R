
# VAR model creation ------------------------------------------------------

# filter all primed and unprimed trials
LFP_stationary_primed <- LFP_stationary[, , primed_ind]
LFP_stationary_unprimed <- LFP_stationary[, , unprimed_ind]

# look at the first primed trial
LFP_stationary_trial1 <- LFP_stationary_primed[, , 1]
LFP_stationary_trial1 <- t(LFP_stationary_trial1) 
dim(LFP_stationary_trial1) # time, channel

# select lag for VAR model using AIC
select_lag <- VARselect(LFP_stationary_trial1, type = "const", lag.max = 30) 
select_lag$selection
plot(select_lag$criteria[1, ]) # 16 seems to be enough already

# create the VAR model
VAR_model_trial1 <- VAR(LFP_stationary_trial1, type = "const", p = 16)



# Granger causality -------------------------------------------------------

# H0: time series X does not cause time series Y to Granger-cause itself

GC <- granger_causality(VAR_model_trial1)
GC_trial1 <- GC$result

# generate matrix of channel combinations
combinations <- matrix(nrow = 15 * 15 - 15, ncol = 2)
k <- 1
for (i in 1:15) {
  for (j in 1:15) {
    if (i != j) {  # Exclude identical pairs
      combinations[k, 1] <- i
      combinations[k, 2] <- j
      k <- k + 1
    }
  }
}
combinations <- combinations[, c(2, 1)]

# extract p-values from F-Test and Chi^2-Test
p_values_F <- round(GC_trial1$p.F[-(seq(from = 15, to = 225, by = 15))], 6)
p_values_Chisq <- round(GC_trial1$p.Chisq[-(seq(from = 15, to = 225, by = 15))], 6)
GC_causality_p_values <- cbind(combinations, p_values_F, p_values_Chisq)

# extract significant GC combis
sign_ind <- which((p_values_F < 0.05) &  (p_values_Chisq < 0.05))
sign_causalities <- GC_causality_p_values[sign_ind, c(1, 2)]
sign_causalities <- c(t(sign_causalities))

# visualisation
g <- graph(sign_causalities, directed = TRUE)
plot(
  g,
  layout = layout.circle,
  vertex.label.cex = 1.5,
  vertex.size = 30,
  edge.arrow.size = 0.5,
  main = "GC Relationships Between Cortical Layers for One Primed Trial"
)

# thoughts:
# incorporating the strength of causality (larger F- or Chi^2-values)?



# Reverse GC test ---------------------------------------------------------

# # some tests
# temp <- as.data.frame(cbind(LFP_stationary[1, , 1], LFP_stationary[1, , 2]))
# colnames(temp) <- c("ts1", "ts2")
# reverse_GC_test <- granger_test(ts1 ~ ts2, data = temp, test.reverse = TRUE)
# 
# temp <- as.data.frame(VAR_model_trial1$y[, 1:2])
# reverse_GC_test <- granger_test(y1 ~ y2, data = temp, test.reverse = TRUE)
#
# not sure yet how to implement this



# GC analysis for all trials ----------------------------------------------

GC_analysis <- function(LFP) {
  
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
    select_lag <- VARselect(LFP_trial_i, type = "const") 
    VAR_model_trial_i <- VAR(LFP_trial_i, type = "const", p = select_lag$selection[1])
    
    # Granger causality
    GC <- granger_causality(VAR_model_trial_i)
    GC_trial_i <- GC$result
    
    # extract p-values from F-Test
    p_values_F <- round(GC_trial_i$p.F[-(seq(from = 15, to = 225, by = 15))], 6)
    
    GC_p_values[, i + 2] <- p_values_F
    
    print(paste("GC calculation of trial", i, "of", dim(LFP)[3], "succeeded."))
    
  }
  
  return(GC_p_values)
  
}

# GC analysis for some sample trials (making sure that the number of primed and unprimed trials is equal)
num_samples <- 30
set.seed(42)
subset_trials_unprimed_ind <- sort(sample(unprimed_ind, num_samples / 2))
subset_trials_primed_ind <- sort(sample(primed_ind, num_samples / 2))
subset_trials_ind <- sort(c(subset_trials_unprimed_ind, subset_trials_primed_ind))
start_time <- Sys.time()
GC_p_values <- GC_analysis(LFP = LFP_stationary[, , subset_trials_ind])
end_time <- Sys.time()
end_time - start_time

# estimation of time to run all
((end_time - start_time) / num_samples * dim(LFP_stationary)[3]) / 60 # would take 22 hours (for just one session)

# divide in primed and unprimed trials
subset_unprimed_ind <- which(subset_trials_ind %in% unprimed_ind)
subset_primed_ind <- which(subset_trials_ind %in% primed_ind)
GC_p_values_unprimed <- GC_p_values[, c(1, 2, subset_unprimed_ind + 2)]
GC_p_values_primed <- GC_p_values[, c(1, 2, subset_primed_ind + 2)]

# calculate the percentage of significant GC values for each channel combi (primed and unprimed condition)
percentage_sign_unprimed <- rowSums(GC_p_values_unprimed[, 3:ncol(GC_p_values_unprimed)] < 0.05) / 
  (ncol(GC_p_values_unprimed) - 2)
percentage_sign_primed <- rowSums(GC_p_values_primed[, 3:ncol(GC_p_values_primed)] < 0.05) / 
  (ncol(GC_p_values_primed) - 2)

percentage_sign <- cbind(GC_p_values_unprimed[, 1:2], percentage_sign_primed, percentage_sign_unprimed)
percentage_sign



# Network visualisation of GC influences ----------------------------------

plot_networks <- function(percentage_sign) {
  
  # for header
  un_primed <- c(NA, NA, "Unprimed", "Primed")
  
  # for both primed and unprimed trials
  for (i in c(3, 4)) {
    
    # assign colours to vertices based on influence
    min_influence <- min(percentage_sign[, i])
    max_influence <- max(percentage_sign[, i])
    color_range <- cm.colors(10) 
    influence_colors <- color_range[cut(
      percentage_sign[, i],
      breaks = seq(min_influence, max_influence, length.out = length(color_range) - 1),
      include.lowest = TRUE
    )]
    
    # create a graph object
    g <- graph_from_data_frame(GC_causality_p_values[, 1:2], directed = TRUE)
    
    # plot the network with coloured nodes based on influence
    plot(
      g,
      layout = layout.circle,
      vertex.label.cex = 1.5,
      vertex.size = 30,
      edge.arrow.size = 0.5,
      edge.width = 2,
      edge.color = influence_colors,
      vertex.color = "white",
      main = paste0("GC Relationships Between Cortical Layers for a Subset of ", un_primed[i], " Trials")
    )
    legend(
      "bottomright",
      legend = c("Smaller Influence", "Bigger Influence"),
      col = c(color_range[1], color_range[length(color_range)]),
      lty = 1,
      lwd = 2,
      cex = 0.8,
      title = "Colour Legend"
    )
    
  }
  
}

plot_networks(percentage_sign = percentage_sign)





