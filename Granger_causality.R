
# VAR model creation ------------------------------------------------------

# filter all primed and unprimed trials
LFP_stationary_primed <- LFP_stationary[, , primed_ind]
LFP_stationary_unprimed <- LFP_stationary[, , unprimed_ind]

# look at the first primed trial as an example
LFP_stationary_trial1 <- LFP_stationary_primed[, , 1]
LFP_stationary_trial1 <- t(LFP_stationary_trial1) 
dim(LFP_stationary_trial1) # time, channel

# select lag for VAR model using AIC
select_lag <- VARselect(LFP_stationary_trial1, type = "const", lag.max = 30) 
select_lag$selection
plot(select_lag$criteria[1, ])

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

# Bonferroni correct p-values from F-Test and Chi^2-Test for multiple comparisons
p_values_F <- round(GC_trial1$p.F[-(seq(from = 15, to = 225, by = 15))], 6)
adj_p_values_F <- p.adjust(p_values_F, method = "bonferroni")
p_values_Chisq <- round(GC_trial1$p.Chisq[-(seq(from = 15, to = 225, by = 15))], 6)
adj_p_values_Chisq <- p.adjust(p_values_Chisq, method = "bonferroni")
GC_causality_p_values <- cbind(combinations, adj_p_values_F, adj_p_values_Chisq)

# extract significant GC combis
sign_ind <- which((adj_p_values_F < 0.05) &  (adj_p_values_Chisq < 0.05))
sign_causalities <- GC_causality_p_values[sign_ind, c(1, 2)]
sign_causalities <- c(t(sign_causalities))



# Visualisation incl. all channels ----------------------------------------

plot(
  graph(sign_causalities, directed = TRUE),
  layout = layout.circle,
  vertex.label.cex = 1.5,
  vertex.size = 30,
  edge.arrow.size = 0.5,
  main = "GC Relationships Between Cortical Layers for One Primed Trial"
)



# Visualisation incl. only layer levels -----------------------------------

# prepare matrix with layer information and p values
layer_influencing <- rep(NA, nrow(GC_causality_p_values))
layer_influencing[GC_causality_p_values[, 1] %in% upper_channels] <- "upper"
layer_influencing[GC_causality_p_values[, 1] %in% middle_channels] <- "middle"
layer_influencing[GC_causality_p_values[, 1] %in% deep_channels] <- "deep"
layer_influenced <- rep(NA, nrow(GC_causality_p_values))
layer_influenced[GC_causality_p_values[, 2] %in% upper_channels] <- "upper"
layer_influenced[GC_causality_p_values[, 2] %in% middle_channels] <- "middle"
layer_influenced[GC_causality_p_values[, 2] %in% deep_channels] <- "deep"
GC_causality_p_values_layer_levels <- cbind(layer_influencing, layer_influenced, GC_causality_p_values)
GC_causality_p_values_layer_levels <- as.data.frame(GC_causality_p_values_layer_levels)

# calculate the percentage of significant GC combis
percentage_sign <- cbind(
  c("upper", "upper", "upper", "middle", "middle", "middle", "deep", "deep", "deep"),
  c("upper", "middle", "deep", "upper", "middle", "deep", "upper", "middle", "deep")
)
percentage_sign_layer_level <- rep(NA, nrow(percentage_sign))
for (i in 1:nrow(percentage_sign)) {
  ind <- GC_causality_p_values_layer_levels[, 1] %in% percentage_sign[i, 1] & GC_causality_p_values_layer_levels[, 2] %in% percentage_sign[i, 2]
  p_values_layer_level <- GC_causality_p_values_layer_levels[ind, ]$adj_p_values_F 
  percentage_sign_layer_level[i] <- sum(p_values_layer_level < 0.05) / length(p_values_layer_level)
}
percentage_sign <- cbind(percentage_sign, percentage_sign_layer_level)

# assign colours to vertices based on influence
min_influence <- min(percentage_sign[, 3])
max_influence <- max(percentage_sign[, 3])
colour_range <- cm.colors(10) 
influence_colours <- colour_range[cut(
  as.numeric(percentage_sign[, 3]),
  breaks = seq(min_influence, max_influence, length.out = length(colour_range) - 1),
  include.lowest = TRUE
)]
  
# plot the network with coloured nodes based on influence
sign_causalities_layer_levels <- c(percentage_sign[, c(1, 2)])
plot(
  graph(sign_causalities_layer_levels, directed = TRUE),
  layout = layout.circle,
  vertex.label.cex = 1.5,
  vertex.size = 30,
  edge.arrow.size = 0.5,
  edge.width = 3,
  edge.color = influence_colours,
  vertex.color = "white",
  main = paste0("GC Relationships Between Cortical Layer Levels for an Example Trial")
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
    select_lag <- VARselect(LFP_trial_i, lag.max = 30, type = "const") 
    VAR_model_trial_i <- VAR(LFP_trial_i, type = "const", p = select_lag$selection[1])
    
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

# GC analysis for some sample trials (making sure that the number of primed and unprimed trials is equal)
num_samples <- 4
set.seed(42)
subset_trials_unprimed_ind <- sort(sample(unprimed_ind, num_samples / 2))
subset_trials_primed_ind <- sort(sample(primed_ind, num_samples / 2))
subset_trials_ind <- sort(c(subset_trials_unprimed_ind, subset_trials_primed_ind))
start_time <- Sys.time()
GC_p_values <- GC_analysis(LFP = LFP_stationary[, , subset_trials_ind])
end_time <- Sys.time()
end_time - start_time

# estimation of time to run all
((end_time - start_time) / num_samples * dim(LFP_stationary)[3]) / 60 # would take 16 hours (for just one session)

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



# Network visualisation of GC influences ----------------------------------

plot_channel_influences <- function(percentage_sign) {
  
  # for header
  un_primed <- c(NA, NA, "Unprimed", "Primed")
  
  # for both primed and unprimed trials
  for (i in c(3, 4)) {
    
    # assign colours to vertices based on influence
    min_influence <- min(percentage_sign[, i])
    max_influence <- max(percentage_sign[, i])
    colour_range <- cm.colors(10) 
    influence_colours <- colour_range[cut(
      percentage_sign[, i],
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
      edge.width = 2,
      edge.color = influence_colours,
      vertex.color = "white",
      main = paste0("GC Relationships Between Cortical Layers for a Subset of ", un_primed[i], " Trials")
    )
    legend(
      "bottomright",
      legend = c("Smaller Influence", "Bigger Influence"),
      col = c(color_range[1], color_range[length(color_range)]),
      lty = 1,
      lwd = 2,
      cex = 0.8
    )
    
  }
  
}

plot_channel_influences(percentage_sign = percentage_sign)



# Network visualisation of GC influences (layer levels) -------------------

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
percentage_sign_layer_level_unprimed <- rep(NA, nrow(percentage_sign))
for (i in 1:nrow(percentage_sign)) {
  ind <- GC_p_values_layer_levels[, 1] %in% percentage_sign[i, 1] & GC_p_values_layer_levels[, 2] %in% percentage_sign[i, 2]
  p_values_layer_level <- percentage_sign_unprimed[ind]
  percentage_sign_layer_level_unprimed[i] <- sum(p_values_layer_level < 0.05) / length(p_values_layer_level)
}
percentage_sign_layer_level_primed <- rep(NA, nrow(percentage_sign))
for (i in 1:nrow(percentage_sign)) {
  ind <- GC_p_values_layer_levels[, 1] %in% percentage_sign[i, 1] & GC_p_values_layer_levels[, 2] %in% percentage_sign[i, 2]
  p_values_layer_level <- percentage_sign_primed[ind]
  percentage_sign_layer_level_primed[i] <- sum(p_values_layer_level < 0.05) / length(p_values_layer_level)
}
percentage_sign <- cbind(percentage_sign, percentage_sign_layer_level_unprimed, percentage_sign_layer_level_primed)

# plot the influences
plot_layer_influences <- function(percentage_sign) {
  
  # for header
  un_primed <- c(NA, NA, "Unprimed", "Primed")
  
  # for both primed and unprimed trials
  for (i in c(3, 4)) {
    
    # assign colours to vertices based on influence
    min_influence <- min(as.numeric(percentage_sign[, i]))
    max_influence <- max(as.numeric(percentage_sign[, i]))
    colour_range <- cm.colors(10) 
    influence_colours <- colour_range[cut(
      as.numeric(percentage_sign[, i]),
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
      main = paste0("GC Relationships Between Cortical Layers for a Subset of ", un_primed[i], " Trials")
    )
    legend(
      "bottomright",
      legend = c("Smaller Influence", "Bigger Influence"),
      col = c(color_range[1], color_range[length(color_range)]),
      lty = 1,
      lwd = 2,
      cex = 0.8
    )
    
  }
  
}

plot_layer_influences(percentage_sign = percentage_sign)

