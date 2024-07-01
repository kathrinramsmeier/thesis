
# VAR Model Creation ------------------------------------------------------

# look at a sample trial as an example
LFP_sample_trial <- LFP_stationary[[333]]
LFP_sample_trial <- do.call(cbind, lapply(LFP_sample_trial, function(x) as.numeric(x)))
dim(LFP_sample_trial) # time, channel
colnames(LFP_sample_trial) <- paste0("channel", 1:15)

# select lag for VAR model using AIC
select_lag <- VARselect(LFP_sample_trial, type = "const", lag.max = 30) 
select_lag$selection
plot(select_lag$criteria[1, ])

# create the VAR model
VAR_model_sample_trial <- vars::VAR(LFP_sample_trial, type = "const", p = 2) # try with lower lag order because the model is to complex otherwise for that "small" amount of observations compared to 15 variables https://www.frontiersin.org/articles/10.3389/fncom.2013.00159/full#F2



# Test Whether the VAR Model Adequately Captures the Correlation  ---------
# Structure in the Data ---------------------------------------------------

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

# 4.) Portmanteau- and Breusch-Godfrey test (not mentioned in the paper)
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

# 5.) investigate residuals (not mentioned in the paper)
plot_residuals <- function(res) {
  par(mfrow = c(3, 5))  
  for (i in 1:15) {  
    plot(res[, i], type = "p", xlab = "Time", ylab = "Residuals")
  }
  par(mfrow = c(1, 1))
}


head(residuals(VAR_model_sample_trial))
BIC(VAR_model_sample_trial)
calculate_adj_sum_square_err(VAR_model = VAR_model_sample_trial) # values less then 0.3 signify that the VAR model may not have captured the data adequately
check_consistency(VAR_model = VAR_model_sample_trial) # values below 80% may give cause fo concern
dW_VAR_test(VAR_model = VAR_model_sample_trial) # if d < 1 there may be cause for concern
P_BG_test(VAR_model = VAR_model_sample_trial) # significant autocorrelation present in the residuals
plot_residuals(res = residuals(VAR_model_sample_trial)) # they show pattern, which is not good but is way better with lag order 3 instead of 1
# -> indication that the VAR model is not capturing the data well



# Granger Causality -------------------------------------------------------

# H0: time series X does not cause time series Y to Granger-cause itself

GC <- granger_causality(VAR_model_sample_trial)
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



# Conditional Granger Causality -------------------------------------------

GC_p_values_sample_trial <- cbind(combinations, adj_p_values_F)

# prepare matrix with layer information and p values
layer_influencing <- rep(NA, length(adj_p_values_F))
layer_influencing[GC_p_values_sample_trial[, 1] %in% upper_channels] <- "upper"
layer_influencing[GC_p_values_sample_trial[, 1] %in% middle_channels] <- "middle"
layer_influencing[GC_p_values_sample_trial[, 1] %in% deep_channels] <- "deep"
layer_influenced <- rep(NA, nrow(GC_p_values_sample_trial))
layer_influenced[GC_p_values_sample_trial[, 2] %in% upper_channels] <- "upper"
layer_influenced[GC_p_values_sample_trial[, 2] %in% middle_channels] <- "middle"
layer_influenced[GC_p_values_sample_trial[, 2] %in% deep_channels] <- "deep"
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
  p_values_layer_level <- adj_p_values_F[ind]
  percentage_sign_layer_level[i] <- sum(p_values_layer_level < 0.05) / length(p_values_layer_level)
}

percentage_sign_layer_level <- cbind(percentage_sign, percentage_sign_layer_level)
percentage_sign_layer_level

causal_signals <- LFP_sample_trial[, deep_channels]
dependent_signals <- LFP_sample_trial[, middle_channels]
conditioning_signals <- LFP_sample_trial[, upper_channels]

cond_GC_analysis <- condGranger(
  data = cbind(causal_signals, dependent_signals, conditioning_signals),
  # number of signals causing the dependent signals
  nx = 3,
  # number of signals that are Granger caused
  ny = 3,
  # autoregressive order (how many past observations to include)
  order = 3
)
cond_GC_analysis$prob
# seems to be suspicious (all are always significant, also the package that contained this function was removed from CRAN)

# other approach:
# for grangers::Granger.conditional(), only 1 ts can be provided so useless in my case



# Visualisation -----------------------------------------------------------

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
  col = c(colour_range[1], colour_range[length(colour_range)]),
  lty = 1,
  lwd = 2,
  cex = 0.8,
  title = "Colour Legend"
)
  


# GC Analysis for All Trials ----------------------------------------------

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

# GC analysis for some sample trials (making sure that the number of primed and unprimed trials is equal)
num_samples <- 20
set.seed(42)
subset_trials_unprimed_ind <- sort(sample(unprimed_ind_hra, num_samples / 2))
subset_trials_primed_ind <- sort(sample(primed_ind_hra, num_samples / 2))
subset_trials_ind <- sort(c(subset_trials_unprimed_ind, subset_trials_primed_ind))
start_time <- Sys.time()
GC_p_values <- GC_analysis(LFP = LFP_stationary_hra[subset_trials_ind], lag_order = 4)
end_time <- Sys.time()
end_time - start_time

# estimation of time to run all (time in h)
((end_time - start_time) / num_samples * length(LFP_stationary_hra)) / 60 / 60 

# checks
GC_p_values$adj_sum_square_err
GC_p_values$consistency
GC_p_values$portmanteu_test_res
GC_p_values$bg_test_res
for (i in 1:num_samples) {
  plot_residuals((GC_p_values$res[[i]])) 
}
par(mfrow = c(1, 1))

# divide in primed and unprimed trials
subset_unprimed_ind <- which(subset_trials_ind %in% unprimed_ind_hra)
subset_primed_ind <- which(subset_trials_ind %in% primed_ind_hra)
GC_p_values_unprimed <- GC_p_values$p_values[, c(1, 2, subset_unprimed_ind + 2)]
GC_p_values_primed <- GC_p_values$p_values[, c(1, 2, subset_primed_ind + 2)]

# calculate the percentage of significant GC values for each channel combi (primed and unprimed condition)
percentage_sign_unprimed <- rowSums(GC_p_values_unprimed[, 3:ncol(GC_p_values_unprimed)] < 0.05) / 
  (ncol(GC_p_values_unprimed) - 2)
percentage_sign_primed <- rowSums(GC_p_values_primed[, 3:ncol(GC_p_values_primed)] < 0.05) / 
  (ncol(GC_p_values_primed) - 2)
percentage_sign <- cbind(GC_p_values_unprimed[, 1:2], percentage_sign_primed, percentage_sign_unprimed)

calculate_GC_p_values <- function(LFP, un_primed_ind, lag_order) {
  
  # filter either primed or unprimed trials
  LFP_trials <- LFP[un_primed_ind]
  
  # GC analysis
  GC_p_values <- GC_analysis(LFP = LFP_trials, lag_order = lag_order)
  
  return(GC_p_values)
  
}



# Network Visualisation of GC Influences ----------------------------------

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
      vertex.label.cex = 1.2,
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
      col = c(colour_range[1], colour_range[length(colour_range)]),
      lty = 1,
      lwd = 2,
      cex = 0.8
    )
    
  }
  
}

plot_channel_influences(percentage_sign = percentage_sign)



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

percentage_sign_layer_level <- transform_layer_levels(GC_p_values)

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
    vertex.label.cex = 1.5,
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


plot_layer_influences(percentage_sign = percentage_sign_layer_level, un_primed = "")



# GC Analysis on Windowed Data --------------------------------------------

# Kamisnki et. al. (helps with non-stationarity)

# (for unequal lengths, i.e., lists)

library(MARSS)

window_size <- 50
overlap <- 25 

# extract windowed sections of the trial data
extract_windows <- function(trial_data, window_size, overlap) {
  
  num_samples <- nrow(trial_data)
  num_channels <- ncol(trial_data)
  num_windows <- floor((num_samples - window_size) / overlap) + 1
  
  windowed_data <- array(NA, dim = c(window_size, num_channels, num_windows))
  
  for (i in 1:num_windows) {
    start <- (i - 1) * overlap + 1
    end <- start + window_size - 1
    windowed_data[, , i] <- trial_data[start:end, ]
  }
  
  return(windowed_data)
}

# extract windowed sections of the trial data
windowed_data <- extract_windows(LFP_sample_trial, window_size, overlap)
dim(windowed_data) # windowed signal, channel, number of windows

# check stationarity (ratio of stationary windowed LFP signals)
p_values <- matrix(nrow = dim(windowed_data)[2], ncol = dim(windowed_data)[3])
for (i in 1:dim(windowed_data)[2]) {
  for (j in 1:dim(windowed_data)[3]) {
    p_values[i, j] <- tseries::adf.test(windowed_data[, i, j])$p.value
  }
}
sum(p_values < 0.05) / length(p_values)

# as comparison (stationarity before)
p_values <- rep(NA, ncol(LFP_sample_trial))
for (i in 1:ncol(LFP_sample_trial)) {
  p_values[i] <- tseries::adf.test(LFP_sample_trial[, i])$p.value
}
sum(p_values < 0.05) / length(p_values)

# test: VAR models on one windowed section
w1 <- windowed_data[, , 1]
dim(w1) # (windowed) time, channel
var_1 <- vars::VAR(w1, type = "const", p = 1)

# checks
head(residuals(var_1)) # the model seems to overfit - the residuals are all very close to 0
AIC(var_1) # model complex compared to the amount of information in the data 
calculate_adj_sum_square_err(VAR_model = var_1) # values less then 0.3 signify that the VAR model may not have captured the data adequately
check_consistency(VAR_model = var_1) # values below 80% may give cause fo concern
dW_VAR_test(VAR_model = var_1) # if d < 1 there may be cause for concern

# Granger causality
GC_1 <- (granger_causality(var_1))$result
p_values_F <- round(GC_1$p.F[-(seq(from = 15, to = 225, by = 15))], 6)

# Bonferroni correct p-values from F-Test and Chi^2-Test for multiple comparisons
adj_p_values_F <- p.adjust(p_values_F, method = "bonferroni")

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

# extract significant GC combis
sign_ind <- which(adj_p_values_F < 0.05)
sign_causalities <- combinations[sign_ind, ]
sign_causalities <- c(t(sign_causalities))



### now for all windowed sections

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

sign_indicators <- matrix(nrow = nrow(combinations), ncol = dim(windowed_data)[3])

# fit VAR models to the windowed data (Kamisnki et. al.)
for (i in 1:dim(windowed_data)[3]) {
  
  # VAR model on windowed section
  w_i <- windowed_data[, , i]
  VAR_model_i <- vars::VAR(w_i, type = "const", p = 1)
  
  # Granger Causality
  GC_sample_trial <- (granger_causality(VAR_model_i))$result
  p_values_F <- round(GC_sample_trial$p.F[-(seq(from = 15, to = 225, by = 15))], 6)
  
  # Bonferroni correct p-values from F-Test and Chi^2-Test for multiple comparisons
  adj_p_values_F <- p.adjust(p_values_F, method = "bonferroni")
  
  # extract significant GC combis
  sign_indicator <- as.numeric(adj_p_values_F < 0.05)
  
  # add each indicator to combinations
  sign_indicators[, i] <- sign_indicator
  
}

# check whether all are significant
all_sign <- rowSums(sign_indicators) == ncol(sign_indicators)
combinations <- cbind(combinations, as.numeric(all_sign))



### now for all trials

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

stationarity_ratio <- rep(NA, length(LFP_stationary))

window_size <- 50  # in ms
overlap <- 25 

# loop over all trials
for (k in 1:length(LFP_stationary)) {
  
  LFP_sample_trial <- do.call(cbind, LFP_stationary[[k]])
  # dim(LFP_sample_trial) # time, channel
  
  # extract windowed sections of the trial data
  windowed_data <- extract_windows(LFP_sample_trial, window_size, overlap)
  # dim(windowed_data) # windowed signal, channel, number of windows
  
  sign_indicators <- matrix(nrow = nrow(combinations), ncol = dim(windowed_data)[3])

  # fit VAR models to the windowed data (Kamisnki et. al.)
  for (i in 1:dim(windowed_data)[3]) {

    # VAR model on windowed section
    w_i <- windowed_data[, , i]
    VAR_model_i <- vars::VAR(w_i, type = "const", p = 1)

    # Granger Causality
    GC_sample_trial <- (granger_causality(VAR_model_i))$result
    p_values_F <- round(GC_sample_trial$p.F[-(seq(from = 15, to = 225, by = 15))], 6)

    # Bonferroni correct p-values from F-Test and Chi^2-Test for multiple comparisons
    adj_p_values_F <- p.adjust(p_values_F, method = "bonferroni")

    # extract significant GC combis
    sign_indicator <- as.numeric(adj_p_values_F < 0.05)

    # add each indicator to combinations
    sign_indicators[, i] <- sign_indicator

  }

  # check whether all are significant
  all_sign <- rowSums(sign_indicators) == ncol(sign_indicators)
  combinations <- cbind(combinations, as.numeric(all_sign))

  # # print(paste("Trial", k, "of", length(LFP_stationary)), "trials completed.")
  
}



