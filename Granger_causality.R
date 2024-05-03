
# look at the first channel
# LFP_stationary_ch1 <- LFP_stationary[1, , ]
LFP_stationary_ch1 <- channel1_differenced
dim(LFP_stationary_ch1) # time, trial

# filter all primed and unprimed trials
channel1_filtered_primed <- LFP_stationary_ch1[, primed_ind]
channel1_filtered_unprimed <- LFP_stationary_ch1[, unprimed_ind]



# VAR model creation ------------------------------------------------------

# select lag for VAR model using AIC
select_lag <- VARselect(channel1_filtered_primed) # not converging
select_lag$selection
plot(select_lag$criteria[1, ]) # try BIC!

# create the VAR model
VAR_model_channel1_primed <- VAR(channel1_filtered_primed, type = "const", p = select_lag$selection)
summary(VAR_model_channel1_primed)



# Granger causality -------------------------------------------------------

# H0: time series X does not cause time series Y to Granger-cause itself

GC <- granger_causality(VAR_model)
GC