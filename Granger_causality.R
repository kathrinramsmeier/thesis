
# VAR model creation ------------------------------------------------------

channel1_filtered <- t(channel1_filtered)

# select lag for VAR model using AIC
select_lag <- VARselect(channel1_filtered, type = "const")
select_lag$selection
plot(select_lag$criteria[1, ]) # try BIC!

# create the VAR model
VAR_model <- VAR(channel1_clipped_demeaned, type = "const", p = select_lag$selection)
summary(VAR_model)



# Granger causality -------------------------------------------------------

GC <- granger_causality(VAR_model)
GC