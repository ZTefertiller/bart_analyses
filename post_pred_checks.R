stl_post_pred <- function(participant_data) {
  ntrials <- nrow(participant_data)
  nmax <- participant_data$optimal_inflations
  vwin <- participant_data$vwin[1]
  vloss <- participant_data$vloss[1]
  omegaone <- participant_data$omegaone[1]
  beta <- participant_data$beta[1]
  
  # Store outputs
  df_out <- data.frame(trial = 1:ntrials,
                       predicted_pumps = numeric(ntrials),
                       predicted_pop = numeric(ntrials),
                       omega = numeric(ntrials))
  
  # Initialize
  omega <- numeric(ntrials)
  predicted_pumps <- numeric(ntrials)
  predicted_pop <- numeric(ntrials)
  
  # First trial
  omega[1] <- omegaone * nmax[1]
  pump_probs <- plogis(-beta * ((1:nmax[1]) - omega[1]))
  predicted_pumps[1] <- rbinom(1, size = nmax[1], prob = mean(pump_probs))
  predicted_pop[1] <- as.integer(predicted_pumps[1] >= participant_data$opportunity[1])
  
  # Rest of the trials
  for (t in 2:ntrials) {
    prev_npump <- predicted_pumps[t - 1]
    prev_outcome <- predicted_pop[t - 1]
    nmax_t <- nmax[t]
    
    if (prev_outcome == 1) {
      omega[t] <- (omega[t - 1] / nmax_t) * 
        nmax_t * (1 - vloss * (1 - (prev_npump / nmax_t)))
    } else {
      omega[t] <- (omega[t - 1] / nmax_t) * 
        nmax_t * (1 + vwin * (prev_npump / nmax_t))
    }
    
    pump_probs <- plogis(-beta * ((1:nmax_t) - omega[t]))
    predicted_pumps[t] <- rbinom(1, size = nmax_t, prob = mean(pump_probs))
    predicted_pop[t] <- as.integer(predicted_pumps[t] >= participant_data$opportunity[t])
  }
  
  # Fill dataframe
  df_out$predicted_pumps <- predicted_pumps
  df_out$predicted_pop <- predicted_pop
  df_out$omega <- omega
  df_out$participant_id <- participant_data$participant_id[1]
  
  return(df_out)
}


main <- function(df) {
  output_list <- list()
  
  for (participant in unique(df$participant_id)) {
    participant_data <- df[df$participant_id == participant, ]
    sim_df <- stl_post_pred(participant_data)
    output_list[[as.character(participant)]] <- sim_df
  }
  
  do.call(rbind, output_list)
}



