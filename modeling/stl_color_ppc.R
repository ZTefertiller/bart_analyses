per_participant <- function(number) {
  sub_index <- number
  
  # Get participant ID and their data
  pid <- unique(df$participant_id)[sub_index]
  sub_data <- df[df$participant_id == pid, ]
  n_trials <- nrow(sub_data)
  
  # Get omega posterior mean
  omega_mean <- colMeans(post$omega_out[, sub_index, 1:n_trials])
  
  # Map balloon_color codes to actual color names
  color_map <- c("b" = "blue", "o" = "orange", "y" = "gold")
  sub_data$balloon_color_code <- color_map[sub_data$balloon_color]
  
  # Dataframe for plotting
  plot_df <- data.frame(
    trial = 1:n_trials,
    omega = omega_mean,
    inflations = sub_data$inflations,
    color = sub_data$balloon_color_code,
    popped = sub_data$popped
  )
  
  # Plot
  ggplot(plot_df, aes(x = trial)) +
    # Lines connecting omega to inflations
    geom_segment(aes(y = omega, yend = inflations, xend = trial, color = color), alpha = 0.7) +
    
    # Omega points
    geom_line(aes(y = omega), color = "black", size = 1.2) +
    
    # Inflation points
    geom_point(aes(y = inflations, fill = color), shape = 21, size = 1, stroke = 0.3, color = "black")+
    # Vertical line at trial 91
    #geom_vline(xintercept = 91, linetype = "dashed", color = "darkred") +
    
    # Popped balloon markers
    geom_point(data = subset(plot_df, popped == 1),
               aes(x = trial, y = inflations),
               shape = 4, size = .5, stroke = 1, color = "darkred") +
    
    scale_fill_identity() +
    scale_color_identity() +
    
    labs(title = paste("Participant", pid),
         x = "Trial", y = "Value") +
    coord_cartesian(ylim = c(0, 100)) +         # <--- add this line
    theme_minimal()
}





run_ppc_all_participants <- function(df, post, output_folder = NULL, n_draws = 200) {
  
  participant_ids <- unique(df$participant_id)
  n_participants <- length(participant_ids)
  
  if (!is.null(output_folder)) {
    if (!dir.exists(output_folder)) {
      dir.create(output_folder)
    }
  }
  
  pb <- txtProgressBar(min = 0, max = n_participants, style = 3)
  for (i in seq_len(n_participants)) {
    
    sub_index <- i
    participant_id <- participant_ids[i]
    participant_data <- df[df$participant_id == participant_id, ]
    participant_data$balloon_color <- match(participant_data$balloon_color, c("b", "o", "y"))
    stopifnot(all(participant_data$balloon_color %in% 1:3))
    n_trials <- nrow(participant_data)
    
    simulated_pumps <- matrix(NA, nrow = n_draws, ncol = n_trials)
    
    for (d in 1:n_draws) {
      
      beta_blue   <- post$b_beta[d, sub_index]
      beta_orange <- post$o_beta[d, sub_index]
      beta_yellow <- post$y_beta[d, sub_index]
      
      omegaone_blue   <- post$b_omegaone[d, sub_index]
      omegaone_orange <- post$o_omegaone[d, sub_index]
      omegaone_yellow <- post$y_omegaone[d, sub_index]
      
      vwin_blue_pre  <- post$b_vwin_pre[d, sub_index]
      vwin_blue_post <- post$b_vwin_post[d, sub_index]
      vloss_blue_pre <- post$b_vloss_pre[d, sub_index]
      vloss_blue_post<- post$b_vloss_post[d, sub_index]
      
      vwin_orange_pre  <- post$o_vwin_pre[d, sub_index]
      vwin_orange_post <- post$o_vwin_post[d, sub_index]
      vloss_orange_pre <- post$o_vloss_pre[d, sub_index]
      vloss_orange_post<- post$o_vloss_post[d, sub_index]
      
      vwin_yellow  <- post$y_vwin[d, sub_index]
      vloss_yellow <- post$y_vloss[d, sub_index]
      
      omega <- matrix(NA_real_, nrow = 3, ncol = n_trials)
      
      for (t in seq_len(n_trials)) {
        
        col <- participant_data$balloon_color[t]
        nmax <- participant_data$color_max[t]
        
        if (t == 1 || participant_data$balloon_color[t-1] != col) {
          # first trial OR color switched -> initialize
          if (col == 1) {
            omega[col, t] <- omegaone_blue * nmax
          } else if (col == 2) {
            omega[col, t] <- omegaone_orange * nmax
          } else if (col == 3) {
            omega[col, t] <- omegaone_yellow * nmax
          }
        } else {
          # normal update
          prev_pumps <- simulated_pumps[d, t-1]
          
          ## --- new safe check ---
          if (is.na(prev_pumps)) {
            stop("Previous pumps are NA at t = ", t)
          }
          
          prev_outcome <- as.integer(prev_pumps >= participant_data$opportunity[t-1])
          prev_omega <- omega[col, t-1]
          
          if (col == 3) { 
            vwin_i <- vwin_yellow
            vloss_i <- vloss_yellow
          } else if (col == 1) {
            vwin_i <- ifelse(t < 91, vwin_blue_pre, vwin_blue_post)
            vloss_i <- ifelse(t < 91, vloss_blue_pre, vloss_blue_post)
          } else if (col == 2) {
            vwin_i <- ifelse(t < 91, vwin_orange_pre, vwin_orange_post)
            vloss_i <- ifelse(t < 91, vloss_orange_pre, vloss_orange_post)
          }
          
          if (prev_outcome == 1) {
            omega[col, t] <- prev_omega * (1 - vloss_i * (1 - prev_pumps / nmax))
          } else {
            omega[col, t] <- prev_omega * (1 + vwin_i * (prev_pumps / nmax))
          }
        }
        
        # simulate pumps after omega is set
        omega_k <- omega[col, t]
        
        beta_i <- c(beta_blue, beta_orange, beta_yellow)[col]
        
        pump_probs <- plogis(-beta_i * ((1:nmax) - omega_k))
        simulated_pumps[d, t] <- rbinom(1, nmax, mean(pump_probs))
      } # end of trial loop
    } # end of draw loop
    
    # --- NOW plot/save after all draws simulated ---
    if (is.null(output_folder)) {
      # Just plot to screen
      ppc_dens_overlay(
        y = participant_data$inflations,
        yrep = simulated_pumps
      )
      title(main = paste("Participant", participant_id))
      
    } else {
      # --- here is the correct full saving block ---
      file_path <- file.path(output_folder, paste0("ppc_participant_", participant_id, ".png"))
      png(file_path, width = 800, height = 600)
      plot.new()  # <--- ADD THIS IMMEDIATELY AFTER png()
      ppc_dens_overlay(
        y = participant_data$inflations,
        yrep = simulated_pumps
      )
      title(main = paste("Participant", participant_id))
      dev.off()
    }
    setTxtProgressBar(pb, i)
  } # end of participant loop
  close(pb)
}


run_ppc_all_participants(
  df = df,
  post = post,
  output_folder = "~/Desktop/participant_checks_ppc",
  n_draws = 200  # or however many you want
)










## correlation

correlation <- function(df, post) {
  participant_ids <- unique(df$participant_id)
  cor_df <- data.frame(participant_id = character(),
                       correlation = numeric(),
                       stringsAsFactors = FALSE)
  
  for (i in seq_along(participant_ids)) {
    pid <- participant_ids[i]
    
    sub_data <- df %>% filter(participant_id == pid)
    omega_mean <- colMeans(post$omega_out[, i, ])  # i assumes post is ordered same as df
    
    # If the omega vector length doesn't match trial count, warn
    if (length(omega_mean) != nrow(sub_data)) {
      warning(paste("Length mismatch for participant", pid))
      next
    }
    
    # Correlation between model-predicted omega and actual inflations
    r <- cor(omega_mean, sub_data$inflations, method = "pearson", use = "complete.obs")
    
    cor_df <- rbind(cor_df, data.frame(participant_id = pid, correlation = r))
  }
  
  cor_df <- cor_df %>% mutate(z = atanh(correlation))
  
  # One-sample t-test
  t_test_result <- t.test(cor_df$z, mu = 0)
  
  print(t_test_result)
  
  hist(cor_df$correlation,
       breaks = 15,
       main = "Per-Participant Correlations (Predicted Omega vs Inflations)",
       xlab = "Pearson Correlation Coefficient")
  
}
library(ggplot2)









