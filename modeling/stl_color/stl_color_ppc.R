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



run_ppc_all_participants <- function(df, post, output_folder = NULL, n_draws = 200, use_base_plot = FALSE) {
  # Load required packages
  if (!require("bayesplot")) {
    message("Installing bayesplot package...")
    install.packages("bayesplot")
    library(bayesplot)
  }
  
  # Ensure that graphics devices are clean
  graphics.off()
  
  participant_ids <- unique(df$participant_id)
  n_participants <- length(participant_ids)
  
  if (!is.null(output_folder)) {
    if (!dir.exists(output_folder)) {
      dir.create(output_folder, recursive = TRUE)
    }
  }
  
  # Print diagnostics about posterior samples
  cat("Checking posterior sample availability:\n")
  cat("b_beta_pre available:", "b_beta_pre" %in% names(post), "\n")
  cat("b_beta_post available:", "b_beta_post" %in% names(post), "\n")
  cat("o_beta_pre available:", "o_beta_pre" %in% names(post), "\n")
  cat("o_beta_post available:", "o_beta_post" %in% names(post), "\n")
  
  # Check expected parameters in Stan model output
  required_params <- c(
    "b_beta_pre", "b_beta_post", "o_beta_pre", "o_beta_post", "y_beta",
    "b_omegaone", "o_omegaone", "y_omegaone",
    "b_vwin_pre", "b_vwin_post", "b_vloss_pre", "b_vloss_post",
    "o_vwin_pre", "o_vwin_post", "o_vloss_pre", "o_vloss_post",
    "y_vwin", "y_vloss"
  )
  
  missing_params <- required_params[!required_params %in% names(post)]
  if (length(missing_params) > 0) {
    warning("Missing parameters in posterior samples: ", paste(missing_params, collapse=", "))
    cat("Available parameters: ", paste(names(post), collapse=", "), "\n")
    
    # Map parameters if necessary - this is a fallback if parameter names don't match
    if ("b_beta" %in% names(post) && !"b_beta_pre" %in% names(post)) {
      post$b_beta_pre <- post$b_beta
      post$b_beta_post <- post$b_beta
      cat("Using b_beta for both pre and post parameters\n")
    }
    if ("o_beta" %in% names(post) && !"o_beta_pre" %in% names(post)) {
      post$o_beta_pre <- post$o_beta
      post$o_beta_post <- post$o_beta
      cat("Using o_beta for both pre and post parameters\n")
    }
  }
  
  # Check data structure
  cat("Data structure check:\n")
  cat("Columns in df:", paste(names(df), collapse=", "), "\n")
  cat("Checking required columns...\n")
  required_cols <- c("participant_id", "balloon_color", "color_max", "opportunity", "inflations")
  for (col in required_cols) {
    cat(col, "exists:", col %in% names(df), "\n")
  }
  
  pb <- txtProgressBar(min = 0, max = n_participants, style = 3)
  for (i in seq_len(n_participants)) {
    
    sub_index <- i
    participant_id <- participant_ids[i]
    participant_data <- df[df$participant_id == participant_id, ]
    
    # Convert balloon colors if not already numeric
    if (is.character(participant_data$balloon_color)) {
      participant_data$balloon_color <- match(participant_data$balloon_color, c("b", "o", "y"))
      cat("Converted balloon colors for participant", participant_id, "\n")
    }
    
    stopifnot(all(participant_data$balloon_color %in% 1:3))
    n_trials <- nrow(participant_data)
    
    # Create a log file for debugging
    if (!is.null(output_folder)) {
      log_file <- file.path(output_folder, paste0("participant_", participant_id, "_log.txt"))
      write(paste("Starting simulation for participant", participant_id), log_file)
      write(paste("Number of trials:", n_trials), log_file, append = TRUE)
    }
    
    simulated_pumps <- matrix(NA, nrow = n_draws, ncol = n_trials)
    
    for (d in 1:n_draws) {
      # Get parameters from posterior samples
      # First check if we have pre/post parameters or just single parameters
      has_pre_post <- all(c("b_beta_pre", "b_beta_post") %in% names(post))
      
      if (has_pre_post) {
        beta_blue_pre   <- post$b_beta_pre[d, sub_index]
        beta_blue_post  <- post$b_beta_post[d, sub_index]
        beta_orange_pre <- post$o_beta_pre[d, sub_index]
        beta_orange_post<- post$o_beta_post[d, sub_index]
      } else {
        # Fallback if no pre/post distinction
        beta_blue_pre   <- post$b_beta[d, sub_index]
        beta_blue_post  <- post$b_beta[d, sub_index]
        beta_orange_pre <- post$o_beta[d, sub_index]
        beta_orange_post<- post$o_beta[d, sub_index]
      }
      
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
      
      # Track last observed trial for each color
      last_trial_blue <- 0
      last_trial_orange <- 0
      last_trial_yellow <- 0
      
      # Initialize omega for each color
      omega <- matrix(NA_real_, nrow = 3, ncol = n_trials)
      
      for (t in seq_len(n_trials)) {
        # Get current trial info
        col <- participant_data$balloon_color[t]
        nmax <- participant_data$color_max[t]
        
        # Safety check for nmax
        if (is.na(nmax) || nmax <= 0) {
          warning("Invalid nmax at trial ", t, " for participant ", participant_id)
          nmax <- 100  # Default fallback
        }
        
        # Determine if this is the first occurrence of this color
        is_first_occurrence <- FALSE
        if (col == 1 && last_trial_blue == 0) is_first_occurrence <- TRUE
        if (col == 2 && last_trial_orange == 0) is_first_occurrence <- TRUE
        if (col == 3 && last_trial_yellow == 0) is_first_occurrence <- TRUE
        
        if (is_first_occurrence) {
          # First occurrence of this color
          if (col == 1) {
            omega[col, t] <- omegaone_blue * nmax
          } else if (col == 2) {
            omega[col, t] <- omegaone_orange * nmax
          } else if (col == 3) {
            omega[col, t] <- omegaone_yellow * nmax
          }
          
          # Log this for the first draw
          if (d == 1 && !is.null(output_folder)) {
            write(paste("Trial", t, ": First occurrence of color", col, 
                        "omega =", omega[col, t]), log_file, append = TRUE)
          }
          
        } else {
          # Not the first occurrence - update based on previous trial of this color
          prev_trial <- if (col == 1) last_trial_blue else if (col == 2) last_trial_orange else last_trial_yellow
          
          prev_pumps <- simulated_pumps[d, prev_trial]
          prev_nmax <- participant_data$color_max[prev_trial]
          prev_opportunity <- participant_data$opportunity[prev_trial]
          prev_outcome <- as.integer(prev_pumps >= prev_opportunity)
          prev_omega <- omega[col, prev_trial]
          
          # Safety checks
          if (is.na(prev_nmax) || prev_nmax <= 0) prev_nmax <- 100
          
          if (col == 3) {
            vwin_i <- vwin_yellow
            vloss_i <- vloss_yellow
          } else if (col == 1) {
            vwin_i <- ifelse(t <= 90, vwin_blue_pre, vwin_blue_post)
            vloss_i <- ifelse(t <= 90, vloss_blue_pre, vloss_blue_post)
          } else if (col == 2) {
            vwin_i <- ifelse(t <= 90, vwin_orange_pre, vwin_orange_post)
            vloss_i <- ifelse(t <= 90, vloss_orange_pre, vloss_orange_post)
          }
          
          # Apply learning rule
          if (prev_outcome == 1) {
            omega[col, t] <- nmax * (prev_omega / prev_nmax) * 
              (1 - vloss_i * (1 - prev_pumps / prev_nmax))
          } else {
            omega[col, t] <- nmax * (prev_omega / prev_nmax) * 
              (1 + vwin_i * (prev_pumps / prev_nmax))
          }
          
          # Log updates for the first draw
          if (d == 1 && !is.null(output_folder)) {
            write(paste("Trial", t, ": Updated color", col, 
                        "omega =", omega[col, t],
                        "prev_outcome =", prev_outcome,
                        "prev_pumps =", prev_pumps), log_file, append = TRUE)
          }
        }
        
        # Simulate pumps after omega is set
        omega_k <- omega[col, t]
        
        # Safety check
        if (is.na(omega_k)) {
          warning("omega_k is NA at t = ", t, ", color = ", col, ", participant = ", participant_id)
          omega_k <- nmax/2  # Reasonable default
        }
        
        # Select appropriate beta based on trial number (pre/post reversal)
        if (t <= 90) {
          beta_i <- c(beta_blue_pre, beta_orange_pre, beta_yellow)[col]
        } else {
          beta_i <- c(beta_blue_post, beta_orange_post, beta_yellow)[col]
        }
        
        # Safety check for beta
        if (is.na(beta_i) || beta_i <= 0) {
          warning("Invalid beta at trial ", t, " for participant ", participant_id)
          beta_i <- 1.0  # Default fallback
        }
        
        # Calculate pump probabilities
        pump_probs <- plogis(-beta_i * ((1:nmax) - omega_k))
        
        # Ensure valid probabilities
        pump_probs <- pmin(pmax(pump_probs, 0.00001), 0.99999)
        
        # Simulate pumps
        simulated_pumps[d, t] <- rbinom(1, nmax, mean(pump_probs))
        
        # Update last trial index for this color
        if (col == 1) last_trial_blue <- t
        else if (col == 2) last_trial_orange <- t
        else if (col == 3) last_trial_yellow <- t
        
      } # end of trial loop
    } # end of draw loop
    
    # Check if simulations produced valid results
    if (all(is.na(simulated_pumps))) {
      warning("All simulated pumps are NA for participant ", participant_id)
      next
    }
    
    # Log simulation summary
    if (!is.null(output_folder)) {
      write(paste("Simulation complete. Mean simulated pumps:", mean(simulated_pumps, na.rm=TRUE)), 
            log_file, append = TRUE)
      write(paste("Mean observed pumps:", mean(participant_data$inflations, na.rm=TRUE)),
            log_file, append = TRUE)
    }
    
    # --- NOW plot/save after all draws simulated ---
    if (is.null(output_folder)) {
      # Just plot to screen
      if (use_base_plot) {
        # Use base R plotting
        observed <- participant_data$inflations
        plot(density(observed, na.rm=TRUE), 
             main=paste("Participant", participant_id),
             xlab="Number of Pumps", ylab="Density",
             lwd=2, col="black")
        
        # Add simulated densities
        for (d in 1:min(n_draws, 50)) {  # Limit to 50 lines for clarity
          lines(density(simulated_pumps[d,], na.rm=TRUE), col="gray70", lwd=0.5)
        }
        
        # Add legend
        legend("topright", legend=c("Observed", "Simulated"), 
               col=c("black", "gray70"), lwd=c(2, 0.5))
      } else {
        # Try using bayesplot
        tryCatch({
          ppc_dens_overlay(
            y = participant_data$inflations,
            yrep = simulated_pumps
          )
          title(main = paste("Participant", participant_id))
        }, error = function(e) {
          message("Error using bayesplot, falling back to base R plotting")
          plot(density(participant_data$inflations, na.rm=TRUE), 
               main=paste("Participant", participant_id),
               xlab="Number of Pumps", ylab="Density")
          for (d in 1:min(n_draws, 50)) {
            lines(density(simulated_pumps[d,], na.rm=TRUE), col="gray70")
          }
        })
      }
      
    } else {
      # --- here is the correct full saving block ---
      file_path <- file.path(output_folder, paste0("ppc_participant_", participant_id, ".png"))
      png(file_path, width = 800, height = 600)
      
      if (use_base_plot) {
        # Use base R plotting
        observed <- participant_data$inflations
        plot(density(observed, na.rm=TRUE), 
             main=paste("Participant", participant_id),
             xlab="Number of Pumps", ylab="Density",
             lwd=2, col="black")
        
        # Add simulated densities
        for (d in 1:min(n_draws, 50)) {
          lines(density(simulated_pumps[d,], na.rm=TRUE), col="gray70", lwd=0.5)
        }
        
        # Add legend
        legend("topright", legend=c("Observed", "Simulated"), 
               col=c("black", "gray70"), lwd=c(2, 0.5))
      } else {
        # Try using bayesplot
        tryCatch({
          ppc_dens_overlay(
            y = participant_data$inflations,
            yrep = simulated_pumps
          )
          title(main = paste("Participant", participant_id))
        }, error = function(e) {
          message("Error using bayesplot, falling back to base R plotting")
          plot(density(participant_data$inflations, na.rm=TRUE), 
               main=paste("Participant", participant_id),
               xlab="Number of Pumps", ylab="Density")
          for (d in 1:min(n_draws, 50)) {
            lines(density(simulated_pumps[d,], na.rm=TRUE), col="gray70")
          }
        })
      }
      
      dev.off()
      
      # Also save a diagnostic plot showing side-by-side boxplots
      diag_file_path <- file.path(output_folder, paste0("diag_participant_", participant_id, ".png"))
      png(diag_file_path, width = 800, height = 600)
      
      # Create side-by-side boxplot comparing observed vs simulated
      par(mfrow=c(1,2))
      boxplot(participant_data$inflations, main="Observed", ylab="Number of Pumps")
      boxplot(as.vector(simulated_pumps), main="Simulated", ylab="Number of Pumps")
      par(mfrow=c(1,1))
      
      dev.off()
    }
    setTxtProgressBar(pb, i)
  } # end of participant loop
  close(pb)
  
  cat("\nPPC analysis complete! Check the output folder for results.\n")
}

run_ppc_all_participants(
  df = df,
  post = post,
  output_folder = "~/Desktop/participant_checks_ppc",
  n_draws = 100  # or however many you want
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









