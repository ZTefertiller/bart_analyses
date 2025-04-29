per_participant <- function(number) {
  
  # Select the participant index
  sub_index <- number
  
  # Number of trials for that participant
  n_trials <- sum(df$participant_id == unique(df$participant_id)[sub_index])
  
  # Get that participant's data
  sub_data <- df[df$participant_id == unique(df$participant_id)[sub_index], ]
  
  # Get omega_out posterior mean for that participant
  omega_mean <- colMeans(post$omega_out[, sub_index, ])
  
  # Set up balloon color coding
  colors <- c("blue", "orange", "yellow")
  balloon_colors <- factor(sub_data$balloon_color, levels = 1:3, labels = colors)
  
  # Start plotting
  plot(omega_mean, type = "l", col = "black", lwd = 2,
       ylab = "Omega", xlab = "Trial",
       ylim = range(c(omega_mean, sub_data$inflations)),
       main = paste("Participant", unique(sub_data$participant_id)))
  
  # Shade background by balloon color
  rects <- rle(as.character(balloon_colors))
  
  start <- 1
  for (i in seq_along(rects$lengths)) {
    end <- start + rects$lengths[i] - 1
    col_shade <- ifelse(rects$values[i] == "blue", rgb(0,0,1,0.1),
                        ifelse(rects$values[i] == "orange", rgb(1,0.5,0,0.1),
                               rgb(1,1,0,0.1)))
    rect(xleft = start, xright = end, ybottom = -Inf, ytop = Inf,
         col = col_shade, border = NA)
    start <- end + 1
  }
  
  # Redraw omega and actual inflations on top
  lines(omega_mean, col = "black", lwd = 2)
  lines(sub_data$inflations, col = "red", lwd = 2, lty = 2)
  
  # Mark reversal (trial 91)
  abline(v = 91, col = "darkred", lty = 3)
  
  # Add a legend (without the reversal line)
  legend("topright", legend = c("Omega", "Inflations"),
         col = c("black", "red"), lty = c(1, 2), bty = "n")
  
  pop_trials <- which(sub_data$popped == 1)  # trials where balloon popped
  points(x = pop_trials,
         y = sub_data$inflations[pop_trials],  # number of pumps when it popped
         pch = 4,  # x-mark symbol
         col = "darkred",
         cex = 0.7,  # size of the x-mark
         lwd = 1.2)  # thickness of the lines in the x
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

