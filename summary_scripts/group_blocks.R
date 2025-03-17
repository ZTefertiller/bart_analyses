library(dplyr)
library(bayesplot)
library(shinystan)
library(ggplot2)
library(rstan)
library(pryr)
library(loo)

# initial setup for my own task design, adapt balloon_color in the code and dataset as needed.

dataset <- read.csv('/Users/zachtefertiller/Desktop/Bart_Jatos_Results_Reversal_180/balloon_task_clean_data.csv', header=TRUE, sep=",", stringsAsFactors=FALSE)

# setup dataframe for reversal
df <- dataset %>%
  mutate(color_max = case_when(
    balloon_color == 'y' ~ 32,
    balloon_color == 'o' & trial_number > 91 ~ 128,
    balloon_color == 'o' & trial_number < 91 ~ 8,
    balloon_color == 'b' & trial_number < 91 ~ 128,
    balloon_color == 'b' & trial_number >= 91 ~ 8
  ))

print_df_info <- function(df) {
  df_name <- deparse(substitute(df))
  number_of_participants <- length(unique(df$participant_id))
  number_of_trials <- nrow(df) / number_of_participants
  cat("Dataframe used:", df_name, "\n")
  cat("Number of trials per participant:", number_of_trials, "\n")
  cat("Number of participants with criteria:", number_of_participants, "\n")
  return(list(
    number_of_participants = number_of_participants,
    number_of_trials = number_of_trials
  ))
}


# run this if you want to include spq, pdi, and caps in model

# dropping people who lack all questionnaires
df <- df %>%
  filter(!is.na(caps_total) & !is.na(pdi_total))


# summarize questionnaires, z-score, put in other dataframe for when we call
df_questionnaires <- df %>%
  group_by(participant_id) %>%
  summarise(
    spq_total  = first(spq_total_score),
    caps_total = first(caps_total),
    pdi_total  = first(pdi_total),
    pdi_distress  = first(pdi_distress),
    pdi_frequency  = first(pdi_frequency),
    pdi_conviction  = first(pdi_conviction),
    caps_distress = first(caps_distress),
    caps_frequency = first(caps_frequency),
    caps_intrusiveness = first(caps_intrusiveness),
    phq_total = first(phq_total),  
    mdq_total = first(mdq_total),
    ppgm_total = first(ppgm_total),
    ipip_total = first(ipip_total_score),
    
  ) %>%
  mutate(
    spq_z  = as.vector(scale(spq_total, center = TRUE, scale = TRUE)),
    caps_z = as.vector(scale(caps_total, center = TRUE, scale = TRUE)),
    pdi_z  = as.vector(scale(pdi_total, center = TRUE, scale = TRUE))
  )

df <- df %>% 
  left_join(df_questionnaires, by = "participant_id")

questionnaires <- df_questionnaires %>% 
  arrange(participant_id)


# separate data for when blue and orange are the high value balloon. yellow is constant
df_blue <- df %>% 
  filter(balloon_color == "b", color_max != 8) %>%
  arrange(participant_id, trial_number)

df_orange <- df %>% 
  filter(balloon_color == "o", color_max != 8) %>%
  arrange(participant_id, trial_number)

df_yellow <- df %>%
  filter(balloon_color == "y") %>%
  arrange(participant_id,trial_number)

# quantiles <- quantile(df_blue$spq_total_score, probs = c(0.1, 0.9), na.rm = TRUE)

# bottom_cutoff <- quantiles[1]
# top_cutoff <- quantiles[2]

mean_spq <- mean(df_blue$spq_total_score, na.rm = TRUE)
sd_spq   <- sd(df_blue$spq_total_score, na.rm = TRUE)

bottom_cutoff <- mean_spq - .80 * sd_spq
top_cutoff    <- mean_spq + .80 * sd_spq


df_blue_high_spq <- df_blue %>%
  # filter(spq_total_score >= 10)
  filter(spq_total_score >= top_cutoff)

df_blue_low_spq <- df_blue %>%
  # filter(spq_total_score <= 10)
  filter(spq_total_score <= bottom_cutoff)

df_orange_high_spq <- df_orange %>%
  # filter(spq_total_score >= 10)
  filter(spq_total_score >= top_cutoff)

df_orange_low_spq <- df_orange %>%
  # filter(spq_total_score <= 10)
  filter(spq_total_score <= bottom_cutoff)

df_yellow_high_spq <- df_yellow %>%
  # filter(spq_total_score <= 10)
  filter(spq_total_score >= top_cutoff)

df_yellow_low_spq <- df_yellow %>%
  # filter(spq_total_score <= 10)
  filter(spq_total_score <= bottom_cutoff)



print(cat("high spq scorers:", length(unique(df_blue_high_spq$participant_id)), "\n"))
print(cat("sanity check # of trials per person: ",
          ((nrow(df_blue_high_spq))/(length(unique(df_blue_high_spq$participant_id))))))

print(cat("low spq scorers:", length(unique(df_blue_low_spq$participant_id)), "\n"))
print(cat("sanity check # of trials per person: ", ((nrow(df_blue_low_spq))/(length(unique(df_blue_low_spq$participant_id))))))


# --- Compute block grouping and participant averages for each group ---

df_high_blue <- df_blue_high_spq %>%
  arrange(participant_id, trial_number) %>%      # ensure proper order
  group_by(participant_id) %>%
  mutate(block = ceiling(row_number() / 10)) %>%   # blocks: 1,2,3
  group_by(participant_id, block) %>%
  summarise(avg_inflations = mean(inflations, na.rm = TRUE), .groups = "drop") %>%
  mutate(spq_group = "High SPQ", balloon = "blue",
         final_block = block)  # keep blocks 1-3

# For Low SPQ:
df_low_blue <- df_blue_low_spq %>%
  arrange(participant_id, trial_number) %>%      
  mutate(block = ceiling(row_number() / 10)) %>%   
  group_by(participant_id, block) %>%
  summarise(avg_inflations = mean(inflations, na.rm = TRUE), .groups = "drop") %>%
  mutate(spq_group = "Low SPQ", balloon = "blue",
         final_block = block)

### --- ORANGE BALLOON PROCESSING (Blocks 4-6) ---
# For High SPQ:
df_high_orange <- df_orange_high_spq %>%
  arrange(participant_id, trial_number) %>%      # ensure proper order
  group_by(participant_id) %>%
  mutate(block = ceiling(row_number() / 10)) %>%   # blocks: 1,2,3 for orange trials
  group_by(participant_id, block) %>%
  summarise(avg_inflations = mean(inflations, na.rm = TRUE), .groups = "drop") %>%
  mutate(spq_group = "High SPQ", balloon = "orange",
         final_block = block + 3)  # shift to 4,5,6

# For Low SPQ:
df_low_orange <- df_orange_low_spq %>%
  arrange(participant_id, trial_number) %>%      
  mutate(block = ceiling(row_number() / 10)) %>%   
  group_by(participant_id, block) %>%
  summarise(avg_inflations = mean(inflations, na.rm = TRUE), .groups = "drop") %>%
  mutate(spq_group = "Low SPQ", balloon = "orange",
         final_block = block + 3)

### --- Combine BLUE and ORANGE participant means ---
participant_means_all <- bind_rows(df_high_blue, df_low_blue,
                                   df_high_orange, df_low_orange)

# Create a new grouping variable that combines balloon and SPQ group
participant_means_all <- participant_means_all %>%
  mutate(group = case_when(
    balloon == "blue" & spq_group == "Low SPQ" ~ "blue_Low",
    balloon == "blue" & spq_group == "High SPQ" ~ "blue_High",
    balloon == "orange" & spq_group == "Low SPQ" ~ "orange_Low",
    balloon == "orange" & spq_group == "High SPQ" ~ "orange_High"
  ))

### --- Compute significance for each final block comparing High vs. Low SPQ ---
sig_by_block <- participant_means_all %>% 
  group_by(final_block) %>% 
  summarise(
    p_value = wilcox.test(avg_inflations ~ spq_group)$p.value,
    max_infl = max(avg_inflations, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  mutate(
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    # Position significance text slightly above the highest average in each block
    y_position = 40
  )

### --- Plot the combined boxplot with significance annotations ---
combined_boxplot <- ggplot(participant_means_all, 
                           aes(x = factor(final_block), y = avg_inflations, 
                               fill = group, pattern = spq_group)) +
  # Patterned boxplots (dodged so that each SPQ group appears side-by-side)
  geom_boxplot_pattern(position = position_dodge(width = 0.75),
                       outlier.shape = NA,
                       color = "black",
                       pattern_fill = "black",
                       pattern_angle = 45,
                       pattern_density = 0.05,
                       pattern_spacing = 0.02) +
  # Add significance text above each final block
  geom_text(data = sig_by_block, 
            aes(x = factor(final_block), y = y_position, label = sig),
            vjust = 0, size = 5) +
  # Use manual scales for fill and pattern; now the fill scale preserves your blue colors
  scale_fill_manual(values = c("blue_Low" = "#4444FF", 
                               "blue_High" = "lightblue",
                               "orange_Low" = "darkorange", 
                               "orange_High" = "orange"),
                    name = "SPQ / Balloon") +
  scale_pattern_manual(values = c("Low SPQ" = "none", "High SPQ" = "stripe"),
                       name = "SPQ Group") +
  # Customize x-axis labels for all 6 final blocks
  scale_x_discrete(labels = c("1" = "random 1", 
                              "2" = "blue 2", 
                              "3" = "blue 3", 
                              "4" = "random 2", 
                              "5" = "orange 5", 
                              "6" = "orange 6")) +
  labs(title = "Participant Average Inflations by Block",
       x = "Block", y = "Average Inflations") +
  theme_minimal() +
  theme(text = element_text(family = "Courier New", face = "bold"),
        legend.position = "right")

print(combined_boxplot)

