library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

dataset <- read.csv("/Users/zachtefertiller/Desktop/Bart_Jatos_Results_Reversal_180/balloon_task_clean_data.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)



# Step 1: Calculate each participant's average task_time
participant_avg <- dataset %>%
  group_by(participant_id) %>%
  summarize(avg_task_time = mean(task_time, na.rm = TRUE)) %>%
  ungroup()

# Step 2: Compute overall mean and standard deviation of these averages
overall_mean <- mean(participant_avg$avg_task_time)
overall_sd   <- sd(participant_avg$avg_task_time)

# Step 3: Identify participant_ids that are NOT 2 SD below the overall average
valid_participants <- participant_avg %>%
  filter(avg_task_time >= (overall_mean - (1 * overall_sd))) %>%
  pull(participant_id)

# Step 4: Filter out rows from dataset for participant_ids that don't meet the threshold
filtered_dataset <- dataset %>%
  filter(participant_id %in% valid_participants)






ggplot(dataset, aes(x = spq_total_score)) +
  geom_histogram(binwidth = 1, color = "black", fill = "lightblue") +
  labs(x = "SPQ Total Score", y = "Number of Participants") +
  theme_minimal()


# Remove duplicate rows based on participant_id and spq_total_score
df_unique <- filtered_dataset %>% 
  distinct(participant_id, task_time, cigarettes_count, adjusted_inflations, total_score, ppgm_total, spq_total_score, pdi_total,caps_total,mdq_total,ipip_total_score)

# Count the number of unique participants for each spq_total_score
score_counts <- df_unique %>% 
  group_by(spq_total_score) %>% 
  summarise(count = n())

# Plot the counts as a bar chart
ggplot(score_counts, aes(x = spq_total_score, y = count)) +
  geom_bar(stat = "identity", fill = "navy", color = "black") +
  labs(x = "SPQ Total Score", y = "Number of Participants") +
  scale_y_continuous(breaks = seq(0,25,5), limits=c(0,16)) +
  theme_minimal()


# Count the number of unique participants for each spq_total_score
pdi_score_counts <- df_unique %>% 
  group_by(pdi_total) %>% 
  summarise(count = n())

# Plot the counts as a bar chart
ggplot(pdi_score_counts, aes(x = pdi_total, y = count)) +
  geom_bar(stat = "identity", fill = "darkred", color = "black") +
  labs(x = "PDI Score", y = "Number of Participants") +
  scale_y_continuous(breaks = seq(0,25,5), limits=c(0,16)) +
  theme_minimal()

# Count the number of unique participants for each spq_total_score
caps_score_counts <- df_unique %>% 
  group_by(caps_total) %>% 
  summarise(count = n())

ggplot(caps_score_counts, aes(x = caps_total, y = count)) +
  geom_bar(stat = "identity", fill = '#6f4685') +
  labs(x = "Caps Total Score", y = "Number of Participants") +
  scale_y_continuous(breaks = seq(0,50,5), limits = c(0,16)) +
  theme_minimal()


df_unique <- filtered_dataset %>% 
  distinct(participant_id, task_time, cigarettes_count, adjusted_inflations, total_score, 
           ppgm_total, spq_total_score, pdi_total, caps_total, mdq_total, ipip_total_score) %>%
  filter(!is.na(caps_total) & caps_total != "")


"magenta, violet, indigo, darkturquoise, darkred, navy, #005F6A, #660033, #052f2b, #C35214"
# Plot the counts as a bar chart
ggplot(caps_score_counts, aes(x = caps_total, y = count)) +
  geom_bar(stat = "identity", fill = "", color = "black") +
  labs(x = "Caps Total Score", y = "Number of Participants") +
  scale_y_continuous(breaks = seq(0,50,5), limits=c(0,16)) +
  theme_minimal()


ggplot(score_counts, aes(x = spq_total_score, y = count)) +
  geom_violin(fill = "navy") +
  labs(x = "", y = "SPQ Total Score") +
  theme_minimal() +
  coord_flip()  # flips so spq_total appears horizontally


# center and scale the spq
df_spq <- df %>%
  group_by(participant_id) %>%
  summarise(spq_total = first(spq_total_score))

df_spq <- df_spq %>%
  mutate(spq_z = as.vector(scale(spq_total, center = TRUE, scale = TRUE)))

df <- df %>% 
  left_join(df_spq, by = "participant_id")

library(dplyr)
library(ggplot2)

# Replace "your_participant_id" with the actual participant id
specific_id <- "NE03SUNS"

# Filter the dataframe for balloon_color "b" and the specific participant
df_filtered <- data %>%
  filter(balloon_color == "b", participant_id == specific_id)

# Create a line plot of inflations vs trial number
ggplot(df_filtered, aes(x = trial_number, y = inflations)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  labs(title = paste("Inflations for Participant", specific_id),
       x = "Trial Number",
       y = "Number of Inflations") +
  theme_minimal()



# 
# library(dplyr)
# library(rlang)
# 
# reg <- function(dataset, outcome, predictor) {
#   # Capture the column names as symbols
#   outcome_sym <- ensym(outcome)
#   predictor_sym <- ensym(predictor)
#   
#   # Optionally, convert symbols to strings for plotting labels
#   outcome_name <- as_string(outcome_sym)
#   predictor_name <- as_string(predictor_sym)
#   
#   # Count observations grouped by the outcome variable
#   score_counts <- dataset %>%
#     group_by(!!outcome_sym) %>%
#     summarise(count = n(), .groups = "drop")
#   
#   # Fit the linear model: outcome ~ predictor
#   # We build the formula dynamically using the column names.
#   model_formula <- as.formula(paste(outcome_name, "~", predictor_name))
#   model <- lm(model_formula, data = dataset)
#   
#   # Create a scatterplot:
#   # x-axis: predictor; y-axis: outcome
#   plot(dataset[[predictor_name]],
#        dataset[[outcome_name]],
#        xlab = predictor_name,
#        ylab = outcome_name)
#   
#   # Add the regression line
#   abline(model, col = "blue", lwd = 2)
#   
#   # Return the model summary
#   return(summary(model))
# }
# 
# 
# reg(df_unique, spq_total_score, adjusted_inflations)
# 
# 








# 
# 
# 
# 
# 
# 
# 
# clustering bs
# 
# # Select the numeric columns you want to cluster on
# data_to_cluster <- dataset %>% select(adjusted_inflations, task_time, spq_total_score)
# 
# # Scale the data so that all features contribute equally
# data_scaled <- scale(data_to_cluster)
# 
# # Set a seed for reproducibility and run k-means clustering (e.g., 3 clusters)
# set.seed(123)
# km_result <- kmeans(data_scaled, centers = 3, nstart = 25)
# 
# # Append the cluster assignment to the original dataset
# dataset$cluster <- as.factor(km_result$cluster)
# 
# # To visualize the clusters, we can reduce the dimensions using PCA:
# pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
# pca_df <- as.data.frame(pca_result$x)
# pca_df$cluster <- dataset$cluster
# 
# # Plot the first two principal components
# ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
#   geom_point(size = 3) +
#   labs(title = "Multidimensional K-Means Clustering (PCA Projection)",
#        x = "Principal Component 1", y = "Principal Component 2") +
#   theme_minimal()
# 
# # Assuming you have performed PCA on your scaled data:
# pca_res <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
# 
# # Create a biplot
# biplot(pca_res, scale = 0)
# 
# library(ggfortify)
# 
# autoplot(pca_res, data = dataset, colour = 'cluster',
#          loadings = TRUE, loadings.label = TRUE, 
#          loadings.label.size = 3) +
#   ggtitle("PCA Plot with Cluster Coloring")
# 
# 
# 
# 
# 
# 
# 
# 
# # Check for NAs
# sum(is.na(data_to_cluster))
# 
# data_to_cluster <- dataset %>% select(adjusted_inflations, pdi_total, caps_total, 
#                                       spq_total_score, ipip_total_score, ppgm_total,
#                                       phq_total, mdq_total)
# 
# 
# # Check for infinite values (after converting to a matrix)
# sum(!is.finite(as.matrix(data_to_cluster)))
# 
# # Replace infinite values with NA (if any)
# data_to_cluster[!is.finite(data_to_cluster)] <- NA
# 
# # Remove rows with any NA values
# data_to_cluster <- data_to_cluster %>% filter(complete.cases(.))
# 
# data_scaled <- scale(data_to_cluster)
# 
# 
# set.seed(123)
# km_result <- kmeans(data_scaled, centers = 3, nstart = 25)
# 
# # Append cluster assignment back to your (cleaned) dataset.
# # (Make sure your 'dataset' aligns with 'data_to_cluster' if you've removed rows)
# dataset$cluster <- as.factor(km_result$cluster)
# 
# # Perform PCA
# pca_result <- prcomp(data_scaled, center = TRUE, scale. = TRUE)
# library(ggfortify)
# 
# autoplot(pca_result, data = dataset, colour = 'cluster',
#          loadings = TRUE, loadings.label = TRUE, 
#          loadings.label.size = 3) +
#   ggtitle("PCA Plot with Cluster Coloring")
# 
# 


# Compute summary statistics for PDI Score
pdi_mean <- mean(df_unique$pdi_total, na.rm = TRUE)
pdi_sd <- sd(df_unique$pdi_total, na.rm = TRUE)

ggplot(pdi_score_counts, aes(x = pdi_total, y = count)) +
  geom_bar(stat = "identity", fill = "darkred", color = "black") +
  # Overlay the mean (dashed line)
  geom_vline(xintercept = pdi_mean, linetype = "dashed", color = "black", size = 1) +
  # Overlay one standard deviation boundaries (dotted lines)
  geom_vline(xintercept = pdi_mean - pdi_sd, linetype = "dotted", color = "black", size = 1) +
  geom_vline(xintercept = pdi_mean + pdi_sd, linetype = "dotted", color = "black", size = 1) +
  # Overlay two standard deviation boundaries (dotted lines)
  geom_vline(xintercept = pdi_mean - 2 * pdi_sd, linetype = "dotted", color = "black", size = 1) +
  geom_vline(xintercept = pdi_mean + 2 * pdi_sd, linetype = "dotted", color = "black", size = 1) +
  labs(x = "PDI Score", y = "Number of Participants") +
  scale_y_continuous(breaks = seq(0, 25, 5), limits = c(0, 16)) +
  theme_minimal()


# Compute summary statistics for CAPS Total Score
caps_mean <- mean(df_unique$caps_total, na.rm = TRUE)
caps_sd <- sd(df_unique$caps_total, na.rm = TRUE)

ggplot(caps_score_counts, aes(x = caps_total, y = count)) +
  geom_bar(stat = "identity", fill = "#6f4680", color = "black") +
  # Overlay the mean (dashed line)
  geom_vline(xintercept = caps_mean, linetype = "dashed", color = "black", size = 1) +
  # Overlay one standard deviation boundaries (dotted lines)
  geom_vline(xintercept = caps_mean - caps_sd, linetype = "dotted", color = "black", size = 1) +
  geom_vline(xintercept = caps_mean + caps_sd, linetype = "dotted", color = "black", size = 1) +
  # Overlay two standard deviation boundaries (dotted lines)
  geom_vline(xintercept = caps_mean - 2 * caps_sd, linetype = "dotted", color = "black", size = 1) +
  geom_vline(xintercept = caps_mean + 2 * caps_sd, linetype = "dotted", color = "black", size = 1) +
  labs(x = "Caps Total Score", y = "Number of Participants") +
  scale_y_continuous(breaks = seq(0, 50, 5), limits = c(0, 16)) +
  theme_minimal()


# Compute summary statistics for SPQ Total Score
spq_mean <- mean(df_unique$spq_total_score, na.rm = TRUE)
spq_sd <- sd(df_unique$spq_total_score, na.rm = TRUE)

ggplot(score_counts, aes(x = spq_total_score, y = count)) +
  geom_bar(stat = "identity", fill = "navy", color = "black") +
  # Overlay the mean (dashed line)
  geom_vline(xintercept = spq_mean, linetype = "dashed", color = "black", size = 1) +
  # Overlay one standard deviation boundaries (dotted lines)
  geom_vline(xintercept = spq_mean - spq_sd, linetype = "dotted", color = "black", size = 1) +
  geom_vline(xintercept = spq_mean + spq_sd, linetype = "dotted", color = "black", size = 1) +
  # Overlay two standard deviation boundaries (dotted lines)
  geom_vline(xintercept = spq_mean - 2 * spq_sd, linetype = "dotted", color = "black", size = 1) +
  geom_vline(xintercept = spq_mean + 2 * spq_sd, linetype = "dotted", color = "black", size = 1) +
  labs(x = "SPQ Total Score", y = "Number of Participants") +
  scale_y_continuous(breaks = seq(0, 25, 5), limits = c(0, 16)) +
  theme_minimal()

