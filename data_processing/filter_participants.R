## we cleaned the data but lets exclude participants who had technical issues ie not enough trials
library(jsonlite)
library(dplyr)
library(readr)
library(stringr)
library(tibble)

df <- read.csv("/Users/zachtefertiller/Desktop/balloon_task_clean_data.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

df %>%
  group_by(participant_id) %>%
  summarise(num_trials = n_distinct(trial_number)) %>%
  arrange(num_trials)


df_filtered <- df %>%
  group_by(participant_id) %>%
  filter(
    n_distinct(trial_number) == 180,                      
    all(sort(unique(trial_number)) == 1:180)              
  ) %>%
  ungroup()             





trial_counts <- df %>%
  group_by(participant_id) %>%
  summarise(num_trials = n_distinct(trial_number)) %>%
  arrange(num_trials)
print(trial_counts)


df_filtered <- df %>%
  group_by(participant_id) %>%
  filter(
    n_distinct(trial_number) == 180,                     
    all(sort(unique(trial_number)) == 1:180)             
  ) %>%
  ungroup()


complete_participants <- trial_counts %>%
  filter(num_trials == 180) %>%
  pull(participant_id)

df_filtered <- df %>%
  filter(participant_id %in% complete_participants)

## remove rows without a trial number or participant_id value
df_clean <- df_filtered %>%
  filter(
    !is.na(participant_id),
    !is.na(trial_number)
  )

remaining_duplicates <- df_clean %>%
  group_by(participant_id, trial_number) %>%
  filter(n() > 1) %>%
  ungroup()

if(nrow(remaining_duplicates) > 0){
  cat("There are still duplicate trials:\n")
  print(remaining_duplicates)
} else {
  cat("All duplicate trials have been successfully removed.\n")
}

# duplicate trials
duplicate_trials <- df_clean %>%
  group_by(participant_id, trial_number) %>%
  filter(n() > 1) %>%
  arrange(participant_id, trial_number) %>%
  ungroup()

if(nrow(duplicate_trials) > 0){
  cat("Duplicate trials found:\n")
  print(duplicate_trials)
} else {
  cat("No duplicate trials found.\n")
}


df_clean_unique <- df_clean %>%
  distinct(participant_id, trial_number, .keep_all = TRUE)


df_clean_unique <- df_clean %>%
  group_by(participant_id, trial_number) %>%
  slice(1) %>%  # Keeps the first occurrence in each group
  ungroup()






# attention check exclusion -- when programming questionnaires I should have
# been a lot more careful and kept the naming / boolean conventions the same
# i hope this does not make you go insane

# ac rejection criteria:
# if value > 1 caps_attention_fails, pdi_attention_fails, phq_attention_fails, spq_attention_fails
# if FALSE ipip_attention_passed
# if "attention" not in mdq_q1 if mdq_1 is taken as a string or list of strings 

ac_exclusion <- TRUE

if(ac_exclusion == TRUE) {
  cat("Pre-exclusion number of participants:\n")
  print(length(unique(df_clean_unique$participant_id)))
  df_clean_unique <- df_clean_unique %>%
    group_by(participant_id) %>%
    filter(!(caps_attention_fails > 0 | spq_attention_fails > 0 | pdi_attention_fails > 0 | 
               phq_attention_fails > 0 | ipip_attention_passed == FALSE | !str_detect(mdq_q1, regex("attention", ignore_case = TRUE))))
  cat("Post-exclusion number of participants:\n")
  print(length(unique(df_clean_unique$participant_id)))
}  




df_clean_unique <- df_clean_unique %>%
  group_by(participant_id) %>%
  mutate(adjusted_inflations = mean(inflations[popped == 0 & balloon_color == "b"], na.rm = TRUE)) %>%
  ungroup()



write.csv(df_clean_unique %>% mutate(across(where(is.list), ~ sapply(., toString))), "/Users/zachtefertiller/Desktop/balloon_task_clean_data.csv", row.names = FALSE)


## or down here more strict if we need to remove the people that had duplicate trials.


