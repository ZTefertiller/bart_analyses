library(jsonlite)
library(dplyr)
library(readr)
library(stringr)
library(tibble)
library(purrr)
library(tidyr)

# read in jatos data, but not in typical json so convert to 
# "ndjson" (newline separated instead of comma separated)
raw_text <- read_file("/Users/zachtefertiller/Desktop/bart_only.json")
ndjson_text <- str_replace_all(raw_text, "\\}\\{", "\\}\n\\{")
write_file(ndjson_text, "/Users/zachtefertiller/Desktop/bart.ndjson")
# read in line by line (each line is its own json)
lines <- readLines("/Users/zachtefertiller/Desktop/bart.ndjson")

# turn each line into list
parsed <- lapply(seq_along(lines), function(i) {
  line <- lines[i]
  parsed_line <- tryCatch({
    fromJSON(line)
  }, error = function(e) {
    warning(paste("Line", i, "failed to parse:", e$message))
    NULL
  })
  return(parsed_line)
})

parsed <- Filter(Negate(is.null), parsed)
is_trial <- function(x) "trial_number" %in% names(x)
parsed_trials <- parsed[sapply(parsed, is_trial)]   # trials
parsed_meta   <- parsed[!sapply(parsed, is_trial)]  # non‐trials (questionnaire data)

# define the required fields for trials
required_fields <- c(
  "participant_id",
  "trial_number",
  "balloon_color",
  "inflations",
  "optimal_inflations",
  "popped",
  "points_earned",
  "inactivity_pop",
  "total_points_so_far",
  "average_inflation_rt",
  "inflation_rts"
)

df_trials <- bind_rows(lapply(seq_along(parsed_trials), function(i) {
  x <- parsed_trials[[i]]
  
  trial_data <- setNames(as.list(rep(NA, length(required_fields))), required_fields)
  
  # populate trial_data with available fields from x
  for (field in names(x)) {
    if (field == "inflation_rts") {
      trial_data[[field]] <- list(x[[field]])
    } else if (field == "average_inflation_rt") {
      trial_data[[field]] <- if (!is.null(x[[field]])) as.numeric(x[[field]]) else NA_real_
    } else {
      trial_data[[field]] <- x[[field]]
    }
  }
  
  as_tibble(trial_data)
}))

df_meta <- bind_rows(lapply(parsed_meta, function(x) {
  as_tibble(x)
}))

df_final <- left_join(df_trials, df_meta, by = "participant_id")
write_csv(df_final, "/Users/zachtefertiller/Desktop/balloon_task_all_data.csv")
glimpse(df_final)

# Process questionnaire data
json_file <- "/Users/zachtefertiller/Desktop/jatos_results_all.json"
lines <- readLines(json_file)
filtered_lines <- lines[!grepl('"prolific_id"', lines) & !grepl('trial_number', lines)]
filtered_json_file <- "/Users/zachtefertiller/Desktop/jatos_results_filtered.json"
writeLines(filtered_lines, filtered_json_file)
raw_text <- read_file(filtered_json_file)
ndjson_text <- str_replace_all(raw_text, "\\}\\{", "\\}\n\\{")
write_file(ndjson_text, "/Users/zachtefertiller/Desktop/jatos_results_filtered.ndjson")
df <- stream_in(file("/Users/zachtefertiller/Desktop/jatos_results_filtered.ndjson"))

# Extract SPQ data properly
df_spq_only <- df %>%
  filter(!map_lgl(spq_data, is.null)) %>%
  group_by(participant_id) %>%
  slice(1) %>%
  distinct(participant_id, spq_data) %>%
  ungroup()

df_spq <- df_spq_only %>%
  rowwise() %>%
  mutate(spq_diag = list(diag(as.matrix(spq_data)))) %>%
  ungroup() %>%
  unnest_wider(col = spq_diag, names_sep = "_") %>%
  rename_with(
    .cols = starts_with("spq_diag_"),
    .fn   = ~ paste0("spq_q", seq_along(.))
  )

# ─────────────────────────────────────────────────────────────────────────
# Extract PDI individual question answers (one row per participant,
# pdi_q1, pdi_q2, … as numeric 1/0)
# ─────────────────────────────────────────────────────────────────────────
pdi_cols <- names(df)[grepl("^pdi_q\\d+$", names(df))]

if (length(pdi_cols) > 0) {
  df_pdi_answers <- df %>%
    rowwise() %>%
    # keep only rows where at least one pdi_q* is non‐NULL
    filter(any(c_across(all_of(pdi_cols)) %>% map_lgl(~ !is.null(.x)))) %>%
    # extract "answer" from each nested pdi_q* list‐column,
    # then convert "Yes"/"No" → 1/0 (NA if neither)
    mutate(across(
      all_of(pdi_cols),
      ~ {
        if (is.list(.x) && "answer" %in% names(.x)) {
          ans <- .x$answer
          case_when(
            ans == "Yes" ~ 1L,
            ans == "No"  ~ 0L,
            TRUE         ~ NA_integer_
          )
        } else {
          NA_integer_
        }
      }
    )) %>%
    ungroup() %>%
    # collapse to one row per participant_id, taking the first non‐NA per question
    group_by(participant_id) %>%
    summarise(
      across(
        all_of(pdi_cols),
        ~ {
          tmp <- na.omit(.x)
          if (length(tmp)) tmp[1] else NA_integer_
        }
      ),
      .groups = "drop"
    )
}

# ─────────────────────────────────────────────────────────────────────────
# Extract CAPS individual question answers (one row per participant,
# caps_q1, caps_q2, … as numeric 1/0)
# ─────────────────────────────────────────────────────────────────────────
caps_cols <- names(df)[grepl("^caps_q\\d+$", names(df))]

if (length(caps_cols) > 0) {
  df_caps_answers <- df %>%
    rowwise() %>%
    # keep only rows where at least one caps_q* is non‐NULL
    filter(any(c_across(all_of(caps_cols)) %>% map_lgl(~ !is.null(.x)))) %>%
    # extract "answer" from each nested caps_q* list‐column,
    # then convert "Yes"/"No" → 1/0 (NA if neither)
    mutate(across(
      all_of(caps_cols),
      ~ {
        if (is.list(.x) && "answer" %in% names(.x)) {
          ans <- .x$answer
          case_when(
            ans == "Yes" ~ 1L,
            ans == "No"  ~ 0L,
            TRUE         ~ NA_integer_
          )
        } else {
          NA_integer_
        }
      }
    )) %>%
    ungroup() %>%
    # collapse to one row per participant_id, taking the first non‐NA per question
    group_by(participant_id) %>%
    summarise(
      across(
        all_of(caps_cols),
        ~ {
          tmp <- na.omit(.x)
          if (length(tmp)) tmp[1] else NA_integer_
        }
      ),
      .groups = "drop"
    )
}

# ─────────────────────────────────────────────────────────────────────────
# Join SPQ, PDI, and CAPS into a single “questionnaire” table
# ─────────────────────────────────────────────────────────────────────────
if (exists("df_spq")) {
  df_questionnaire_all <- df_spq
  if (exists("df_pdi_answers")) {
    df_questionnaire_all <- left_join(df_questionnaire_all, df_pdi_answers, by = "participant_id")
  }
  if (exists("df_caps_answers")) {
    df_questionnaire_all <- left_join(df_questionnaire_all, df_caps_answers, by = "participant_id")
  }
} else {
  df_questionnaire_all <- tibble(participant_id = character())
  if (exists("df_pdi_answers"))  df_questionnaire_all <- left_join(df_questionnaire_all, df_pdi_answers, by = "participant_id")
  if (exists("df_caps_answers")) df_questionnaire_all <- left_join(df_questionnaire_all, df_caps_answers, by = "participant_id")
}

# Remove the original nested pdi_q* and caps_q* columns from df
df_reduced <- df %>%
  select(-all_of(pdi_cols), -all_of(caps_cols))

# Join questionnaire data back to the reduced df (nested columns dropped)
df_joined <- left_join(df_reduced, df_questionnaire_all, by = "participant_id")

# Create wide format, removing any remaining nested columns and duplicates
df_wide <-
  jsonlite::flatten(df_joined) %>%
  group_by(participant_id) %>%
  summarise(
    across(
      everything(),
      ~ {
        non_na_vals <- .x[!is.na(.x) & .x != "NULL"]
        if (length(non_na_vals) == 0) {
          NA_character_
        } else {
          paste(unique(non_na_vals), collapse = ", ")
        }
      },
      .names = "{col}"
    )
  )

df_final <- left_join(df_final, df_wide, by = "participant_id")

# ─────────────────────────────────────────────────────────────────────────
# Sort df_final by participant_id, then by trial_number
# ─────────────────────────────────────────────────────────────────────────
df_final <- df_final %>%
  arrange(participant_id, trial_number)

# Remove problematic columns and duplicates
df_final <- subset(df_final, select = -c(
  # Remove duplicate PPGM columns
  ppgm1a, ppgm1b, ppgm2, ppgm3a, ppgm3b, ppgm4, ppgm5a, ppgm5b, ppgm6,
  ppgm7, ppgm8, ppgm9, ppgm10a, ppgm11, ppgm12, ppgm13, ppgm14,
  # Remove timestamp columns
  timestamp_start, timestamp_end, start_timestamp, end_timestamp,
  # Remove problematic SPQ data columns
  spq_data.x, spq_data.y,
  # Remove event column
  event
))

# Clean up list columns
df_final$mdq_q1 <- gsub("^c\\(|\\)$", "", df_final$mdq_q1)
df_final$mdq_q1 <- gsub("\"", "", df_final$mdq_q1)
df_final$gambling_types <- gsub("^c\\(|\\)$", "", df_final$gambling_types)
df_final$gambling_types <- gsub("\"", "", df_final$gambling_types)

# Reorder columns: required fields first, then alphabetized
all_columns <- colnames(df_final)
alphabetized_cols <- setdiff(all_columns, required_fields)
alphabetized_cols <- sort(alphabetized_cols)
df_final <- df_final %>%
  select(all_of(required_fields), all_of(alphabetized_cols))

# Save raw data (converting lists to strings for CSV compatibility)
write.csv(
  df_final %>% mutate(across(where(is.list), ~ sapply(., toString))),
  "/Users/zachtefertiller/Desktop/balloon_task_all_data.csv",
  row.names = FALSE
)

# ─────────────────────────────────────────────────────────────────────────
# Create clean dataset with only essential columns
# ─────────────────────────────────────────────────────────────────────────
df_clean <- df_final %>%
  dplyr::select(
    # Trial data
    participant_id, trial_number, balloon_color, inflations,
    optimal_inflations, popped, points_earned, inactivity_pop,
    total_points_so_far, average_inflation_rt, inflation_rts,
    # Summary scores
    total_score, total_pops, total_inactive, total_money, task_time,
    # Questionnaire totals
    spq_total_score, pdi_total, pdi_distress, pdi_frequency, pdi_conviction,
    caps_total, caps_distress, caps_intrusiveness, caps_frequency,
    phq_total, ipip_total_score, ppgm_total, mdq_total,
    # Individual questionnaire questions (numeric 1/0)
    starts_with("spq_q"), starts_with("pdi_q"), starts_with("caps_q"),
    # Demographics and other variables
    antipsychotics, family_bipolar, family_schizophrenia,
    drink_frequency, cigarettes_count, glasses_per_day
  ) %>%
  # Sort df_clean by participant_id, then by trial_number
  arrange(participant_id, trial_number)

# Save clean dataset
write.csv(
  df_clean %>% mutate(across(where(is.list), ~ sapply(., toString))),
  "/Users/zachtefertiller/Desktop/balloon_task_clean_data.csv",
  row.names = FALSE
)

# Print summary of what was processed
cat("Processing complete!\n")
cat("Total participants:", length(unique(df_final$participant_id)), "\n")
cat("Total trials:", nrow(df_final), "\n")
cat("SPQ questions extracted:", sum(grepl("^spq_q\\d+$", names(df_final))), "\n")
cat("PDI questions extracted:", sum(grepl("^pdi_q\\d+$", names(df_final))), "\n")
cat("CAPS questions extracted:", sum(grepl("^caps_q\\d+$", names(df_final))), "\n")
cat("Columns in final dataset:", ncol(df_final), "\n")
cat("Columns in clean dataset:", ncol(df_clean), "\n")

