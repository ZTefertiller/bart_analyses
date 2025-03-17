library(dplyr)
library(rjson)

# read in file
summary_data <- readLines('/Users/zachtefertiller/Desktop/jatos_results_data_20250116122322.json')

# parse line by line
parsed_data <- lapply(summary_data, fromJSON)

# turn into df
summary_df <- do.call(rbind, lapply(parsed_data, as.data.frame))
print(summary_df)

#average bonus
bonuses <- as.numeric(summary_df$bonus_payment)
average_bonus <- mean(bonuses)
print(average_bonus)

# match prolific ID to bonus payment
payments <- paste0(summary_df$prolific_id, ",", summary_df$bonus_payment)

# prompt for unique file name ending 
user_input <- readline("Enter a suffix (or leave blank for default): ")

# concatenate file name
if (nzchar(user_input)) {
  # nzchar() checks if it's a non-empty string
  output_file <- paste0("bonus_payments_", user_input, ".txt")
} else {
  output_file <- "bonus_payments_14_01_2025.txt"
}

writeLines(payments, output_file)
cat("Wrote data to", output_file, "\n")
