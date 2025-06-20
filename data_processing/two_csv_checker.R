# Read the CSV files
df1 <- read.csv("/Users/zachtefertiller/Desktop/bart_rl_180_data/balloon_task_clean_data.csv", stringsAsFactors = FALSE)
df2 <- read.csv("/Users/zachtefertiller/Desktop/balloon_task_clean_data.csv", stringsAsFactors = FALSE)

# Create a composite key for each data frame using trial_number and participant_id.
df1$composite_key <- paste(df1$trial_number, df1$participant_id, sep = "_")
df2$composite_key <- paste(df2$trial_number, df2$participant_id, sep = "_")

# Set the composite key as row names.
rownames(df1) <- df1$composite_key
rownames(df2) <- df2$composite_key

# Identify common columns (excluding the composite key if present)
common_cols <- intersect(setdiff(names(df1), "composite_key"), setdiff(names(df2), "composite_key"))

# Identify common rows based on the composite row names.
common_rows <- intersect(rownames(df1), rownames(df2))

# Subset each data frame to only the common rows and common columns.
df1_common <- df1[common_rows, common_cols, drop = FALSE]
df2_common <- df2[common_rows, common_cols, drop = FALSE]

# Compare the common parts for exact equality.
if (isTRUE(all.equal(df1_common, df2_common))) {
  print("All common rows and columns match exactly!")
} else {
  print("Differences found in the common rows/columns!")
  # Optionally, you can inspect the differences:
  print(all.equal(df1_common, df2_common))
}
