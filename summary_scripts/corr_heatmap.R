# ======== 1.  LOAD DEPENDENCIES ========
library(Hmisc)   # rcorr()
library(dplyr)
library(tidyr)
library(ggplot2)

# ======== 2.  FUNCTION ========
cor_heatmap <- function(data,
                        drop_cols     = NULL,
                        digits_r      = 2,
                        show_triangle = c("both", "lower"),   # <- default is BOTH
                        fill_palette  = c("steelblue", "white", "firebrick"),
                        label_base    = 12,        # controls inside-tile text size
                        title         = "Balloon Task",
                        font_family   = "Courier New",
                        return_plot   = FALSE,
                        ...) {
  
  # ---- bookkeeping ----
  data_name     <- deparse(substitute(data))
  drop_cols     <- drop_cols %||% character(0)
  show_triangle <- match.arg(show_triangle)
  
  message("Data-frame   : ", data_name)
  message("Dropped cols : ", if (length(drop_cols)) paste(drop_cols, collapse = ", ") else "none")
  
  # ---- numeric columns only ----
  num_df <- data[ , !(names(data) %in% drop_cols)] |>
    dplyr::select(where(is.numeric))
  if (ncol(num_df) < 2L)
    stop("Need at least two numeric columns after dropping 'drop_cols'.")
  
  n_vars     <- ncol(num_df)
  label_size <- label_base / sqrt(n_vars)      # auto-scale inside-tile text
  
  # ---- r and p ----
  rc   <- Hmisc::rcorr(as.matrix(num_df), type = "pearson")
  rmat <- rc$r
  pmat <- rc$P
  
  # ---- long format for ggplot ----
  long <- as.data.frame(as.table(rmat)) |>
    rename(r = Freq) |>
    mutate(p     = as.vector(pmat),
           stars = case_when(
             p < .001 ~ "***",
             p < .01  ~ "**",
             p < .05  ~ "*",
             TRUE     ~ ""
           ),
           label = ifelse(
             stars == "",
             sprintf(paste0("%.", digits_r, "f"), r),           # just r
             sprintf("%s\n%.*f", stars, digits_r, r)            # stars above r
           ),
           Var1  = factor(Var1, levels = colnames(rmat)),
           Var2  = factor(Var2, levels = colnames(rmat)))
  
  # ---- drop nothing unless asked for lower only ----
  if (show_triangle == "lower") {
    long <- long |> filter(as.integer(Var1) > as.integer(Var2))
  }
  
  # ---- plot ----
  p <- ggplot(long, aes(Var2, Var1, fill = r)) +
    geom_tile(color = "grey90") +
    geom_text(aes(label = label),
              size     = label_size,
              family   = font_family,
              fontface = "bold") +
    scale_fill_gradient2(low       = fill_palette[1],
                         mid       = fill_palette[2],
                         high      = fill_palette[3],
                         midpoint  = 0,
                         limits    = c(-1, 1),
                         name      = "Pearson r") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_family = font_family) +
    theme(
      axis.text.x = element_text(angle  = 45,
                                 hjust  = 1,
                                 vjust  = 1,
                                 face   = "bold",
                                 size   = 11,
                                 margin = margin(t = -2)),
      axis.text.y = element_text(face   = "bold",
                                 size   = 11),
      axis.title  = element_blank(),
      panel.grid  = element_blank(),
      plot.title  = element_text(hjust = .5, face = "bold", size = 16)
    ) +
    coord_fixed(...)
  
  if (return_plot) return(p) else print(p)
}

# ======== 3.  BUILD + SAVE (whole matrix) ========
plt <- cor_heatmap(summary_dataframe,
                   drop_cols     = c("participant_id"),
                   show_triangle = "both",   # <- KEEP BOTH HALVES
                   return_plot   = TRUE)

ggsave("~/Desktop/balloon_task_corr.png",
       plot   = plt,
       width  = 14,
       height = 12,
       units  = "in",
       dpi    = 600,
       bg     = "white")

