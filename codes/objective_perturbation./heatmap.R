# ==============================================================================
# Script: heatmap.R (exp3)
# Description: Generates publication-ready heatmaps of empirical coverage and CI 
#              widths across sample sizes n and privacy parameters epsilon for 
#              the logistic regression objective perturbation experiments.
#
# Inputs:
#   summary/Repro_logistic.csv from each configuration folder
# Outputs:
#   heatmap_coverage.pdf, heatmap_width.pdf, figure6_operatorII.pdf
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Path Initialization
# ------------------------------------------------------------------------------

ROOT_DIR <- path.expand("~/R_Simuls/exp3_objectiveperturb/eff")
setwd(ROOT_DIR)

# ------------------------------------------------------------------------------
# 2. Grid Dimensions & Label Formatter
# ------------------------------------------------------------------------------

n_values   <- c(100, 200, 500, 1000)
eps_values <- c(0.1, 0.3, 1, 3)

format_eps <- function(eps) {
  if (eps == 0.1) return("0.1")
  if (eps == 0.3) return("0.3")
  if (eps == 1.0) return("1")
  if (eps == 3.0) return("3")
  format(eps, scientific = FALSE, trim = TRUE)
}

# ------------------------------------------------------------------------------
# 3. Read Summary Data Files
# ------------------------------------------------------------------------------

results <- data.frame(
  n                   = integer(0),
  epsilon             = numeric(0),
  coverage            = numeric(0),
  width_all           = numeric(0),
  width_real          = numeric(0),
  n_full_range        = integer(0),
  n_no_accepted_point = integer(0),
  mean_seconds        = numeric(0)
)

for (n_cur in n_values) {
  for (eps_cur in eps_values) {
    eps_string  <- format_eps(eps_cur)
    folder_name <- sprintf("n_%d_eps_%s_eff", n_cur, eps_string)
    csv_file    <- file.path(ROOT_DIR, folder_name, "summary", "Repro_logistic.csv")
    
    cat(sprintf("Reading: %s\n", csv_file))
    
    if (!file.exists(csv_file)) {
      warning(sprintf("FILE NOT FOUND: %s", csv_file))
      next
    }
    
    tmp <- read.csv(csv_file, header = TRUE)
    required_columns <- c(
      "coverage", "width_all", "width_real", "n_full_range",
      "n_no_accepted_point", "mean_seconds"
    )
    
    missing_columns <- setdiff(required_columns, colnames(tmp))
    if (length(missing_columns) > 0L) {
      stop(
        sprintf(
          "\nMissing column(s) in:\n%s\nMissing: %s",
          csv_file, paste(missing_columns, collapse = ", ")
        )
      )
    }
    
    results <- rbind(
      results,
      data.frame(
        n                   = n_cur,
        epsilon             = eps_cur,
        coverage            = tmp$coverage[1],
        width_all           = tmp$width_all[1],
        width_real          = tmp$width_real[1],
        n_full_range        = tmp$n_full_range[1],
        n_no_accepted_point = tmp$n_no_accepted_point[1],
        mean_seconds        = tmp$mean_seconds[1]
      )
    )
  }
}

cat("\n============================================================\n")
cat("RESULTS READ FROM CSV FILES\n")
cat("============================================================\n")
print(results, row.names = FALSE)
cat("============================================================\n\n")

if (nrow(results) != 16L) {
  warning(sprintf("Expected 16 configurations but found %d.", nrow(results)))
}

write.csv(results, file.path(ROOT_DIR, "heatmap_data.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Matrix Structuring
# ------------------------------------------------------------------------------

coverage_matrix <- matrix(
  NA_real_,
  nrow = length(eps_values),
  ncol = length(n_values),
  dimnames = list(paste0("\u03b5 = ", eps_values), n_values)
)

width_matrix      <- coverage_matrix
width_real_matrix <- coverage_matrix
fallback_matrix   <- coverage_matrix

for (i in seq_along(eps_values)) {
  for (j in seq_along(n_values)) {
    row_match <- which(
      results$n == n_values[j] & abs(results$epsilon - eps_values[i]) < 1e-12
    )
    if (length(row_match) == 1L) {
      coverage_matrix[i, j]   <- results$coverage[row_match]
      width_matrix[i, j]      <- results$width_all[row_match]
      width_real_matrix[i, j] <- results$width_real[row_match]
      fallback_matrix[i, j]   <- results$n_full_range[row_match]
    }
  }
}

cat("\nCOVERAGE MATRIX\n")
print(round(coverage_matrix, 3))

cat("\nWIDTH MATRIX -- ALL INTERVALS\n")
print(round(width_matrix, 2))

# ------------------------------------------------------------------------------
# 5. Heatmap Plotting Function
# ------------------------------------------------------------------------------

draw_heatmap <- function(
    Z,
    main_title,
    legend_title,
    digits = 2,
    palette_function
) {
  nr <- nrow(Z)
  nc <- ncol(Z)
  
  # Orient so epsilon = 0.1 is at the top
  Z_plot <- Z[nr:1, , drop = FALSE]
  
  image(
    x = seq_len(nc),
    y = seq_len(nr),
    z = t(Z_plot),
    col = palette_function(100),
    axes = FALSE,
    xlab = "Number of Samples (N)",
    ylab = expression(paste("Privacy Parameter (", epsilon, ")")),
    main = main_title
  )
  
  axis(side = 1, at = seq_len(nc), labels = n_values)
  axis(side = 2, at = seq_len(nr), labels = rev(eps_values), las = 1)
  
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      original_i <- nr - i + 1L
      value <- Z[original_i, j]
      
      if (is.finite(value)) {
        label <- formatC(value, format = "f", digits = digits)
        text(x = j, y = i, labels = label, cex = 1.15)
      }
    }
  }
  
  box()
  
  # Color scale bar
  usr     <- par("usr")
  x_left  <- usr[2] + 0.15
  x_right <- usr[2] + 0.30
  y_seq   <- seq(usr[3], usr[4], length.out = 101)
  z_range <- range(Z, finite = TRUE)
  cols    <- palette_function(100)
  
  for (k in 1:100) {
    rect(
      xleft = x_left, ybottom = y_seq[k], xright = x_right, ytop = y_seq[k + 1L],
      col = cols[k], border = NA, xpd = TRUE
    )
  }
  
  axis(
    side = 4,
    at = seq(usr[3], usr[4], length.out = 5),
    labels = formatC(seq(z_range[1], z_range[2], length.out = 5), format = "f", digits = digits),
    las = 1,
    line = -1,
    xpd = TRUE
  )
  
  mtext(legend_title, side = 4, line = 3.3, xpd = TRUE)
}

# ------------------------------------------------------------------------------
# 6. Color Palettes
# ------------------------------------------------------------------------------

coverage_palette <- function(n) {
  colorRampPalette(c("#fff5f0", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d"))(n)
}

width_palette <- function(n) {
  colorRampPalette(c("#fff5f0", "#fcbba1", "#fb6a4a", "#cb181d", "#67000d"))(n)
}

# ------------------------------------------------------------------------------
# 7. Render & Export Figures
# ------------------------------------------------------------------------------

# Coverage Heatmap
pdf(file.path(ROOT_DIR, "heatmap_coverage.pdf"), width = 7, height = 5.5)
par(mar = c(5, 5, 4, 7))
draw_heatmap(
  Z = coverage_matrix, main_title = "Coverage of 90% CI\nEfficient Operator II",
  legend_title = "Coverage Rate", digits = 3, palette_function = coverage_palette
)
dev.off()

png(file.path(ROOT_DIR, "heatmap_coverage.png"), width = 1600, height = 1200, res = 200)
par(mar = c(5, 5, 4, 7))
draw_heatmap(
  Z = coverage_matrix, main_title = "Coverage of 90% CI\nEfficient Operator II",
  legend_title = "Coverage Rate", digits = 3, palette_function = coverage_palette
)
dev.off()

# Interval Width Heatmap
pdf(file.path(ROOT_DIR, "heatmap_width.pdf"), width = 7, height = 5.5)
par(mar = c(5, 5, 4, 7))
draw_heatmap(
  Z = width_matrix, main_title = "Width of 90% CI\nEfficient Operator II",
  legend_title = "CI Width", digits = 2, palette_function = width_palette
)
dev.off()

png(file.path(ROOT_DIR, "heatmap_width.png"), width = 1600, height = 1200, res = 200)
par(mar = c(5, 5, 4, 7))
draw_heatmap(
  Z = width_matrix, main_title = "Width of 90% CI\nEfficient Operator II",
  legend_title = "CI Width", digits = 2, palette_function = width_palette
)
dev.off()

# Combined Figure
pdf(file.path(ROOT_DIR, "figure6_operatorII.pdf"), width = 12, height = 5.5)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 6))
draw_heatmap(
  Z = coverage_matrix, main_title = "Coverage of 90% CI\n(Efficient Operator II)",
  legend_title = "Coverage Rate", digits = 3, palette_function = coverage_palette
)
draw_heatmap(
  Z = width_matrix, main_title = "Width of 90% CI\n(Efficient Operator II)",
  legend_title = "CI Width", digits = 2, palette_function = width_palette
)
dev.off()

png(file.path(ROOT_DIR, "figure6_operatorII.png"), width = 2400, height = 1000, res = 200)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 6))
draw_heatmap(
  Z = coverage_matrix, main_title = "Coverage of 90% CI\n(Efficient Operator II)",
  legend_title = "Coverage Rate", digits = 3, palette_function = coverage_palette
)
draw_heatmap(
  Z = width_matrix, main_title = "Width of 90% CI\n(Efficient Operator II)",
  legend_title = "CI Width", digits = 2, palette_function = width_palette
)
dev.off()

cat("\n============================================================\n")
cat("HEATMAPS CREATED SUCCESSFULLY\n")
cat("============================================================\n")
cat(sprintf("Coverage PDF : %s\n", file.path(ROOT_DIR, "heatmap_coverage.pdf")))
cat(sprintf("Width PDF    : %s\n", file.path(ROOT_DIR, "heatmap_width.pdf")))
cat(sprintf("Combined PDF : %s\n", file.path(ROOT_DIR, "figure6_operatorII.pdf")))
cat("============================================================\n")
