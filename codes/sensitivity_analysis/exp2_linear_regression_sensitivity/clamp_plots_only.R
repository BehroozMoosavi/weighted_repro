# ==============================================================================
# Script: clamp_plots_only.R
# Description: Generates publication-ready heatmaps of empirical power 
#              (beta1 = 1) and Type-I error (beta1 = 0) across clamping thresholds 
#              Delta and sample sizes n using precomputed results.
#
# Inputs:
#   exp2_delta_sensitivity_long.csv
# Outputs:
#   exp2_delta_sensitivity_Figure6_clean.pdf
#   exp2_delta_sensitivity_Figure6_clean.png
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Environment Setup & File Loading
# ------------------------------------------------------------------------------

graphics.off()

PROJECT_DIR  <- path.expand("~/R_Simuls/sensitivity_analysis/exp2")
results_file <- file.path(PROJECT_DIR, "exp2_delta_sensitivity_long.csv")

if (!file.exists(results_file)) {
  stop(paste0("Cannot find the completed CSV:\n", results_file))
}

results <- read.csv(results_file, stringsAsFactors = FALSE)
cat("\nLoaded:\n", results_file, "\nRows: ", nrow(results), "\n\n")

# ------------------------------------------------------------------------------
# 1. Sensitivity Grid Specifications
# ------------------------------------------------------------------------------

n_list     <- c(100L, 200L, 500L, 1000L)
Delta_list <- c(0.5, 1.0, 2.0, 5.0, 10.0)

required_cols <- c("beta1_true", "Delta", "n", "rejection_probability")
missing_cols <- setdiff(required_cols, names(results))
if (length(missing_cols) > 0L) {
  stop(paste("Missing column(s):", paste(missing_cols, collapse = ", ")))
}

# ------------------------------------------------------------------------------
# 2. Matrix Structuring
# ------------------------------------------------------------------------------

make_result_matrix <- function(beta_value) {
  mat <- matrix(
    NA_real_,
    nrow = length(Delta_list),
    ncol = length(n_list),
    dimnames = list(as.character(Delta_list), as.character(n_list))
  )
  
  for (i in seq_along(Delta_list)) {
    for (j in seq_along(n_list)) {
      z <- results[
        abs(results$beta1_true - beta_value) < 1e-12 &
          abs(results$Delta - Delta_list[i]) < 1e-12 &
          results$n == n_list[j],
        ,
        drop = FALSE
      ]
      
      if (nrow(z) == 1L) {
        mat[i, j] <- z[["rejection_probability"]]
      } else if (nrow(z) == 0L) {
        warning(sprintf("Missing cell: beta1=%.1f, Delta=%.1f, n=%d", beta_value, Delta_list[i], n_list[j]))
      } else {
        stop(sprintf("Duplicate cell: beta1=%.1f, Delta=%.1f, n=%d", beta_value, Delta_list[i], n_list[j]))
      }
    }
  }
  mat
}

power_mat <- make_result_matrix(beta_value = 1)
size_mat  <- make_result_matrix(beta_value = 0)

cat("\n============================================================\n")
cat("POWER -- beta1 = 1\nRows = Delta, columns = n\n")
cat("============================================================\n")
print(round(power_mat, 3))

cat("\n============================================================\n")
cat("TYPE-I ERROR -- beta1 = 0\nRows = Delta, columns = n\n")
cat("============================================================\n")
print(round(size_mat, 3))

if (anyNA(power_mat)) warning("Power matrix contains missing values.")
if (anyNA(size_mat))  warning("Type-I error matrix contains missing values.")

# ------------------------------------------------------------------------------
# 3. Palette & Display Functions
# ------------------------------------------------------------------------------

heat_cols <- grDevices::colorRampPalette(
  c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08519C", "#08306B")
)(256)

prob_color <- function(x) {
  if (!is.finite(x)) return("grey90")
  x <- min(1, max(0, x))
  idx <- 1L + floor(x * 255)
  idx <- max(1L, min(256L, idx))
  heat_cols[idx]
}

prob_text_color <- function(x) {
  if (!is.finite(x)) return("black")
  if (x >= 0.55) "white" else "black"
}

prob_label <- function(x) {
  if (!is.finite(x)) return("--")
  sprintf("%.3f", x)
}

# ------------------------------------------------------------------------------
# 4. Canvas Drawing Routines (Base Graphics)
# ------------------------------------------------------------------------------

draw_heatmap_canvas <- function(
    mat,
    x_left,
    x_right,
    y_bottom,
    y_top,
    panel_label,
    title_expr,
    show_y_labels = TRUE
) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  cell_w <- (x_right - x_left) / nc
  cell_h <- (y_top - y_bottom) / nr
  
  # Orient Delta = 0.5 at top, Delta = 10 at bottom
  for (r_display in seq_len(nr)) {
    original_row <- r_display
    y1 <- y_top - r_display * cell_h
    y2 <- y_top - (r_display - 1) * cell_h
    
    for (j in seq_len(nc)) {
      x1 <- x_left + (j - 1) * cell_w
      x2 <- x_left + j * cell_w
      val <- mat[original_row, j]
      
      rect(xleft = x1, ybottom = y1, xright = x2, ytop = y2,
           col = prob_color(val), border = "white", lwd = 1.7)
      text(x = 0.5 * (x1 + x2), y = 0.5 * (y1 + y2), labels = prob_label(val),
           col = prob_text_color(val), cex = 0.92, font = 2, family = "serif")
    }
  }
  
  # Outer boundary
  rect(xleft = x_left, ybottom = y_bottom, xright = x_right, ytop = y_top,
       border = "black", lwd = 0.9)
  
  # X-axis ticks & labels
  for (j in seq_len(nc)) {
    x_center <- x_left + (j - 0.5) * cell_w
    text(x = x_center, y = y_bottom - 0.030, labels = n_list[j], cex = 0.88, family = "serif")
  }
  text(x = 0.5 * (x_left + x_right), y = y_bottom - 0.080,
       labels = expression(paste("Sample size, ", n)), cex = 0.98, family = "serif")
  
  # Y-axis labels & title
  if (show_y_labels) {
    for (i in seq_len(nr)) {
      y_center <- y_top - (i - 0.5) * cell_h
      text(x = x_left - 0.025, y = y_center, labels = Delta_list[i],
           adj = 1, cex = 0.88, family = "serif")
    }
    text(x = x_left - 0.092, y = 0.5 * (y_bottom + y_top),
         labels = expression(paste("Clamping threshold, ", Delta)),
         srt = 90, cex = 0.98, family = "serif")
  }
  
  # Panel labels
  text(x = x_left, y = y_top + 0.060, labels = panel_label,
       adj = c(0, 0.5), font = 2, cex = 1.02, family = "serif")
  text(x = 0.5 * (x_left + x_right), y = y_top + 0.060, labels = title_expr,
       font = 2, cex = 1.02, family = "serif")
}

draw_colorbar_canvas <- function(
    x_left,
    x_right,
    y_bottom,
    y_top
) {
  n_segments <- 256L
  segment_h <- (y_top - y_bottom) / n_segments
  
  for (k in seq_len(n_segments)) {
    y1 <- y_bottom + (k - 1) * segment_h
    y2 <- y_bottom + k * segment_h
    val <- (k - 1) / (n_segments - 1)
    rect(xleft = x_left, ybottom = y1, xright = x_right, ytop = y2,
         col = prob_color(val), border = NA)
  }
  rect(xleft = x_left, ybottom = y_bottom, xright = x_right, ytop = y_top,
       border = "black", lwd = 0.8)
  
  ticks <- c(0, 0.25, 0.50, 0.75, 1)
  for (tick in ticks) {
    yy <- y_bottom + tick * (y_top - y_bottom)
    segments(x0 = x_right, y0 = yy, x1 = x_right + 0.008, y1 = yy, lwd = 0.8)
    text(x = x_right + 0.014, y = yy,
         labels = if (tick == 0) "0" else if (tick == 1) "1" else sprintf("%.2f", tick),
         adj = 0, cex = 0.80, family = "serif")
  }
  text(x = x_right + 0.075, y = 0.5 * (y_bottom + y_top),
       labels = "Rejection probability", srt = 90, cex = 0.90, family = "serif")
}

draw_complete_figure <- function() {
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0), xaxs = "i", yaxs = "i", family = "serif")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  
  y_bottom <- 0.18
  y_top    <- 0.82
  left_A   <- 0.10
  right_A  <- 0.45
  left_B   <- 0.51
  right_B  <- 0.86
  cb_left  <- 0.895
  cb_right <- 0.915
  
  # Panel A: Power
  draw_heatmap_canvas(
    mat = power_mat, x_left = left_A, x_right = right_A,
    y_bottom = y_bottom, y_top = y_top, panel_label = "(a)",
    title_expr = expression(paste("Power  (", beta[1], " = 1)")),
    show_y_labels = TRUE
  )
  
  # Panel B: Type-I Error
  draw_heatmap_canvas(
    mat = size_mat, x_left = left_B, x_right = right_B,
    y_bottom = y_bottom, y_top = y_top, panel_label = "(b)",
    title_expr = expression(paste("Type-I error  (", beta[1], " = 0)")),
    show_y_labels = FALSE
  )
  
  # Shared Color Scale Bar
  draw_colorbar_canvas(x_left = cb_left, x_right = cb_right, y_bottom = y_bottom, y_top = y_top)
}

# ------------------------------------------------------------------------------
# 5. Output Graphics Export
# ------------------------------------------------------------------------------

pdf_file <- file.path(PROJECT_DIR, "exp2_delta_sensitivity_Figure6_clean.pdf")
grDevices::pdf(file = pdf_file, width = 11, height = 5.5, family = "Times", useDingbats = FALSE)
draw_complete_figure()
invisible(dev.off())

png_file <- file.path(PROJECT_DIR, "exp2_delta_sensitivity_Figure6_clean.png")
grDevices::png(filename = png_file, width = 3300, height = 1650, res = 300, bg = "white")
draw_complete_figure()
invisible(dev.off())

cat("\n============================================================\n")
cat("FIGURES CREATED SUCCESSFULLY\n")
cat("============================================================\n\n")
cat("PDF: ", pdf_file, "\n")
cat("PNG: ", png_file, "\n")
