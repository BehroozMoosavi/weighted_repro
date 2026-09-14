# ==============================================================================
# Script: heatmap.R
# Description: Generates publication-ready heatmaps from precomputed power 
#              matrices (Efficient_INT.csv) for the private linear regression 
#              simulation study. Does not rerun underlying simulations.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Paths & File Verification
# ------------------------------------------------------------------------------

PROJECT_DIR <- path.expand("~/R_Simuls/Linearegression")
RESULTS_DIR <- file.path(PROJECT_DIR, "results_eff")
CSV_FILE    <- file.path(RESULTS_DIR, "Efficient_INT.csv")
OUTPUT_PNG  <- file.path(RESULTS_DIR, "Efficient_power_heatmap.png")
OUTPUT_PDF  <- file.path(RESULTS_DIR, "Efficient_power_heatmap.pdf")

if (!file.exists(CSV_FILE)) {
  stop(paste0("CSV file not found:\n", CSV_FILE))
}

power_df <- read.csv(CSV_FILE, row.names = 1, check.names = FALSE)
cat("\nData read from:\n", CSV_FILE, "\n\n")
print(power_df)

# ------------------------------------------------------------------------------
# 2. Data Preparation
# ------------------------------------------------------------------------------

power_matrix <- as.matrix(power_df)
storage.mode(power_matrix) <- "numeric"

n_values    <- as.numeric(colnames(power_matrix))
beta_values <- as.numeric(rownames(power_matrix))

cat("\nSample sizes:\n")
print(n_values)
cat("\nTrue beta1 values:\n")
print(beta_values)
cat("\nPower matrix:\n")
print(power_matrix)

# ------------------------------------------------------------------------------
# 3. Heatmap Rendering Routine
# ------------------------------------------------------------------------------

draw_power_heatmap <- function(
    power_matrix,
    n_values,
    beta_values
) {
  nr <- nrow(power_matrix)
  nc <- ncol(power_matrix)
  
  # Reverse row ordering so beta1 = 1 is placed at the top
  Z <- power_matrix[nr:1, , drop = FALSE]
  beta_plot <- rev(beta_values)
  
  cols <- colorRampPalette(
    c(
      "#fff5f0",
      "#fee0d2",
      "#fcbba1",
      "#fc9272",
      "#fb6a4a",
      "#ef3b2c",
      "#cb181d",
      "#99000d"
    )
  )(100)
  
  # Render heat grid
  image(
    x = seq_len(nc),
    y = seq_len(nr),
    z = t(Z),
    col = cols,
    zlim = c(0, 1),
    axes = FALSE,
    xlab = "Number of Samples (N)",
    ylab = expression(beta[1]),
    main = "Efficient Method",
    useRaster = TRUE
  )
  
  axis(side = 1, at = seq_len(nc), labels = n_values, cex.axis = 1.1)
  axis(side = 2, at = seq_len(nr), labels = beta_plot, las = 1, cex.axis = 1.1)
  
  # Annotate cell values with contrast-adjusted text color
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      value <- Z[i, j]
      if (is.finite(value)) {
        label <- sprintf("%.3f", value)
        text_col <- if (value >= 0.65) "white" else "black"
        text(
          x = j,
          y = i,
          labels = label,
          cex = 1.05,
          font = 2,
          col = text_col
        )
      }
    }
  }
  
  box()
  
  # Render dedicated color bar
  usr <- par("usr")
  xleft  <- usr[2] + 0.25
  xright <- usr[2] + 0.45
  y_breaks <- seq(usr[3], usr[4], length.out = 101)
  
  for (k in seq_len(100)) {
    rect(
      xleft = xleft,
      ybottom = y_breaks[k],
      xright = xright,
      ytop = y_breaks[k + 1],
      col = cols[k],
      border = NA,
      xpd = TRUE
    )
  }
  
  power_ticks <- seq(0, 1, by = 0.2)
  y_ticks <- usr[3] + power_ticks * (usr[4] - usr[3])
  
  axis(
    side = 4,
    at = y_ticks,
    labels = sprintf("%.1f", power_ticks),
    las = 1,
    line = -1,
    xpd = TRUE
  )
  
  mtext("Power", side = 4, line = 3.5, xpd = TRUE)
}

# ------------------------------------------------------------------------------
# 4. Output Export
# ------------------------------------------------------------------------------

png(OUTPUT_PNG, width = 1700, height = 1200, res = 200)
par(mar = c(5, 5, 4, 8))
draw_power_heatmap(power_matrix, n_values, beta_values)
dev.off()

pdf(OUTPUT_PDF, width = 8, height = 5.8)
par(mar = c(5, 5, 4, 8))
draw_power_heatmap(power_matrix, n_values, beta_values)
dev.off()

cat("\n============================================================\n")
cat("POWER HEATMAP CREATED\n")
cat("============================================================\n")
cat(sprintf("Input CSV : %s\n", CSV_FILE))
cat(sprintf("PNG       : %s\n", OUTPUT_PNG))
cat(sprintf("PDF       : %s\n", OUTPUT_PDF))
cat("============================================================\n")
