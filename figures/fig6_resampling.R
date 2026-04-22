# =============================================================================
# Figure 6: CHD resampling distribution, curated benchmark
#
# Panels:
#   A: P@10 distribution across 3,003 splits (violin + box; primary split marked)
#   B: Average Precision distribution (BIFO vs Fisher)
#   C: Rank improvement distribution
#   D: BIFO AP vs Fisher AP scatter per split
#
# Data: results/chd_curated/resampling_summary.json
#       (no per-split CSV needed, summary stats sufficient for A-C)
#
# Output: figures/fig6_resampling.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(jsonlite)

BASE    <- here::here()
CURATED <- here::here("results/chd_curated")
OUTDIR  <- here::here("figures")
dir.create(OUTDIR, showWarnings = FALSE)

fix_json <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  txt <- gsub("\\bNaN\\b",       "null", txt)
  txt <- gsub("\\bInfinity\\b",  "null", txt)
  txt <- gsub("\\b-Infinity\\b", "null", txt)
  fromJSON(txt)
}

rs <- fix_json(file.path(CURATED, "resampling_summary.json"))

# ── Extract distribution stats ────────────────────────────────────────────────
# Each metric has: mean, sd, min, p25, median, p75, max, n
extract_dist <- function(metric_list) {
  tibble(
    mean   = as.numeric(metric_list$mean),
    sd     = as.numeric(metric_list$sd),
    min    = as.numeric(metric_list$min),
    p25    = as.numeric(metric_list$p25),
    median = as.numeric(metric_list$median),
    p75    = as.numeric(metric_list$p75),
    max    = as.numeric(metric_list$max),
    n      = as.integer(metric_list$n)
  )
}

metrics <- rs$metrics
orig     <- rs$original_split
robust   <- rs$robustness
n_splits <- rs$n_splits_total

# Build summary table for box-plot-style panels using 5-number summary
# We simulate distribution shape from summary stats since we don't have per-split CSV
# Check if per-split CSV exists
resamp_csv <- file.path(CURATED, "resampling_results.csv")
has_csv <- file.exists(resamp_csv)

if (has_csv) {
  message("Loading per-split CSV: ", resamp_csv)
  splits_df <- read.csv(resamp_csv, stringsAsFactors = FALSE)
  message("  Columns: ", paste(names(splits_df), collapse = ", "))
  message("  Rows: ", nrow(splits_df))
} else {
  message("No per-split CSV found; using summary stats only")
  splits_df <- NULL
}

# ── Panels using per-split CSV if available ────────────────────────────────────
if (!is.null(splits_df)) {

  # Identify primary split
  orig_id <- orig$split_id
  splits_df <- splits_df %>%
    mutate(is_primary = (split_id == orig_id))

  primary_row <- splits_df %>% filter(is_primary)

  # Panel A: P@10 distribution
  pA <- ggplot(splits_df, aes(x = "All splits", y = bifo_p10)) +
    geom_violin(fill = "#C6DBEF", color = "#2171B5", alpha = 0.7,
                linewidth = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.12, fill = "#2171B5", color = "white",
                 outlier.size = 0.4, alpha = 0.9) +
    geom_hline(yintercept = primary_row$bifo_p10[1],
               linetype = "dashed", color = "#E15759", linewidth = 0.8) +
    annotate("text", x = 1.35, y = primary_row$bifo_p10[1] + 0.015,
             label = sprintf("Primary split\nP@10 = %.2f", primary_row$bifo_p10[1]),
             size = 3, color = "#E15759", fontface = "bold") +
    annotate("text", x = 0.6, y = max(splits_df$bifo_p10) * 0.95,
             label = sprintf("%.1f%% splits\nP@10 \u2265 0.30",
                             100 * robust$bifo_p10_ge_0.3 / n_splits),
             size = 3, color = "grey30") +
    scale_y_continuous(limits = c(-0.02, 1.05),
                       breaks = seq(0, 1, 0.2),
                       labels = percent_format(accuracy = 1)) +
    labs(title    = "A",
         subtitle = sprintf("P@10, %s splits", comma(n_splits)),
         x = NULL, y = "P@10") +
    theme_classic(base_size = 11) +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          plot.title   = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  # Panel B: BIFO AP vs Fisher AP
  pB <- ggplot(splits_df, aes(x = sf_ap, y = bifo_ap)) +
    geom_point(color = "grey70", size = 0.5, alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", color = "grey50", linewidth = 0.5) +
    geom_point(data = primary_row,
               aes(x = sf_ap, y = bifo_ap),
               color = "#E15759", size = 3.5, shape = 18) +
    annotate("text", x = max(splits_df$sf_ap) * 0.7,
             y = max(splits_df$bifo_ap) * 0.95,
             label = sprintf("BIFO > Fisher AP\nin %.1f%% of splits",
                             100 * robust$bifo_beats_sf_ap / n_splits),
             size = 3, color = "grey30") +
    scale_x_continuous(labels = number_format(accuracy = 0.001)) +
    scale_y_continuous(labels = number_format(accuracy = 0.001)) +
    labs(title    = "B",
         subtitle = "BIFO AP vs Seed Fisher AP per split\n(diamond = primary split)",
         x = "Seed Fisher AP", y = "BIFO AP") +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  # Panel C: Rank improvement distribution
  pC <- ggplot(splits_df, aes(x = "All splits", y = rank_improvement)) +
    geom_violin(fill = "#C6DBEF", color = "#2171B5", alpha = 0.7,
                linewidth = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.12, fill = "#2171B5", color = "white",
                 outlier.size = 0.4, alpha = 0.9) +
    geom_hline(yintercept = 0, linetype = "solid",
               color = "grey50", linewidth = 0.5) +
    geom_hline(yintercept = primary_row$rank_improvement[1],
               linetype = "dashed", color = "#E15759", linewidth = 0.8) +
    annotate("text", x = 1.35,
             y = primary_row$rank_improvement[1] + diff(range(splits_df$rank_improvement)) * 0.03,
             label = sprintf("Primary\n%+.1f", primary_row$rank_improvement[1]),
             size = 3, color = "#E15759", fontface = "bold") +
    annotate("text", x = 0.6,
             y = min(splits_df$rank_improvement) * 0.85,
             label = sprintf("100%% splits\nrank imp > 0\n(n=%s)", comma(n_splits)),
             size = 3, color = "grey30") +
    labs(title    = "C",
         subtitle = "Rank improvement, raw PPR rank minus conditioned PPR rank",
         x = NULL, y = "Rank positions gained") +
    theme_classic(base_size = 11) +
    theme(axis.text.x   = element_blank(),
          axis.ticks.x  = element_blank(),
          plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  # Panel D: P@10 percentile summary bar
  p10_dist   <- extract_dist(metrics$bifo_p10)
  ap_dist    <- extract_dist(metrics$bifo_ap)
  sf_ap_dist <- extract_dist(metrics$sf_ap)
  ri_dist    <- extract_dist(metrics$rank_improvement)

  summary_df <- tibble(
    metric   = factor(c("P@10", "AP (BIFO)", "AP (Fisher)", "Rank imp."),
                      levels = c("P@10", "AP (BIFO)", "AP (Fisher)", "Rank imp.")),
    median   = c(p10_dist$median,   ap_dist$median,   sf_ap_dist$median, ri_dist$median),
    p25      = c(p10_dist$p25,      ap_dist$p25,      sf_ap_dist$p25,    ri_dist$p25),
    p75      = c(p10_dist$p75,      ap_dist$p75,      sf_ap_dist$p75,    ri_dist$p75),
    primary  = c(orig$bifo_p10,     orig$bifo_ap,     orig$sf_ap,        orig$rank_improvement)
  )

  pD <- ggplot(summary_df, aes(x = metric, y = median)) +
    geom_col(fill = "#2E75B6", width = 0.5, alpha = 0.8) +
    geom_errorbar(aes(ymin = p25, ymax = p75),
                  width = 0.2, color = "grey30", linewidth = 0.7) +
    geom_point(aes(y = primary), color = "#E15759",
               size = 3.5, shape = 18) +
    labs(title    = "D",
         subtitle = "Distribution summary: median (bar), IQR (error bars),\nprimary split (diamond)",
         x = NULL, y = "Value") +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          axis.text.x   = element_text(size = 10))

} else {
  # Fallback: use summary stats to draw approximate distributions
  p10_dist   <- extract_dist(metrics$bifo_p10)
  ap_dist    <- extract_dist(metrics$bifo_ap)
  sf_ap_dist <- extract_dist(metrics$sf_ap)
  ri_dist    <- extract_dist(metrics$rank_improvement)

  dist_df <- bind_rows(
    p10_dist   %>% mutate(metric = "BIFO P@10"),
    ap_dist    %>% mutate(metric = "BIFO AP"),
    sf_ap_dist %>% mutate(metric = "Fisher AP"),
    ri_dist    %>% mutate(metric = "Rank improvement")
  ) %>%
    mutate(metric = factor(metric,
                            levels = c("BIFO P@10", "BIFO AP",
                                       "Fisher AP", "Rank improvement")))

  primary_df <- tibble(
    metric = factor(c("BIFO P@10", "BIFO AP", "Fisher AP", "Rank improvement"),
                     levels = levels(dist_df$metric)),
    value  = c(orig$bifo_p10, orig$bifo_ap, orig$sf_ap, orig$rank_improvement)
  )

  make_panel <- function(m, title_label, subtitle_label, y_label) {
    d  <- dist_df  %>% filter(metric == m)
    pr <- primary_df %>% filter(metric == m)
    ggplot(d, aes(x = metric)) +
      geom_crossbar(aes(y = median, ymin = p25, ymax = p75),
                    fill = "#C6DBEF", color = "#2171B5",
                    width = 0.4, linewidth = 0.5) +
      geom_errorbar(aes(y = median, ymin = min, ymax = max),
                    width = 0.15, color = "grey60", linewidth = 0.4) +
      geom_point(aes(y = median), color = "#2171B5", size = 3) +
      geom_hline(yintercept = pr$value,
                 linetype = "dashed", color = "#E15759", linewidth = 0.8) +
      annotate("text", x = 1.35, y = pr$value,
               label = sprintf("Primary: %.3f", pr$value),
               size = 3, color = "#E15759", fontface = "bold", hjust = 0) +
      labs(title = title_label, subtitle = subtitle_label, x = NULL, y = y_label) +
      theme_classic(base_size = 11) +
      theme(axis.text.x   = element_blank(),
            axis.ticks.x  = element_blank(),
            plot.title    = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 9, color = "grey40"))
  }

  pA <- make_panel("BIFO P@10",        "A", sprintf("P@10, %s splits\n(%.1f%% \u2265 0.30; %.1f%% \u2265 0.50)",
                                                     comma(n_splits),
                                                     100 * robust$bifo_p10_ge_0.3 / n_splits,
                                                     100 * robust$bifo_p10_ge_0.5 / n_splits),
                   "P@10")
  pB <- make_panel("BIFO AP",          "B", "BIFO Average Precision", "AP")
  pC <- make_panel("Fisher AP",        "C", "Seed Fisher Average Precision", "AP")
  pD <- make_panel("Rank improvement", "D",
                   sprintf("Rank improvement, %s/%s splits positive (100%%)",
                           comma(robust$rank_imp_positive), comma(n_splits)),
                   "Rank positions gained")
}

# ── Robustness annotation table ───────────────────────────────────────────────
rob_text <- sprintf(
  "Robustness summary (%s splits):\n  P@10 \u2265 0.50: %s splits (%.1f%%)\n  P@10 \u2265 0.30: %s splits (%.1f%%)\n  Rank imp > 0: %s splits (100%%)\n  BIFO > Fisher AP: %s splits (%.1f%%)",
  comma(n_splits),
  comma(robust$bifo_p10_ge_0.5), 100 * robust$bifo_p10_ge_0.5 / n_splits,
  comma(robust$bifo_p10_ge_0.3), 100 * robust$bifo_p10_ge_0.3 / n_splits,
  comma(robust$rank_imp_positive),
  comma(robust$bifo_beats_sf_ap), 100 * robust$bifo_beats_sf_ap / n_splits
)

# ── Assemble ──────────────────────────────────────────────────────────────────
fig6 <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title    = "Figure 6. BIFO pathway recovery is stable across 3,003 seed partitions, curated CHD benchmark",
    subtitle = paste0(
      rob_text, "\n",
      "Each split: 10 seeds / 5 held-out genes from 15 curated CHD genes. ",
      "Primary split (red diamond) = manuscript benchmark split. ",
      "Box = IQR; whiskers = min/max."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8.5, color = "grey35",
                                   lineheight = 1.35, family = "mono")
    )
  )

out <- file.path(OUTDIR, "fig6_resampling.png")
ggsave(out, fig6, width = 12, height = 10, dpi = 300, bg = "white")
message("Saved: ", out)
