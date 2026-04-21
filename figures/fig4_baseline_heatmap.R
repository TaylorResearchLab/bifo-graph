# =============================================================================
# Figure 4: Baseline method comparison heatmap — KF-CHD and KF-NBL
#
# Panels:
#   A: KF-CHD heatmap — method × metric
#   B: KF-NBL heatmap — method × metric
#
# Data: kf_{cohort}_results/baseline_comparison.json
# Output: figures/fig4_baseline_heatmap.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(jsonlite)

BASE   <- "/mnt/isilon/taylor_lab/data/projects/BIFO_2026"
OUTDIR <- file.path(BASE, "figures")
dir.create(OUTDIR, showWarnings = FALSE)

fix_json <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  txt <- gsub("\\bNaN\\b",       "null", txt)
  txt <- gsub("\\bInfinity\\b",  "null", txt)
  txt <- gsub("\\b-Infinity\\b", "null", txt)
  fromJSON(txt)
}

# ── Load baseline JSON for both cohorts ───────────────────────────────────────
load_baseline_json <- function(path, cohort) {
  d <- fix_json(path)
  m <- d$methods  # data.frame with one row per method
  tibble(
    cohort           = cohort,
    method           = m$method,
    p10              = as.numeric(m$precision_at_10),
    p20              = as.numeric(m$precision_at_20),
    p50              = as.numeric(m$precision_at_50),
    ap               = as.numeric(m$average_precision),
    mean_rank        = as.numeric(m$mean_ref_rank)
  )
}

chd_bl <- load_baseline_json(
  file.path(BASE, "bifo-graph/results/kf_chd/baseline_comparison.json"), "KF-CHD")
nbl_bl <- load_baseline_json(
  file.path(BASE, "bifo-graph/results/kf_nbl/baseline_comparison.json"), "KF-NBL")

df <- bind_rows(chd_bl, nbl_bl)

# ── Method labels ─────────────────────────────────────────────────────────────
method_map <- c(
  seed_fisher          = "Seed Fisher\n(B1, corrected)",
  neighborhood_fisher  = "Neighborhood\nFisher (B2)",
  raw_ppr_gsea         = "Raw PPR\nGSEA (B3)",
  cond_ppr_gsea        = "Cond PPR\nGSEA (B3b)",
  bifo_full            = "BIFO\nFull-arm (B4)"
)
df <- df %>%
  mutate(method_label = factor(recode(method, !!!method_map),
                                levels = method_map))

# ── Long format for heatmap ───────────────────────────────────────────────────
df_long <- df %>%
  pivot_longer(cols = c(p10, p20, p50, ap, mean_rank),
               names_to = "metric", values_to = "value") %>%
  mutate(
    metric_label = factor(
      recode(metric,
             p10       = "P@10",
             p20       = "P@20",
             p50       = "P@50",
             ap        = "Avg\nPrecision",
             mean_rank = "Mean Ref\nRank"),
      levels = c("P@10", "P@20", "P@50", "Avg\nPrecision", "Mean Ref\nRank")
    )
  )

# Normalize within cohort×metric for color scaling
# For mean_rank: invert (lower rank = better = higher scaled value)
df_long <- df_long %>%
  group_by(cohort, metric) %>%
  mutate(
    norm = case_when(
      metric == "mean_rank" ~ 1 - (value / max(value, na.rm = TRUE)),
      max(value, na.rm = TRUE) == 0 ~ 0,
      TRUE ~ value / max(value, na.rm = TRUE)
    )
  ) %>%
  ungroup() %>%
  mutate(
    disp_label = case_when(
      metric == "mean_rank" ~ comma(round(value, 0)),
      TRUE                  ~ sprintf("%.4f", value)
    )
  )

# ── Plot function ─────────────────────────────────────────────────────────────
make_heatmap <- function(data, cohort_name) {
  d <- data %>% filter(cohort == cohort_name)

  # Annotation: note WP_CILIOPATHIES rank=1 for seed_fisher and bifo_full
  d <- d %>%
    mutate(note = ifelse(
      metric %in% c("p10", "mean_rank") &
        method %in% c("seed_fisher", "bifo_full"),
      "\u2605", ""   # star for cilia-rank-1 methods
    ))

  ggplot(d, aes(x = metric_label, y = method_label, fill = norm)) +
    geom_tile(color = "white", linewidth = 1.0) +
    geom_text(aes(label = disp_label),
              size = 2.9, fontface = "bold", color = "grey10") +
    scale_fill_gradientn(
      colors   = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#084594"),
      na.value = "grey92",
      limits   = c(0, 1),
      name     = "Scaled\nperformance\n(within metric)"
    ) +
    scale_x_discrete(position = "top") +
    labs(title = cohort_name, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid    = element_blank(),
      axis.text.x   = element_text(face = "bold", size = 9.5),
      axis.text.y   = element_text(size = 9.5),
      plot.title    = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.title  = element_text(size = 8),
      legend.text   = element_text(size = 8)
    )
}

pCHD <- make_heatmap(df_long, "KF-CHD")
pNBL <- make_heatmap(df_long, "KF-NBL")

fig4 <- pCHD / pNBL +
  plot_annotation(
    title    = "Figure 4. Baseline enrichment method comparison — KF-CHD and KF-NBL cohorts",
    subtitle = paste0(
      "Heatmap of five enrichment methods across five performance metrics ",
      "(colour scaled within each metric; darker = better).\n",
      "Metrics evaluated against the 17-pathway cilia reference set. ",
      "Mean Ref Rank: mean rank of reference pathways (lower = better, colour inverted).\n",
      "Note: WP_CILIOPATHIES ranks 1st under Seed Fisher (corrected) in both cohorts; BIFO ranks 43rd/2,130 (CHD) and 3rd/2,196 (NBL). ",
      "in both cohorts; see Table 8.1 for per-pathway ranks."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8.5, color = "grey35", lineheight = 1.35)
    )
  )

out <- file.path(OUTDIR, "fig4_baseline_heatmap.png")
ggsave(out, fig4, width = 10, height = 9, dpi = 300, bg = "white")
message("Saved: ", out)
