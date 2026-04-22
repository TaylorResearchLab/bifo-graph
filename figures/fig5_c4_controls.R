# =============================================================================
# Figure 5: C4 pathway-split controls, Notch and MAPK
#
# Panels:
#   A: P@10, Notch vs MAPK
#   B: Enrichment@10, Notch vs MAPK
#   C: Mean CHD pathway rank, conditioned vs raw (both controls)
#   D: Rank improvement, Notch vs MAPK
#
# Data: results/chd_curated/c4_notch/pathway_metrics.json
#       results/chd_curated/c4_mapk/pathway_metrics.json
# Output: figures/fig5_c4_controls.png
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

pm_notch <- fix_json(file.path(CURATED, "c4_notch/pathway_metrics.json"))
pm_mapk  <- fix_json(file.path(CURATED, "c4_mapk/pathway_metrics.json"))

# Load CHD benchmark metrics for cross-reference in subtitle
CHD_BENCH <- here::here("bifo-graph/results/chd_benchmark")
pm_chd    <- fix_json(file.path(CHD_BENCH, "pathway_metrics_full.json"))
chd_rank_imp <- sprintf("+%.1f", pm_chd$metrics$rank_improvement_cond_vs_raw)

# Also load c4 results.json files for gene-level stats
rn <- fix_json(file.path(CURATED, "c4_notch/results.json"))
rm <- fix_json(file.path(CURATED, "c4_mapk/results.json"))

# ── Build metrics data frame ──────────────────────────────────────────────────
c4_df <- bind_rows(
  tibble(
    control      = "C4/Notch\n(recovery)",
    task         = "Recovery",
    n_seeds      = rn$parameters$n_seeds_provided,
    n_heldout    = rn$parameters$n_heldout_provided,
    n_reference  = pm_notch$metrics$n_chd_in_reference,
    p10          = pm_notch$metrics$top10_precision,
    enrich10     = pm_notch$metrics$top10_enrichment_ratio,
    mean_cond    = pm_notch$metrics$mean_rank_chd_cond,
    mean_raw     = pm_notch$metrics$mean_rank_chd_raw,
    rank_imp     = pm_notch$metrics$rank_improvement_cond_vs_raw,
    auroc_cond   = rn$conditioned$auroc,
    auprc_cond   = rn$conditioned$auprc
  ),
  tibble(
    control      = "C4/MAPK\n(orthogonal)",
    task         = "Orthogonal",
    n_seeds      = rm$parameters$n_seeds_provided,
    n_heldout    = rm$parameters$n_heldout_provided,
    n_reference  = pm_mapk$metrics$n_chd_in_reference,
    p10          = pm_mapk$metrics$top10_precision,
    enrich10     = pm_mapk$metrics$top10_enrichment_ratio,
    mean_cond    = pm_mapk$metrics$mean_rank_chd_cond,
    mean_raw     = pm_mapk$metrics$mean_rank_chd_raw,
    rank_imp     = pm_mapk$metrics$rank_improvement_cond_vs_raw,
    auroc_cond   = rm$conditioned$auroc,
    auprc_cond   = rm$conditioned$auprc
  )
)

c4_df$control <- factor(c4_df$control,
                         levels = c("C4/Notch\n(recovery)", "C4/MAPK\n(orthogonal)"))

ctrl_colors <- c("C4/Notch\n(recovery)"   = "#2E75B6",
                 "C4/MAPK\n(orthogonal)"  = "#ED7D31")

theme_c4 <- theme_classic(base_size = 12) +
  theme(legend.position  = "none",
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 9, color = "grey40"),
        axis.text.x      = element_text(size = 10.5))

# Panel A: P@10
pA <- ggplot(c4_df, aes(x = control, y = p10, fill = control)) +
  geom_col(width = 0.5, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", p10)),
            vjust = -0.4, size = 4.0, fontface = "bold") +
  scale_fill_manual(values = ctrl_colors) +
  scale_y_continuous(limits = c(0, 1.0),
                     breaks = seq(0, 1, 0.2),
                     labels = percent_format(accuracy = 1)) +
  labs(title = "A", subtitle = "Precision@10",
       x = NULL, y = "P@10") +
  theme_c4

# Panel B: Enrichment@10
pB <- ggplot(c4_df, aes(x = control, y = enrich10, fill = control)) +
  geom_col(width = 0.5, color = "white", linewidth = 0.4) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "grey60", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f\u00d7", enrich10)),
            vjust = -0.4, size = 4.0, fontface = "bold") +
  scale_fill_manual(values = ctrl_colors) +
  scale_y_continuous(limits = c(0, max(c4_df$enrich10) * 1.25)) +
  labs(title = "B", subtitle = "Enrichment@10 (\u00d7 background rate)",
       x = NULL, y = "Enrichment@10") +
  theme_c4

# Panel C: Mean rank conditioned vs raw, grouped bars
rank_long <- c4_df %>%
  select(control, mean_cond, mean_raw) %>%
  pivot_longer(cols = c(mean_cond, mean_raw),
               names_to = "ppr_arm", values_to = "mean_rank") %>%
  mutate(ppr_arm = recode(ppr_arm,
                          "mean_cond" = "BIFO conditioned",
                          "mean_raw"  = "Raw PPR"))

pC <- ggplot(rank_long,
             aes(x = control, y = mean_rank,
                 fill = ppr_arm, group = ppr_arm)) +
  geom_col(position = position_dodge(0.6), width = 0.5,
           color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", mean_rank)),
            position = position_dodge(0.6),
            vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c("BIFO conditioned" = "#2E75B6",
                               "Raw PPR"          = "#C8C8C8"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, max(rank_long$mean_rank) * 1.25),
                     labels = comma) +
  labs(title    = "C",
       subtitle = "Mean CHD pathway rank, conditioned vs raw PPR\n(lower = better)",
       x = NULL, y = "Mean rank") +
  theme_classic(base_size = 12) +
  theme(legend.position  = "bottom",
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 9, color = "grey40"),
        axis.text.x      = element_text(size = 10.5))

# Panel D: Rank improvement
ri_abs_max <- max(abs(c4_df$rank_imp), na.rm = TRUE)
pD <- ggplot(c4_df, aes(x = control, y = rank_imp, fill = control)) +
  geom_col(width = 0.5, color = "white", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%+.1f", rank_imp)),
            vjust = ifelse(c4_df$rank_imp >= 0, -0.4, 1.3),
            size = 4.0, fontface = "bold") +
  scale_fill_manual(values = ctrl_colors) +
  scale_y_continuous(limits = c(-ri_abs_max * 1.35,
                                 ri_abs_max * 0.25)) +
  labs(title    = "D",
       subtitle = "Rank improvement (raw PPR rank \u2212 conditioned PPR rank)",
       x = NULL, y = "Rank positions") +
  theme_c4

# ── Summary annotation ────────────────────────────────────────────────────────
notch_seeds <- c4_df$n_seeds[c4_df$task == "Recovery"]
notch_ref   <- c4_df$n_reference[c4_df$task == "Recovery"]
mapk_seeds  <- c4_df$n_seeds[c4_df$task == "Orthogonal"]
mapk_ref    <- c4_df$n_reference[c4_df$task == "Orthogonal"]

# ── Assemble ──────────────────────────────────────────────────────────────────
fig5 <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title    = "Figure 5. C4 pathway-split controls, curated CHD benchmark",
    subtitle = paste0(
      "Two pathway-family controls using seeds drawn from curated pathway gene lists rather\n",
      "than the CHD variant gene pool.\n",
      "C4/Notch (recovery task): ", notch_seeds, " Notch pathway seeds \u2192 ",
      notch_ref, " Notch/developmental reference pathways.\n",
      "C4/MAPK (orthogonal control): ", mapk_seeds, " MAPK seeds \u2192 ",
      mapk_ref, " MAPK reference pathways.\n",
      "A/B: BIFO recovers within-family signal (Notch 25\u00d7 enrichment) and appropriately\n",
      "attenuates orthogonal signal (MAPK 5.5\u00d7).\n",
      "C/D: With within-family seeds, raw PPR outranks conditioned PPR for the source family\n",
      "(negative rank improvement in both controls). This is the expected mechanism: BIFO\n",
      "conditioning redistributes rank mass toward cross-family bridges, which helps with\n",
      paste0("heterogeneous seed sets (see Fig 3: ", chd_rank_imp, " rank imp for curated CHD) but slightly\n"),
      "dilutes within-family concentration here."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 8.5, color = "grey35", lineheight = 1.35)
    )
  )

out <- file.path(OUTDIR, "fig5_c4_controls.png")
ggsave(out, fig5, width = 11, height = 11, dpi = 300, bg = "white")
message("Saved: ", out)
