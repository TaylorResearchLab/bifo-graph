# =============================================================================
# Figure 2: Four-arm gene-level recovery, curated CHD benchmark
#
# Panels:
#   A: AUROC by arm
#   B: AUPRC by arm
#   C: Localization (mean PPR rank of held-out genes / n_nodes) by arm
#   D: Entropy by arm
#
# Arms: raw, metadata_filtered, conditioned, random_sparsification_control
# Data: results/chd_curated/results_full.json (and ablation, mech for context)
#
# Output: figures/fig2_gene_recovery.png
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

# ── Load all three graph arm results ─────────────────────────────────────────
files <- list(
  Full        = file.path(CURATED, "results_full.json"),
  Ablation    = file.path(CURATED, "results_ablation.json"),
  Mechanistic = file.path(CURATED, "results_mech.json")
)

# Extract gene-level metrics per arm per graph configuration
extract_arms <- function(json_path, graph_label) {
  d <- fix_json(json_path)
  arms <- c("raw", "metadata_filtered", "conditioned", "random_sparsification_control")
  bind_rows(lapply(arms, function(arm) {
    a <- d[[arm]]
    tibble(
      graph_arm  = graph_label,
      ppr_arm    = arm,
      auroc      = as.numeric(a$auroc),
      auprc      = as.numeric(a$auprc),
      localization = as.numeric(a$localization),
      entropy    = as.numeric(a$entropy)
    )
  }))
}

df <- bind_rows(
  extract_arms(files$Full,        "Full"),
  extract_arms(files$Ablation,    "Ablation"),
  extract_arms(files$Mechanistic, "Mechanistic")
)

# Clean up arm labels
arm_labels <- c(
  raw                           = "Raw",
  metadata_filtered             = "Metadata\nfiltered",
  conditioned                   = "BIFO\nconditioned",
  random_sparsification_control = "Random\nsparsification"
)

graph_labels <- c(
  Full        = "Full (all flow classes)",
  Ablation    = "Ablation (no bridge edges)",
  Mechanistic = "Mechanistic only"
)

df <- df %>%
  mutate(
    ppr_label   = factor(recode(ppr_arm,   !!!arm_labels),
                          levels = arm_labels),
    graph_label = factor(recode(graph_arm, !!!graph_labels),
                          levels = graph_labels)
  )

# Colors: highlight BIFO conditioned arm
arm_colors <- c(
  "Raw"                    = "#A9A9A9",
  "Metadata\nfiltered"     = "#C8A882",
  "BIFO\nconditioned"      = "#2E75B6",
  "Random\nsparsification" = "#D9D9D9"
)

theme_bifo <- theme_classic(base_size = 11) +
  theme(
    axis.text.x   = element_text(angle = 20, hjust = 1, size = 9),
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    legend.position = "none",
    strip.text    = element_text(face = "bold", size = 9.5)
  )

# Panel A: AUROC
pA <- ggplot(df, aes(x = ppr_label, y = auroc, fill = ppr_label)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", auroc)),
            vjust = -0.4, size = 3.0, fontface = "bold") +
  facet_wrap(~ graph_label, nrow = 1) +
  scale_fill_manual(values = arm_colors) +
  scale_y_continuous(limits = c(0.97, 1.005),
                     breaks = c(0.97, 0.98, 0.99, 1.00),
                     oob = squish) +
  labs(title = "A", subtitle = "AUROC, held-out gene recovery",
       x = NULL, y = "AUROC") +
  theme_bifo

# Panel B: AUPRC
auprc_max <- max(df$auprc, na.rm = TRUE)
pB <- ggplot(df, aes(x = ppr_label, y = auprc, fill = ppr_label)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.3f", auprc)),
            vjust = -0.4, size = 3.0, fontface = "bold") +
  facet_wrap(~ graph_label, nrow = 1) +
  scale_fill_manual(values = arm_colors) +
  scale_y_continuous(limits = c(0, auprc_max * 1.25)) +
  labs(title = "B", subtitle = "AUPRC, held-out gene recovery",
       x = NULL, y = "AUPRC") +
  theme_bifo

# Panel C: Localization (lower = better, seeds more tightly localized)
pC <- ggplot(df, aes(x = ppr_label, y = localization, fill = ppr_label)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.4f", localization)),
            vjust = -0.4, size = 3.0, fontface = "bold") +
  facet_wrap(~ graph_label, nrow = 1) +
  scale_fill_manual(values = arm_colors) +
  scale_y_continuous(limits = c(0, max(df$localization, na.rm=TRUE) * 1.25)) +
  labs(title = "C",
       subtitle = "Localization, mean held-out rank / n nodes\n(lower = more concentrated near seeds)",
       x = NULL, y = "Localization") +
  theme_bifo

# Panel D: Entropy (lower = more concentrated signal)
pD <- ggplot(df, aes(x = ppr_label, y = entropy, fill = ppr_label)) +
  geom_col(width = 0.6, color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", entropy)),
            vjust = -0.4, size = 3.0, fontface = "bold") +
  facet_wrap(~ graph_label, nrow = 1) +
  scale_fill_manual(values = arm_colors) +
  scale_y_continuous(limits = c(0, max(df$entropy, na.rm=TRUE) * 1.15)) +
  labs(title = "D",
       subtitle = "PPR score entropy, bits\n(lower = more concentrated signal)",
       x = NULL, y = "Entropy (bits)") +
  theme_bifo

# ── Legend panel ──────────────────────────────────────────────────────────────
legend_df <- tibble(
  label = factor(names(arm_colors), levels = names(arm_colors)),
  x = 1:4, y = 1
)
legend_p <- ggplot(legend_df, aes(x = label, y = y, fill = label)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = arm_colors, name = "PPR arm") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.title = element_text(face = "bold", size = 9),
        legend.text  = element_text(size = 9)) +
  guides(fill = guide_legend(nrow = 1))

# Extract legend
get_legend <- function(p) {
  g <- ggplotGrob(p)
  leg <- g$grobs[sapply(g$grobs, function(x) x$name) == "guide-box"][[1]]
  leg
}

# ── Assemble ──────────────────────────────────────────────────────────────────
fig2 <- (pA / pB / pC / pD) +
  plot_annotation(
    title    = "Figure 2. Four-arm gene-level recovery, curated CHD benchmark",
    subtitle = paste0(
      "Each panel shows four PPR propagation arms (columns) across three graph configurations (facets).\n",
      "Graph arms: Full = all BIFO flow classes; Ablation = no Pathway Contribution bridge edges; ",
      "Mechanistic = mechanistic-only edges.\n",
      "PPR arms: Raw (unconditioned), Metadata-filtered, BIFO-conditioned (primary), ",
      "Random sparsification control.\n",
      "AUROC near-ceiling across all arms reflects the strongly connected seed neighborhood; ",
      "AUPRC and entropy are more discriminating."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 8.5, color = "grey35", lineheight = 1.35)
    )
  )

out <- file.path(OUTDIR, "fig2_gene_recovery.png")
ggsave(out, fig2, width = 13, height = 14, dpi = 300, bg = "white")
message("Saved: ", out)
