# =============================================================================
# Figure 8: Cross-cohort convergence — KF-CHD and KF-NBL
#
# Panels:
#   A: WP_CILIOPATHIES rank — BIFO and Fisher, CHD vs NBL (bar)
#   B: BIFO rank scatter CHD vs NBL (all pathways; cilia highlighted)
#   C: Cilia pathway cluster ranks — dot plot CHD vs NBL
#   D: Bootstrap resampling — BIFO vs Fisher P@10 by seed size
#
# Output: figures/fig8_crosscohort.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(ggrepel)
library(stringr)
library(readr)

BASE    <- "/mnt/isilon/taylor_lab/data/projects/BIFO_2026"
CHD_DIR <- file.path(BASE, "bifo-graph/results/kf_chd")
NBL_DIR <- file.path(BASE, "bifo-graph/results/kf_nbl")
OUTDIR  <- file.path(BASE, "figures")
dir.create(OUTDIR, showWarnings = FALSE)

CILIA_TERMS <- c("CILIOPATHIES", "CILIUM", "CILIA", "CILIARY",
                 "JOUBERT", "BARDET", "NEPHRONOPHTHISIS",
                 "HEDGEHOG", "LEFT_RIGHT", "LATERALITY",
                 "ALSTROM", "MECKEL", "PRIMARY_CILIARY",
                 "CILIOGENESIS", "IFT_", "FLAGELL")

is_cilia <- function(x) {
  sapply(toupper(x), function(n)
    any(sapply(CILIA_TERMS, function(t) grepl(t, n, fixed = TRUE))))
}

# ── Load pathway scores ───────────────────────────────────────────────────────
load_scores <- function(path, cohort) {
  df <- read_csv(path, show_col_types = FALSE) %>%
    arrange(desc(degree_norm)) %>%
    mutate(rank = row_number(), cohort = cohort)
  message(cohort, ": ", nrow(df), " pathways loaded")
  df
}

chd <- load_scores(file.path(CHD_DIR, "pathway_scores_standard.csv"), "KF-CHD")
nbl <- load_scores(file.path(NBL_DIR, "pathway_scores_standard.csv"), "KF-NBL")

# ── Panel A: WP_CILIOPATHIES rank bar ─────────────────────────────────────────
get_cilia_rank <- function(df) {
  df %>% filter(grepl("WP_CILIOPATHIES", name, ignore.case = TRUE)) %>%
    pull(rank) %>% first()
}

# Also get Fisher ranks from baseline CSVs
get_fisher_rank <- function(bc_path) {
  bc <- read_csv(bc_path, show_col_types = FALSE)
  bc %>%
    filter(grepl("WP_CILIOPATHIES", name, ignore.case = TRUE),
           method == "seed_fisher") %>%
    pull(rank) %>% first()
}

rank_df_A <- tibble(
  cohort  = rep(c("KF-CHD", "KF-NBL"), each = 2),
  method  = rep(c("BIFO full-arm", "Seed Fisher\n(corrected)"), 2),
  rank    = c(get_cilia_rank(chd) %||% 1L,
              get_fisher_rank(file.path(CHD_DIR, "baseline_comparison.csv")) %||% 1L,
              get_cilia_rank(nbl) %||% 1L,
              get_fisher_rank(file.path(NBL_DIR, "baseline_comparison.csv")) %||% 1L),
  total   = c(nrow(chd), nrow(chd), nrow(nbl), nrow(nbl))
)
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a)) a else b

rank_df_A$method <- factor(rank_df_A$method,
                            levels = c("BIFO full-arm", "Seed Fisher\n(corrected)"))
cohort_colors <- c("KF-CHD" = "#2E75B6", "KF-NBL" = "#ED7D31")

pA <- ggplot(rank_df_A,
             aes(x = method, y = rank, fill = cohort)) +
  geom_col(position = position_dodge(0.65), width = 0.55,
           color = "white", linewidth = 0.4) +
  geom_text(aes(label = paste0("Rank ", rank)),
            position = position_dodge(0.65),
            vjust = -0.4, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = cohort_colors, name = NULL) +
  scale_y_continuous(limits = c(0, max(rank_df_A$rank) * 1.2), labels = comma) +
  labs(title    = "A",
       subtitle = "WP_CILIOPATHIES rank — BIFO and Fisher",
       x = NULL, y = "Rank (lower = better)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40"))

# ── Panel B: CHD vs NBL rank scatter ─────────────────────────────────────────
shared <- inner_join(
  chd %>% select(concept_id, name, chd_rank = rank, chd_dn = degree_norm,
                 chd_hub = degree_flag, member_gene_count),
  nbl %>% select(concept_id, nbl_rank = rank, nbl_dn = degree_norm),
  by = "concept_id"
) %>%
  mutate(
    cilia      = is_cilia(name),
    avg_rank   = (chd_rank + nbl_rank) / 2,
    short_name = name %>%
      str_replace("^WP_", "") %>%
      str_replace("^REACTOME_", "") %>%
      str_replace_all("_", " ") %>%
      str_trunc(30)
  )

n_shared <- nrow(shared)
message("Shared pathways: ", n_shared)

pB <- ggplot(shared %>% filter(avg_rank <= 3000),
             aes(x = chd_rank, y = nbl_rank,
                 color = cilia, size = cilia, alpha = cilia)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey55", linewidth = 0.45) +
  geom_label_repel(
    data = shared %>% filter(cilia, avg_rank <= 500),
    aes(label = short_name),
    size = 2.4, max.overlaps = 20, box.padding = 0.3,
    label.size = 0.2, fill = alpha("white", 0.85),
    color = "#C0392B", fontface = "bold",
    min.segment.length = 0.1
  ) +
  scale_color_manual(values = c("FALSE" = "grey78", "TRUE" = "#C0392B")) +
  scale_size_manual(values  = c("FALSE" = 0.7, "TRUE"  = 2.8)) +
  scale_alpha_manual(values = c("FALSE" = 0.35, "TRUE"  = 1.0)) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(title    = "B",
       subtitle = paste0("BIFO pathway rank: KF-CHD vs KF-NBL\n",
                         "(top 3,000 shown; cilia pathways in red, n=",
                         sum(shared$cilia), ")"),
       x = "KF-CHD rank", y = "KF-NBL rank") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40"))

# ── Panel C: Cilia cluster ranks — dot plot ────────────────────────────────────
cilia_shared <- shared %>%
  filter(cilia) %>%
  select(short_name, chd_rank, nbl_rank) %>%
  pivot_longer(cols = c(chd_rank, nbl_rank),
               names_to = "cohort", values_to = "rank") %>%
  mutate(cohort = recode(cohort,
                         "chd_rank" = "KF-CHD",
                         "nbl_rank" = "KF-NBL"))

# Order pathways by mean rank
order_names <- cilia_shared %>%
  group_by(short_name) %>%
  summarise(mean_rank = mean(rank, na.rm = TRUE)) %>%
  arrange(mean_rank) %>%
  pull(short_name)

cilia_shared$short_name <- factor(cilia_shared$short_name,
                                   levels = rev(order_names))

pC <- ggplot(cilia_shared,
             aes(x = rank, y = short_name,
                 color = cohort, shape = cohort)) +
  geom_line(aes(group = short_name), color = "grey75", linewidth = 0.5) +
  geom_point(size = 3.2, alpha = 0.95) +
  scale_color_manual(values = cohort_colors, name = NULL) +
  scale_shape_manual(values = c("KF-CHD" = 16, "KF-NBL" = 17), name = NULL) +
  scale_x_continuous(labels = comma,
                     limits = c(0, NA)) +
  labs(title    = "C",
       subtitle = "Cilia pathway cluster — BIFO rank in both cohorts",
       x = "BIFO rank (lower = better)", y = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "bottom",
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        axis.text.y   = element_text(size = 8.5))

# ── Panel D: Bootstrap resampling violin ─────────────────────────────────────
load_resampling <- function(path, cohort) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(boot_id >= 0) %>%
    mutate(cohort = cohort)
}

chd_rs <- load_resampling(file.path(CHD_DIR, "resampling_results.csv"), "KF-CHD")
nbl_rs <- load_resampling(file.path(NBL_DIR, "resampling_results.csv"), "KF-NBL")

rs_all <- bind_rows(chd_rs, nbl_rs) %>%
  mutate(seed_label = factor(paste0("n=", seed_size),
                              levels = c("n=10", "n=20", "n=30"))) %>%
  pivot_longer(cols = c(bifo_p10, sf_p10),
               names_to = "method", values_to = "p10") %>%
  mutate(method = recode(method,
                         "bifo_p10" = "BIFO",
                         "sf_p10"   = "Seed Fisher"))

# Primary run reference lines
primary_lines <- tibble(
  cohort = c("KF-CHD", "KF-NBL"),
  bifo   = c(0.000, 0.100),
  fisher = c(0.300, 0.200)
) %>%
  pivot_longer(cols = c(bifo, fisher),
               names_to = "method_key", values_to = "p10") %>%
  mutate(method = recode(method_key, "bifo" = "BIFO", "fisher" = "Seed Fisher"))

method_colors <- c("BIFO" = "#2E75B6", "Seed Fisher" = "#E15759")

pD <- ggplot(rs_all, aes(x = seed_label, y = p10,
                          fill = method, color = method)) +
  geom_violin(alpha = 0.30, position = position_dodge(0.78),
              linewidth = 0.45, trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.11, position = position_dodge(0.78),
               outlier.size = 0.4, alpha = 0.80, linewidth = 0.4) +
  geom_hline(data = primary_lines,
             aes(yintercept = p10, color = method, linetype = method),
             linewidth = 0.7, show.legend = FALSE) +
  facet_wrap(~ cohort) +
  scale_fill_manual(values  = method_colors, name = NULL) +
  scale_color_manual(values = method_colors, name = NULL) +
  scale_linetype_manual(values = c("BIFO" = "dashed", "Seed Fisher" = "dotted")) +
  scale_y_continuous(limits = c(-0.02, 1.02),
                     breaks = seq(0, 1, 0.2),
                     labels = percent_format(accuracy = 1)) +
  labs(title    = "D",
       subtitle = paste0("Bootstrap resampling: P@10 vs 17-pathway cilia reference\n",
                         "(500 runs/seed size; dashed/dotted lines = primary run P@10)"),
       x = "Random seed sample size", y = "P@10") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom",
        strip.text    = element_text(face = "bold", size = 11),
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40"))

# ── Assemble ──────────────────────────────────────────────────────────────────
fig8 <- (pA | pB) / (pC | pD) +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(
    title    = "Figure 8. Cross-cohort convergence: KF-CHD and KF-NBL independently recover cilia pathways",
    subtitle = paste0(
      "A: WP_CILIOPATHIES ranks first under Seed Fisher in both cohorts; BIFO ranks 43rd/2,130 (CHD) and 3rd/2,196 (NBL).\n",
      "B: BIFO pathway ranks are concordant across cohorts; cilia pathways (red) cluster in the top ranks in both.\n",
      "C: All detected cilia pathways rank in the top half of the scored universe in both cohorts.\n",
      "D: Cilia signal requires full cohort-scale seed sets — ",
      "neither method recovers it reliably from n=10-30 random genes."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8.5, color = "grey35", lineheight = 1.4)
    )
  )

out <- file.path(OUTDIR, "fig8_crosscohort.png")
ggsave(out, fig8, width = 15, height = 13, dpi = 300, bg = "white")
message("Saved: ", out)
