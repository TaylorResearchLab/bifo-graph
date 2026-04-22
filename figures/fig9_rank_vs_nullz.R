# =============================================================================
# Figure 9: Rank vs null_z landscape, KF-CHD and KF-NBL
#
# Scatter plot of degree_norm rank vs pathway-node rewiring null_z for all
# scored pathways in each cohort. Cilia reference pathways highlighted.
# Illustrates the rank vs statistical enrichment distinction: WP_CILIOPATHIES
# ranks 43rd by degree_norm but has the highest null_z in KF-CHD.
#
# Data: results/kf_chd/pathway_scores_standard.csv
#       results/kf_nbl/pathway_scores_standard.csv
#       data/cohorts/chd/kf_chd_cilia_reference.txt
#       data/cohorts/nbl/kf_nbl_cilia_reference.txt
#
# Output: figures/fig9_rank_vs_nullz.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)
library(ggrepel)

BASE   <- here::here()
KF_CHD <- here::here("bifo-graph/results/kf_chd")
KF_NBL <- here::here("bifo-graph/results/kf_nbl")
OUTDIR <- here::here("figures")
dir.create(OUTDIR, showWarnings = FALSE)

# ── Load data ─────────────────────────────────────────────────────────────────
chd <- read.csv(file.path(KF_CHD, "pathway_scores_standard.csv"),
                stringsAsFactors = FALSE)
nbl <- read.csv(file.path(KF_NBL, "pathway_scores_standard.csv"),
                stringsAsFactors = FALSE)

chd_ref <- readLines(file.path(BASE, "bifo-graph/data/cohorts/chd/kf_chd_cilia_reference.txt"))
nbl_ref <- readLines(file.path(BASE, "bifo-graph/data/cohorts/nbl/kf_nbl_cilia_reference.txt"))

# Add rank
chd <- chd %>% arrange(desc(degree_norm)) %>% mutate(rank = row_number())
nbl <- nbl %>% arrange(desc(degree_norm)) %>% mutate(rank = row_number())

# Mark cilia reference pathways
chd <- chd %>% mutate(cilia = concept_id %in% chd_ref)
nbl <- nbl %>% mutate(cilia = concept_id %in% nbl_ref)

# Key pathways to label
label_names <- c("WP_CILIOPATHIES", "WP_JOUBERT_SYNDROME",
                 "REACTOME_CILIUM_ASSEMBLY",
                 "WP_GENES_RELATED_TO_PRIMARY_CILIUM_DEVELOPMENT_BASED_ON_CRISPR",
                 "REACTOME_ANCHORING_OF_THE_BASAL_BODY_TO_THE_PLASMA_MEMBRANE")

chd <- chd %>% mutate(label_text = ifelse(name %in% label_names, name, NA))
nbl <- nbl %>% mutate(label_text = ifelse(name %in% label_names, name, NA))

# Shorten long names for labels
shorten <- function(x) {
  x <- gsub("REACTOME_", "R:", x)
  x <- gsub("WP_GENES_RELATED_TO_PRIMARY_CILIUM_DEVELOPMENT_BASED_ON_CRISPR",
            "WP_PRIMARY_CILIUM_CRISPR", x)
  x <- gsub("REACTOME_ANCHORING_OF_THE_BASAL_BODY_TO_THE_PLASMA_MEMBRANE",
            "R:ANCHORING_BASAL_BODY", x)
  x
}
chd <- chd %>% mutate(label_text = shorten(label_text))
nbl <- nbl %>% mutate(label_text = shorten(label_text))

# ── Plot function ──────────────────────────────────────────────────────────────
make_panel <- function(df, cohort_label, n_pathways, top_pathway = "WP_CILIOPATHIES") {
  
  top <- df %>% filter(name == top_pathway)
  
  ggplot(df, aes(x = rank, y = null_z)) +
    # All pathways, grey background
    geom_point(data = df %>% filter(!cilia),
               color = "grey80", size = 0.6, alpha = 0.6) +
    # Cilia pathways, colored
    geom_point(data = df %>% filter(cilia),
               color = "#2E75B6", size = 2.2, alpha = 0.9) +
    # Horizontal null reference
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.4) +
    # Labels for key pathways
    ggrepel::geom_label_repel(
      data        = df %>% filter(!is.na(label_text)),
      aes(label   = label_text),
      size        = 2.6,
      color       = "#1a4f7a",
      fill        = "white",
      label.size  = 0.2,
      box.padding = 0.4,
      max.overlaps = 20,
      segment.color = "#2E75B6",
      segment.size  = 0.3
    ) +
    scale_x_continuous(
      labels = comma,
      limits = c(1, n_pathways),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title    = cohort_label,
      x        = "Degree-norm rank (lower = better)",
      y        = "Pathway-node rewiring null_z"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      axis.title    = element_text(size = 10),
      axis.text     = element_text(size = 9)
    )
}

pA <- make_panel(chd, "KF-CHD  (1,276 seeds, 2,130 pathways)", nrow(chd))
pB <- make_panel(nbl, "KF-NBL  (1,395 seeds, 2,196 pathways)", nrow(nbl))

# ── Assemble ──────────────────────────────────────────────────────────────────
fig9 <- pA | pB +
  plot_annotation(
    title    = "Figure 9. Rank vs statistical enrichment, KF cohort pathway landscape",
    subtitle = paste0(
      "Each point is a scored pathway. Blue points: 16-pathway cilia reference set. ",
      "Degree-norm rank reflects absolute propagated signal; null_z reflects signal ",
      "relative to degree-preserving membership rewiring null (N=1,000 permutations).\n",
      "WP_CILIOPATHIES ranks 43rd/2,130 in KF-CHD (null_z = 41.2) and 3rd/2,196 in KF-NBL (null_z = 18.4) by degree_norm, ",
      "but shows the highest pathway-node null_z among well-calibrated pathways in KF-CHD, illustrating that rank and ",
      "statistical enrichment are complementary, not redundant, metrics."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8.5, color = "grey35", lineheight = 1.35)
    )
  )

out <- file.path(OUTDIR, "fig9_rank_vs_nullz.png")
ggsave(out, fig9, width = 14, height = 6, dpi = 300, bg = "white")
message("Saved: ", out)