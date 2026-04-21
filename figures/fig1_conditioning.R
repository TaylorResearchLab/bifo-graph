# =============================================================================
# Figure 1: BIFO conditioning coverage — edge funnel and flow class distribution
#
# Panels:
#   A: Edge funnel — input → kept → propagating (stacked bar / waterfall)
#   B: Dropped edge breakdown by reason
#   C: Flow class distribution of conditioned propagating edges
#   D: Node resolution coverage
#
# Data: results/chd_curated/results_full.json (coverage + graph_stats blocks)
# Output: figures/fig1_conditioning.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(jsonlite)

BASE    <- "/mnt/isilon/taylor_lab/data/projects/BIFO_2026"
CURATED <- file.path(BASE, "results/chd_curated")
OUTDIR  <- file.path(BASE, "figures")
dir.create(OUTDIR, showWarnings = FALSE)

fix_json <- function(path) {
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  txt <- gsub("\\bNaN\\b",       "null", txt)
  txt <- gsub("\\bInfinity\\b",  "null", txt)
  txt <- gsub("\\b-Infinity\\b", "null", txt)
  fromJSON(txt)
}

rf  <- fix_json(file.path(CURATED, "results_full.json"))
cov <- rf$coverage
gs  <- rf$graph_stats

# ── Panel A: Edge funnel waterfall ────────────────────────────────────────────
funnel_df <- tibble(
  stage = factor(c(
    "Input edges\n(merged)",
    "Resolved\nendpoints",
    "Kept after\nconditioning",
    "Propagating\nedges"
  ), levels = c(
    "Input edges\n(merged)",
    "Resolved\nendpoints",
    "Kept after\nconditioning",
    "Propagating\nedges"
  )),
  count = c(
    cov$total_concept_edges,
    cov$concept_edges_with_resolved_endpoints,
    cov$kept_edges,
    gs$conditioned_propagating_edges
  ),
  pct = c(100,
          100 * cov$concept_edges_with_resolved_endpoints / cov$total_concept_edges,
          100 * cov$kept_edges / cov$total_concept_edges,
          100 * gs$conditioned_propagating_edges / cov$total_concept_edges)
)

pA <- ggplot(funnel_df, aes(x = stage, y = count)) +
  geom_col(fill = "#2E75B6", width = 0.55, alpha = 0.85) +
  geom_text(aes(label = paste0(comma(count), "\n(", sprintf("%.1f", pct), "%)")),
            vjust = -0.3, size = 3.2, fontface = "bold", lineheight = 0.9) +
  scale_y_continuous(labels = comma,
                     limits = c(0, cov$total_concept_edges * 1.18)) +
  labs(title    = "A",
       subtitle = "Edge funnel — input to propagating",
       x = NULL, y = "Edge count") +
  theme_classic(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        axis.text.x   = element_text(size = 9, angle = 15, hjust = 0.9))

# ── Panel B: Dropped edge breakdown ───────────────────────────────────────────
drop_df <- tibble(
  reason = factor(c(
    "Unmapped\npredicate",
    "Unresolved\nentity",
    "Non-flow\nclassification"
  ), levels = c(
    "Unmapped\npredicate",
    "Unresolved\nentity",
    "Non-flow\nclassification"
  )),
  count = c(
    cov$dropped_unmapped_predicate,
    cov$dropped_unresolved_entity,
    cov$dropped_non_flow
  )
)
drop_df <- drop_df %>%
  mutate(pct = 100 * count / cov$total_concept_edges)

pB <- ggplot(drop_df, aes(x = reason, y = count, fill = reason)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.4) +
  geom_text(aes(label = paste0(comma(count), "\n(", sprintf("%.1f", pct), "%)")),
            vjust = -0.3, size = 3.2, fontface = "bold", lineheight = 0.9) +
  scale_fill_manual(values = c(
    "Unmapped\npredicate"      = "#ED7D31",
    "Unresolved\nentity"       = "#FFC000",
    "Non-flow\nclassification" = "#A9A9A9"
  )) +
  scale_y_continuous(labels = comma,
                     limits = c(0, max(drop_df$count) * 1.25)) +
  labs(title    = "B",
       subtitle = "Dropped edges by reason",
       x = NULL, y = "Edge count") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        axis.text.x   = element_text(size = 9, angle = 15, hjust = 0.9))

# ── Panel C: Flow class distribution ─────────────────────────────────────────
flow_df <- tibble(
  flow_class = names(gs$flow_class_distribution),
  count      = unlist(gs$flow_class_distribution)
) %>%
  arrange(desc(count)) %>%
  mutate(
    flow_class = factor(flow_class, levels = rev(flow_class)),
    pct        = 100 * count / gs$conditioned_propagating_edges,
    is_bridge  = grepl("Pathway Contribution", flow_class)
  )

pC <- ggplot(flow_df, aes(x = count, y = flow_class, fill = is_bridge)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0(comma(count), " (", sprintf("%.1f", pct), "%)")),
            hjust = -0.08, size = 3.0, fontface = "bold") +
  scale_fill_manual(values = c("FALSE" = "#6BAED6", "TRUE" = "#2E75B6")) +
  scale_x_continuous(labels = comma,
                     limits = c(0, max(flow_df$count) * 1.35)) +
  labs(title    = "C",
       subtitle = "Flow class distribution — conditioned propagating edges\n(blue = Pathway Contribution bridge edges)",
       x = "Edge count", y = NULL) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        axis.text.y   = element_text(size = 9.5))

# ── Panel D: Node resolution ──────────────────────────────────────────────────
node_df <- tibble(
  category = factor(c("Resolved", "Unresolved"),
                    levels = c("Resolved", "Unresolved")),
  count    = c(cov$resolved_concept_nodes, cov$unresolved_concept_nodes),
  pct      = 100 * c(cov$resolved_concept_nodes,
                      cov$unresolved_concept_nodes) / cov$total_concept_nodes
)

pD <- ggplot(node_df, aes(x = category, y = count, fill = category)) +
  geom_col(width = 0.5, color = "white", linewidth = 0.4) +
  geom_text(aes(label = paste0(comma(count), "\n(", sprintf("%.1f", pct), "%)")),
            vjust = -0.3, size = 3.5, fontface = "bold", lineheight = 0.9) +
  scale_fill_manual(values = c("Resolved" = "#2E75B6", "Unresolved" = "#C8C8C8")) +
  scale_y_continuous(labels = comma,
                     limits = c(0, cov$total_concept_nodes * 1.18)) +
  labs(title    = "D",
       subtitle = paste0("Node entity resolution\n(total nodes: ",
                         comma(cov$total_concept_nodes), ")"),
       x = NULL, y = "Node count") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none",
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40"))

# ── Assemble ──────────────────────────────────────────────────────────────────
fig1 <- ((pA | pB | pD) / pC) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title    = "Figure 1. BIFO conditioning coverage — curated CHD benchmark graph",
    subtitle = paste0(
      "Input: ", comma(cov$total_concept_edges),
      " merged edges (94,790 seed\u2194hop1 + 79,562 pathway membership).\n",
      "Level 1 entity resolution: ", sprintf("%.1f", 100 * cov$level1_concept_coverage),
      "% of concept nodes resolved to source vocabulary (SAB).\n",
      "Level 2 edge conditioning: ", sprintf("%.1f", 100 * cov$level2_edge_coverage),
      "% of edges kept as biologically admissible; ",
      comma(gs$conditioned_propagating_edges), " propagating edges used in PPR operator.\n",
      "Pathway Contribution bridge edges (", comma(gs$flow_class_distribution[["Pathway Contribution"]]),
      "; ", sprintf("%.1f", 100 * gs$flow_class_distribution[["Pathway Contribution"]] /
                      gs$conditioned_propagating_edges), "% of propagating edges) ",
      "are the architectural finding enabling gene\u2192pathway signal transfer."
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 8.5, color = "grey35", lineheight = 1.35)
    )
  )

out <- file.path(OUTDIR, "fig1_conditioning.png")
ggsave(out, fig1, width = 13, height = 10, dpi = 300, bg = "white")
message("Saved: ", out)
