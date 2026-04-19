#!/usr/bin/env Rscript

library(tidyverse)
library(ggpubr)
library(rstatix)
library(patchwork)

# Load saved alpha objects
load("alpha_loop_results.RData")

# Consistent order across all panels
method_order <- c("benthic trawl",
                  "gillnet",
                  "midwater trawl",
                  "net",
                  "rod and reel",
                  "spear")

# -----------------------------
# Helper function
# -----------------------------
make_bw_alpha_plot <- function(df, body_site_name, method_order) {
  
# Clean data
  plot_df <- df %>%
    filter(!is.na(method_gear),
           method_gear != "not applicable",
           method_gear != "NA") %>%
    mutate(method_gear = factor(method_gear, levels = method_order)) %>%
    filter(!is.na(method_gear))
  
# Pairwise Wilcoxon (BH corrected)
pw_res <- plot_df %>%
    pairwise_wilcox_test(
      Shannon ~ method_gear,
      p.adjust.method = "BH"
    )
  
pw_sig <- pw_res %>%
    filter(p.adj < 0.05)
  
# Add bracket positions
  if (nrow(pw_sig) > 0) {
    pw_sig <- pw_sig %>%
      add_xy_position(x = "method_gear", fun = "max") %>%
      mutate(y.position = pmin(y.position, 4.8))  # keep inside 0–5
  }
  
  # Plot
  p <- ggplot(plot_df, aes(x = method_gear, y = Shannon)) +
    geom_boxplot(
      width = 0.65,
      fill = "white",
      color = "black",
      linewidth = 0.7,
      outlier.shape = 16,
      outlier.size = 2
    ) +
    scale_y_continuous(breaks = seq(0, 5, 1)) +
    coord_cartesian(ylim = c(0, 5), clip = "off") +
    labs(
      x = "Method gear",
      y = "Shannon index"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(hjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.line = element_line(linewidth = 0.7, color = "black"),
      axis.ticks = element_line(linewidth = 0.6, color = "black"),
      plot.margin = margin(10, 15, 10, 10)
    )
  
# Add significance brackets
  if (nrow(pw_sig) > 0) {
    p <- p +
      stat_pvalue_manual(
        pw_sig,
        label = "p.adj.signif",
        tip.length = 0.01,
        bracket.size = 0.4,
        size = 4,
        hide.ns = TRUE
      )
  }
  
  return(list(
    plot = p,
    pairwise_all = pw_res,
    pairwise_sig = pw_sig
  ))
}

# Create plots
gill_alpha <- make_bw_alpha_plot(gill_div_num, "Gill", method_order)
midgut_alpha <- make_bw_alpha_plot(midgut_div_num, "Midgut", method_order)
hindgut_alpha <- make_bw_alpha_plot(hindgut_div_num, "Hindgut", method_order)

# Show plots
gill_alpha$plot
midgut_alpha$plot
hindgut_alpha$plot

# Save figures

ggsave("Figure1a_gill_alpha_bw.png",
       gill_alpha$plot, width = 8, height = 4.5, dpi = 600, bg = "white")

ggsave("Figure1b_midgut_alpha_bw.png",
       midgut_alpha$plot, width = 8, height = 4.5, dpi = 600, bg = "white")

ggsave("Figure1c_hindgut_alpha_bw.png",
       hindgut_alpha$plot, width = 8, height = 4.5, dpi = 600, bg = "white")

# Save stats tables

write.csv(gill_alpha$pairwise_all, "gill_pairwise_all.csv", row.names = FALSE)
write.csv(gill_alpha$pairwise_sig, "gill_pairwise_significant.csv", row.names = FALSE)

write.csv(midgut_alpha$pairwise_all, "midgut_pairwise_all.csv", row.names = FALSE)
write.csv(midgut_alpha$pairwise_sig, "midgut_pairwise_significant.csv", row.names = FALSE)

write.csv(hindgut_alpha$pairwise_all, "hindgut_pairwise_all.csv", row.names = FALSE)
write.csv(hindgut_alpha$pairwise_sig, "hindgut_pairwise_significant.csv", row.names = FALSE)

