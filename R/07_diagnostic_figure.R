# 01. Reset --------------------------------------------------------------------
rm(list = ls())

# 02. Libraries ------------------------------------------------------------------
library("easypackages")
libraries('ggplot2', 'patchwork', 'dplyr')

# 03. What this script does -----------------------------------------------------
# Builds a two-panel diagnostic figure summarising 06's two-tier validation.
# (A) Predicted vs. true simulated breeding value, both methods, the same 315
# genotypes: the artifact-free check that the pipeline recovers genetic signal,
# possible only because the population is simulated. (B) The bootstrap
# distribution of the accuracy difference for Table 5's literal operational gate
# (vs. the observed composite, the only outcome a real programme has), with the
# 95% CI and the zero line. Panel B's subtitle carries the MGIDI in-sample
# caveat from 06's header note.

IN  <- "../outputs"
OUT <- "../figures"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

gate <- readRDS(file.path(IN, "stage_gate_decision.rds"))
cmp  <- gate$comparison_table

# 04. Panel A: predicted vs. TRUE breeding value, both methods (oracle validation) ---
panel_data <- bind_rows(
  cmp %>% transmute(Composite_TrueBV, Predicted = Composite_GEBV, Method = "Stage 3: GBLUP (composite, out-of-fold)"),
  cmp %>% transmute(Composite_TrueBV, Predicted = -MGIDI, Method = "Stage 1: MGIDI (negated, whole-population)")
)

acc_labels <- data.frame(
  Method = c("Stage 3: GBLUP (composite, out-of-fold)", "Stage 1: MGIDI (negated, whole-population)"),
  label  = c(sprintf("r = %.3f", gate$acc_gblup_oracle), sprintf("r = %.3f", gate$acc_stage1_oracle))
)

panel_a <- ggplot(panel_data, aes(Composite_TrueBV, Predicted)) +
  geom_point(alpha = 0.5, size = 1.4, color = "#2C6E7F") +
  geom_smooth(method = "lm", se = FALSE, color = "#B24C36", linewidth = 0.6) +
  geom_text(data = acc_labels, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.15, vjust = 1.5, size = 3.2, inherit.aes = FALSE) +
  facet_wrap(~ Method, ncol = 1, scales = "free_y") +
  labs(x = "TRUE composite breeding value (simulation ground truth)",
       y = "Predicted merit (composite GEBV or negated MGIDI)",
       title = "A. Oracle validation: predicted vs. TRUE breeding value",
       subtitle = sprintf("Noise ceiling (observed composite vs. true value): r = %.3f", gate$noise_ceiling)) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold", size = 8.5),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey30"))

# 05. Panel B: bootstrap distribution for Table 5's literal operational gate -----
boot_df <- data.frame(diff = gate$boot_diff)
ci <- gate$ci
verdict_short <- sub("^(GO|NO-GO).*$", "\\1", gate$verdict)  # drop the long parenthetical for title width

panel_b <- ggplot(boot_df, aes(diff)) +
  geom_histogram(bins = 60, fill = "#2C6E7F", alpha = 0.75, color = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_vline(xintercept = ci, color = "#B24C36", linewidth = 0.6) +
  annotate("text", x = 0, y = Inf, label = "zero", vjust = 1.4, hjust = -0.15,
           size = 3, color = "grey30") +
  labs(x = "Bootstrap (GBLUP accuracy − Stage-1 accuracy), vs. observed composite",
       y = sprintf("Count (of %s resamples)", format(gate$n_boot, big.mark = ",")),
       title = sprintf("B. Table 5's literal operational gate: %s", verdict_short),
       subtitle = sprintf("%.0f%% CI: [%.3f, %.3f]; MGIDI not held out here (see README)",
                           gate$ci_level * 100, ci[1], ci[2])) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "grey30"))

# 06. Combine and write -----------------------------------------------------------
fig <- panel_a / panel_b + plot_layout(heights = c(1.4, 1))

ggsave(file.path(OUT, "Fig_stage3_gate_diagnostic.png"), fig,
       width = 6.5, height = 8.5, dpi = 300, bg = "white")
ggsave(file.path(OUT, "Fig_stage3_gate_diagnostic.pdf"), fig,
       width = 6.5, height = 8.5, bg = "white")

cat("Wrote:", file.path(OUT, "Fig_stage3_gate_diagnostic.png"), "\n")
