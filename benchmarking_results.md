---
title: "Analysis for 'Dawnn: single-cell differential abundance with neural networks'"
author: "George T. Hall"
date: "Compiled on 21 April 2023"
output:
    html_document:
        toc: true
        toc_depth: 2
        toc_float:
            collapsed: false
        number_sections: true
        code_folding: none
        keep_md: yes
---

<!---
(C) University College London 2023. This software is licenced under the terms
of the GNU GENERAL PUBLIC LICENSE Version 3. See COPYING.txt for the licence
details.
--->



This notebook generates the figures in the paper _Dawnn: single-cell
differential abundance with neural networks_. It is my intention that all
results shoud be completely reproducible. To this end, all objects and code
used to generate the results files read in this notebook are available
[here](this is the link). If you are unable to reproduce any result, please let
me know at `george.hall@ucl.ac.uk`.


```r
library(tidyverse)
library(ggplot2)
library(plotly)
library(ggpubr)
library(patchwork)
library(rstatix)
```

# Benchmarking TPR and FDR

We start by defining some functions to read results files and plot the figures.
For tidyness, these functions are hidden by default. Click on the title to
reveal the code.

<details>
<summary>`read_results` Read results file</summary>

```r
read_results <- function(out_file_name) {
    results <- read.csv(out_file_name, header = FALSE)
    colnames(results) <- c("method", "family", "tpr", "fpr", "fdr", "time",
                           "upreg_pc1", "upreg_pop", "rep")

    return(results)
}
```
</details>

<details>
<summary>`parse_results` Parse results dataframe</summary>

```r
parse_results <- function(results) {
    results <- tidyr::pivot_longer(results, 3:5, names_to = "stat")
    results$hline <- 0
    results[results$stat == "tpr", ]$hline <- 10 # Deliberately outside of graph
    results[results$stat == "fdr", ]$hline <- 0.1
    results$method <- factor(results$method,
                             levels = c("daseq", "milo", "dawnn"),
                             labels = c("DA-seq", "Milo", "Dawnn"))
    results$stat <- factor(results$stat, levels = c("tpr", "fpr", "fdr"),
                           labels = c("True Positive Rate",
                                      "False Positive Rate",
                                      "False Discovery Rate"))
    # Milo paper used natural logarithm here (i.e. log()), but I think log2 is
    # more natural (e.g. more common with DE analysis).
    lfcs <- log2(results$upreg_pc1 / (1 - results$upreg_pc1))
    results$upreg_lfc <- as.factor(round(lfcs, 1))

    results <- subset(results, subset = stat %in% c("True Positive Rate",
                                                    "False Discovery Rate"))
    results <- results %>% arrange(upreg_pop, rep)

    return(results)
}
```
</details>

<details>
<summary>`plot_accuracy_results` Plot accuracy results</summary>

```r
plot_accuracy_results <- function(results) {
    results <- parse_results(results)

    p <- ggplot(results, aes(x = upreg_lfc, y = value)) +
        geom_boxplot(aes(color = method), outlier.shape = NA) +
        geom_jitter(aes(color = method), size = 0.5, width = 0.2, alpha = 0.5) +
        geom_hline(aes(yintercept = hline, color = method), alpha = 0.75) +
        scale_color_manual(values = c("#1E88E5", "#FFC107", "#D81B60")) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none", axis.title.y = element_blank(),
              strip.background = element_blank()) +
        xlab(bquote(Maximum~log[2] * `-fold` ~ change)) + ylim(0, 1) +
        facet_grid(cols = vars(method), rows = vars(stat))

    return(p)
}
```
</details>

<details>
<summary>`plot_accuracy_results_signif` Plot accuracy results with statistics</summary>

```r
plot_accuracy_results_signif <- function(results) {

    # Built following
    # www.datanovia.com/en/blog/how-to-add-p-values-to-ggplot-facets

    results <- parse_results(results)
    bxp <- ggboxplot(results, x = "method", y = "value", color = "method",
                     add = "jitter", add.params = list(size = 0.2), size = 0.3,
                     facet.by = c("stat", "upreg_lfc")) +
           scale_color_manual(values = c("#1E88E5", "#FFC107", "#D81B60"))

    stat.test <- results %>%
                 group_by(stat, upreg_lfc) %>%
                 wilcox_test(value ~ method, p.adjust.method = "bonferroni",
                             exact = TRUE, paired = TRUE) %>%
                 add_xy_position(x = "method", dodge = 0.8)

    eff_size <- results %>%
                group_by(stat, upreg_lfc) %>%
                wilcox_effsize(value ~ method, exact = TRUE, paired = TRUE)
    stat.test.detailed <- results %>%
                          group_by(stat, upreg_lfc) %>%
                          wilcox_test(value ~ method, p.adjust.method =
                                      "bonferroni", exact = TRUE, paired =
                                      TRUE, detailed = TRUE) %>%
                          add_xy_position(x = "method", dodge = 0.8)

    eff_size_symbol <- ifelse(stat.test.detailed$estimate < 0, "<", ">")
    eff_size_vect <- c()
    for (i in seq_along(eff_size$effsize)) {
        eff_size_vect <- c(eff_size_vect,
                           paste0(rep(eff_size_symbol[i],
                                      round(4 * eff_size$effsize[i]) + 1),
                                  collapse = ""))
    }
    stat.test$eff_size <- eff_size_vect

    stat.test$y.position <- rep(c(1.120, 1.300, 1.480), 12)
    stat.test$y.position.effect_size <- rep(c(1.170, 1.350, 1.530), 12)

    p <- bxp +
         stat_pvalue_manual(stat.test, label = "p.adj.signif",
                            tip.length = 0.01, hide.ns = FALSE) +
         stat_pvalue_manual(stat.test, label = "eff_size",
                            tip.length = 0.01, hide.ns = TRUE,
                            y.position = "y.position.effect_size",
                            remove.bracket = TRUE) +
         scale_y_continuous(breaks = seq(0, 1, 0.2),
                            expand = expansion(mult = c(0.01, 0.1))) &
         theme(axis.text.x = element_text(angle = 45, hjust = 1),
               legend.position = "none", axis.title = element_blank(),
               strip.background = element_blank())

    return(p)
}
```
</details>

<details>
<summary>`ggsave_plots` Save ggplots</summary>

```r
ggsave_plots <- function(non_stat_plot, stat_plot, file_name) {
    ggsave(filename = paste0(file_name, ".png", collapse = ""),
           plot = non_stat_plot, dpi = 300, units = "cm", width = 12,
           height = 10)
    ggsave(filename = paste0(file_name, "_stat.png", collapse = ""),
           plot = stat_plot, dpi = 300, units = "cm", width = 12,
           height = 10)
}
```
</details>

<br>

We now plot the accuracy results (with and without tests for statistical
significance) for each benchmarking dataset.

## Simulated discrete clusters


```r
# Generated as in data_for_benchmarking_results/collect_results_all_sim_dat.R
results <- read_results("data_for_benchmarking_results/tpr_fdr_results_discrete_clusters_rerun.csv")
(p_disc_clusts_acc <- plot_accuracy_results(results))
```

<img src="benchmarking_results_files/figure-html/discrete_clusters_acc-1.png" style="display: block; margin: auto;" />

```r
(p_disc_clusts_acc_stats <- plot_accuracy_results_signif(results))
```

<img src="benchmarking_results_files/figure-html/discrete_clusters_acc-2.png" style="display: block; margin: auto;" />

```r
ggsave_plots(p_disc_clusts_acc, p_disc_clusts_acc_stats,
             "manuscript/images/sim_discrete_clusters_benchmarking")
```

## Simulated linear trajectory


```r
# Generated as in data_for_benchmarking_results/collect_results_all_sim_dat.R
results <- read_results("data_for_benchmarking_results/tpr_fdr_results_linear_traj_rerun.csv")
(p_linear_traj_acc <- plot_accuracy_results(results))
```

<img src="benchmarking_results_files/figure-html/linear_traj_acc-1.png" style="display: block; margin: auto;" />

```r
(p_linear_traj_acc_stats <- plot_accuracy_results_signif(results))
```

<img src="benchmarking_results_files/figure-html/linear_traj_acc-2.png" style="display: block; margin: auto;" />

```r
ggsave_plots(p_linear_traj_acc, p_linear_traj_acc_stats,
             "manuscript/images/sim_linear_traj_benchmarking")
```

## Simulated branching trajectory


```r
# Generated as in data_for_benchmarking_results/collect_results_all_sim_dat.R
results <- read_results("data_for_benchmarking_results/tpr_fdr_results_branch_traj_rerun.csv")

# The effect size code was crashing here since, at upreg_pc1 = 0.7, both Milo
# and Dawnn had only fdrs of 0, which led to all differences being 0. So below
# I artificially increase one fdr datapoint for each method by 0.00000000001,
# to introduce a non-zero difference. This increase is not reflected in the
# plots.

milo_upreg_pc1_07_idxs <- intersect(which(results$upreg_pc1 == 0.7),
                                    which(results$method == "milo"))
results[milo_upreg_pc1_07_idxs[1], ]$fdr <- 0.00000000001
daseq_upreg_pc1_07_idxs <- intersect(which(results$upreg_pc1 == 0.7),
                                     which(results$method == "daseq"))
results[daseq_upreg_pc1_07_idxs[2], ]$fdr <- 0.00000000001

(p_branch_traj_acc <- plot_accuracy_results(results))
```

<img src="benchmarking_results_files/figure-html/branch_traj_acc-1.png" style="display: block; margin: auto;" />

```r
(p_branch_traj_acc_stats <- plot_accuracy_results_signif(results))
```

<img src="benchmarking_results_files/figure-html/branch_traj_acc-2.png" style="display: block; margin: auto;" />

```r
ggsave_plots(p_branch_traj_acc, p_branch_traj_acc_stats,
             "manuscript/images/sim_branching_traj_benchmarking")
```

## Mouse gastrulation dataset


```r
# Generated as in data_for_benchmarking_results/collecting_results_mouse.sh
results <- read_results("data_for_benchmarking_results/tpr_fdr_results_mouse_regen.csv")
(p_mouse_acc <- plot_accuracy_results(results))
```

<img src="benchmarking_results_files/figure-html/mouse_acc-1.png" style="display: block; margin: auto;" />

```r
(p_mouse_acc_stats <- plot_accuracy_results_signif(results))
```

<img src="benchmarking_results_files/figure-html/mouse_acc-2.png" style="display: block; margin: auto;" />

```r
ggsave_plots(p_mouse_acc, p_mouse_acc_stats,
             "manuscript/images/mouse_gastrulation_benchmarking")
```

## Skin dataset


```r
# Generated as in data_for_benchmarking_results/collecting_results_skin.sh
results <- read_results("data_for_benchmarking_results/tpr_fdr_results_skin_regen.csv")
(p_skin_acc <- plot_accuracy_results(results))
```

<img src="benchmarking_results_files/figure-html/skin_acc-1.png" style="display: block; margin: auto;" />

```r
(p_skin_acc_stats <- plot_accuracy_results_signif(results))
```

<img src="benchmarking_results_files/figure-html/skin_acc-2.png" style="display: block; margin: auto;" />

```r
ggsave_plots(p_skin_acc, p_skin_acc_stats,
             "manuscript/images/keratinocyte_benchmarking")
```

## Organoid dataset


```r
# Generated as in data_for_benchmarking_results/collecting_results_organoid.sh
results <- read_results("data_for_benchmarking_results/tpr_fdr_results_organoid_regen.csv")
(p_organoid_acc <- plot_accuracy_results(results))
```

<img src="benchmarking_results_files/figure-html/organoid_acc-1.png" style="display: block; margin: auto;" />

```r
(p_organoid_acc_stats <- plot_accuracy_results_signif(results))
```

<img src="benchmarking_results_files/figure-html/organoid_acc-2.png" style="display: block; margin: auto;" />

```r
ggsave_plots(p_organoid_acc, p_organoid_acc_stats,
             "manuscript/images/bile_duct_organoids_benchmarking")
```

## Heart dataset


```r
# Generated as in data_for_benchmarking_results/collecting_results_heart.sh
results <- read_results("data_for_benchmarking_results/tpr_fdr_results_heart_regen.csv")

(p_heart_acc <- plot_accuracy_results(results))
```

<img src="benchmarking_results_files/figure-html/heart_acc-1.png" style="display: block; margin: auto;" />

```r
(p_heart_acc_stats <- plot_accuracy_results_signif(results))
```

<img src="benchmarking_results_files/figure-html/heart_acc-2.png" style="display: block; margin: auto;" />

```r
ggsave_plots(p_heart_acc, p_heart_acc_stats,
             "manuscript/images/heart_benchmarking")
```

# Recovering liver cirrhosis findings

We follow Dann _et al._ in reproducing published results using Dawnn (see
[Figure 5d](https://www.nature.com/articles/s41587-021-01033-z/figures/5) in
their paper). This section is a very lightly adapted version of [their
code](https://github.com/MarioniLab/milo_analysis_2020/blob/main/notebooks/Fig5_liver_cirrhosis.Rmd#L261).
We demonstrate that DA-seq is also able to recover these published findings.


```r
# Label cell types as DA or not, as per dataset's paper
group.by <- "annotation_indepth"
paper_DA <- list(cirrhotic = c("MPs (4)", "MPs (5)", "Endothelia (6)",
                               "Endothelia (7)", "Mes (3)", "Tcells (2)",
                               "Myofibroblasts"),
                 healthy = c("MPs (7)", "Endothelia (1)", "Tcells (1)",
                             "Tcells (3)", "Tcells (1)", "ILCs (1)"))
expDA_df <- bind_rows(
  data.frame(annotation_indepth = paper_DA[["cirrhotic"]],
             pred_DA = "cirrhotic"),
  data.frame(annotation_indepth = paper_DA[["healthy"]],
             pred_DA = "healthy")
  )

create_beeswarm_plot <- function(res, expDA_df) {
    p <- res %>%
      left_join(expDA_df) %>%
      mutate(logFC_color = ifelse(is_da == 1, lfc, NA)) %>%
      arrange(annotation_lineage) %>%
      filter(!is.na(annotation_lineage)) %>%
      ggplot(aes(annotation_indepth, lfc, color = logFC_color)) +
      scale_color_gradient2() +
      guides(color = "none") +
      scale_y_continuous(breaks = scales::pretty_breaks(n = 14)) +
      xlab(group.by) + ylab("Log Fold Change") +
      ggbeeswarm::geom_quasirandom(alpha = 1, size = 0.2) +
      coord_flip() +
      facet_grid(annotation_lineage~., scales = "free", space = "free") +
      theme_bw(base_size = 22) +
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "pt")) +
      ylim(-6.5, 6.5)
}
```


```r
# Results generated by "benchmarking_liver_cirrhosis_analysis.R".
liver_cirrhosis_results <- read.csv("data_for_benchmarking_results/liver_cirrhosis_results_rerun.csv")
milo_cirrhosis_df <- subset(liver_cirrhosis_results, method == "milo")
daseq_cirrhosis_df <- subset(liver_cirrhosis_results, method == "daseq")
dawnn_cirrhosis_df <- subset(liver_cirrhosis_results, method == "dawnn")

# Subset out only the cell types also caught by Milo
cell_types_in_milo <- unique(milo_cirrhosis_df$annotation_indepth)
daseq_cirrhosis_df <- subset(daseq_cirrhosis_df,
                             annotation_indepth %in% cell_types_in_milo)
dawnn_cirrhosis_df <- subset(dawnn_cirrhosis_df,
                             annotation_indepth %in% cell_types_in_milo)

# Generate Milo plot
p_milo <- create_beeswarm_plot(milo_cirrhosis_df, expDA_df)
p_milo <- p_milo + theme(strip.text.y = element_blank(),
                         axis.title.x = element_blank()) + ggtitle("Milo")

# Generate DAseq plot
p_daseq <- create_beeswarm_plot(daseq_cirrhosis_df, expDA_df)
p_daseq <- p_daseq + theme(strip.text.y = element_blank()) + ggtitle("DA-seq")

# Generate Dawnn plot
p_dawnn <- create_beeswarm_plot(dawnn_cirrhosis_df, expDA_df)
p_dawnn <- p_dawnn + theme(strip.text.y = element_text(angle = 0),
                           axis.title.x = element_blank()) + ggtitle("Dawnn")

lhs <- dawnn_cirrhosis_df %>%
  left_join(expDA_df) %>%
  group_by(annotation_indepth) %>%
  summarise(pred_DA = dplyr::first(pred_DA),
            annotation_lineage = dplyr::first(annotation_lineage)) %>%
  mutate(end = ifelse(pred_DA == "healthy", 0, 1),
         start = ifelse(pred_DA == "healthy", 1, 0)) %>%
  filter(!is.na(annotation_lineage)) %>%
  ggplot(aes(annotation_indepth, start, xend = annotation_indepth, yend = end,
             color = pred_DA)) +
  geom_segment(size = 1,
               arrow = arrow(length = unit(0.1, "npc"), type = "closed")) +
  coord_flip() +
  xlab("Cell type annotation") +
  facet_grid(annotation_lineage~., scales = "free", space = "free") +
  scale_color_brewer(palette = "Set1", direction = -1,
                     labels = c("enriched in cirrhotic", "enriched in healthy"),
                     na.translate = FALSE,
                     name = "Ramachandran et al.\nDA predictions") +
  guides(color = guide_legend(ncol = 1)) +
  theme_bw(base_size = 22) +
  ylim(-0.1, 1.1) +
  theme(strip.text.y = element_blank(),
        strip.text.x = element_text(angle = 90),
        plot.margin = unit(c(0, 0, 0, 0), "cm"), panel.grid = element_blank(),
        axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), legend.position = "bottom")

p_all_beeswarms <- lhs + p_milo + p_daseq + p_dawnn +
     plot_layout(widths = c(1, 10, 10, 10), guides = "collect") &
     theme(legend.position = "bottom", legend.justification = 0,
           plot.title = element_text(hjust = 0.5))
p_all_beeswarms
```

<img src="benchmarking_results_files/figure-html/generate_liver_beeswarms-1.png" style="display: block; margin: auto;" />

```r
ggsave(filename = "manuscript/images/cirrhosis_all_methods_beeswarms.png",
       plot = p_all_beeswarms, dpi = 300, units = "cm", width = 30, height = 40)
```

# P(Condition1) estimation accuracy

To compare against MELD - which does not call cells as in regions of
differential abundance - we again follow Dann _et al._ and compare the mean
squared errors of their estimates of P(Condition1). We also do so for DA-seq.


```r
plot_pc1_estimates <- function(results, other_method, other_mses) {
    mses_df <- data.frame(x = 0.8, y = 0.35, lab = other_mses,
                          cell_type = factor(sort(rep(unique(results$cell_type), 3))),
                          upreg_pc1 = factor(rep(unique(results$upreg_pc1), 3)))
    min_x <- min(results$truth)
    max_x <- max(results$truth)
    min_y <- min(c(results$meld, results$dawnn, results$daseq))
    max_y <- max(c(results$meld, results$dawnn, results$daseq))

    p <- ggplot(results, aes(x = truth, y = .data[[other_method]])) +
            geom_point(size = 0, alpha = 0.03) +
            xlab("Ground truth") + ylab("Estimated value") +
            xlim(c(min_x, max_x)) + ylim(c(min_y, max_y)) +
            geom_abline(intercept = 0, slope = 1, color = "red") +
            geom_text(data = mses_df, aes(x = x, y = y, label = lab)) +
            facet_grid(upreg_pc1 ~ cell_type) & theme_bw() &
            theme(strip.background = element_blank(),
                  plot.title = element_text(hjust = 0.5))

    return(p)
}

meld_out <- data.frame(data.table::fread("data_for_benchmarking_results/meld_pc1_ests.csv",
                                         header = FALSE))
meld_out <- meld_out[((3 * (1:9)) - 2), ]
meld_out <- cbind(meld_out, upreg_pc1 = c(rep(c(0.7, 0.8, 0.9), 3)))
meld_out <- cbind(meld_out, cell_type = c(rep("Erythroid1", 3),
                                          rep("Gut", 3),
                                          rep("Somitic mesoderm", 3)))
meld_out_vector <- as.numeric(c(t(as.matrix(meld_out[, 2:64019]))))

dawnn_out <- data.frame(data.table::fread("data_for_benchmarking_results/dawnn_mse_benchmarking_regen_final_model_dawnn_rerun.csv",
                                          header = FALSE))
dawnn_out_vector <- as.numeric(c(t(as.matrix(dawnn_out[, 1:64018]))))

daseq_out <- data.frame(data.table::fread("data_for_benchmarking_results/daseq_pc1_ests.csv",
                                          header = FALSE))
daseq_out_vector <- as.numeric(c(t(as.matrix(daseq_out[, 1:64018]))))

ground_truth <- data.table::fread("data_for_benchmarking_results/pc1_ground_truth_regen.csv",
                                  header = FALSE)
ground_truth <- ground_truth[order(ground_truth$V2), ]
ground_truth <- ground_truth[, 5:(4 + 64018)]
ground_truth_vector <- as.numeric(c(as.matrix(t(ground_truth))))

cell_type_vector <- c(sapply(meld_out$cell_type, function(x) {rep(x, 64018)}))
upreg_pc1_vector <- c(sapply(meld_out$upreg_pc1, function(x) {rep(x, 64018)}))

results <- data.frame(cell_name = rep(paste0("cell_", 1:64018), 9),
                      truth = ground_truth_vector, meld = meld_out_vector,
                      dawnn = dawnn_out_vector,
                      daseq = daseq_out_vector,
                      cell_type = cell_type_vector,
                      upreg_pc1 = upreg_pc1_vector)

results$meld_se <- (results$truth - results$meld)**2
results$dawnn_se <- (results$truth - results$dawnn)**2
results$daseq_se <- (results$truth - results$daseq)**2

results$cell_type_and_upreg_pc1 <- paste(results$cell_type, results$upreg_pc1,
                                         sep = "_")
meld_mses <- aggregate(results$meld_se, list(results$cell_type_and_upreg_pc1),
                       FUN = mean)
dawnn_mses <- aggregate(results$dawnn_se, list(results$cell_type_and_upreg_pc1),
                       FUN = mean)
daseq_mses <- aggregate(results$daseq_se, list(results$cell_type_and_upreg_pc1),
                       FUN = mean)

p_meld <- plot_pc1_estimates(results, "meld", round(meld_mses$x, 5))
ggsave(filename = "manuscript/images/mses_meld.png", plot = p_meld, dpi = 300,
       units = "cm", width = 12, height = 10)
p_dawnn <- plot_pc1_estimates(results, "dawnn", round(dawnn_mses$x, 5))
ggsave(filename = "manuscript/images/mses_dawnn.png", plot = p_dawnn, dpi = 300,
       units = "cm", width = 12, height = 10)
p_daseq <- plot_pc1_estimates(results, "daseq", round(daseq_mses$x, 5))
ggsave(filename = "manuscript/images/mses_daseq.png", plot = p_daseq, dpi = 300,
       units = "cm", width = 12, height = 10)

layout <- "
#AA#
BBCC
"
p_mses <- (p_dawnn + ggtitle("Dawnn")) + (p_meld + ggtitle("MELD")) +
          (p_daseq + ggtitle("DA-seq")) + plot_layout(design = layout)
p_mses
```

<img src="benchmarking_results_files/figure-html/pc1_mses-1.png" style="display: block; margin: auto;" />

```r
ggsave(filename = "manuscript/images/mse_comparisons.png", plot = p_mses,
       dpi = 300, units = "cm", width = 25, height = 20)
```

# Runtime analysis

Finally, we compare the runtimes of Dawnn, Milo, and DA-seq.


```r
runtime_results <- read.csv("data_for_benchmarking_results/runtime_results_regen.csv",
                            header = FALSE)
colnames(runtime_results) <- c("Method", "n", "t", "rep")
runtime_results$Method[runtime_results$Method == "dawnn"] <- "Dawnn"
runtime_results$Method[runtime_results$Method == "milo"] <- "Milo"
runtime_results$Method[runtime_results$Method == "daseq"] <- "DA-seq"
runtime_results$Method <- factor(runtime_results$Method,
                                 levels = c("DA-seq", "Milo", "Dawnn"))

runtime_plot <- ggplot(runtime_results, aes(x = n, y = t, color = Method)) +
                geom_point() + xlab("Number of cells") + ylab("Time (secs)") +
                scale_color_manual(values = c("#1E88E5", "#FFC107", "#D81B60")) +
                theme_classic() + theme(legend.title = element_blank())
runtime_plot
```

<img src="benchmarking_results_files/figure-html/runtime_results-1.png" style="display: block; margin: auto;" />

```r
ggsave(filename = "manuscript/images/runtime_plot.pdf", plot = runtime_plot,
       dpi = 300, units = "cm", width = 20, height = 10)
```

# Session info


```r
sessionInfo()
```

```
## R version 4.0.3 (2020-10-10)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Big Sur 10.16
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] rstatix_0.7.0        patchwork_1.1.0.9000 ggpubr_0.4.0        
##  [4] plotly_4.10.0        forcats_0.5.1        stringr_1.4.0       
##  [7] dplyr_1.1.0          purrr_0.3.4          readr_2.1.3         
## [10] tidyr_1.2.0          tibble_3.1.8         ggplot2_3.4.0       
## [13] tidyverse_1.3.1      rmarkdown_2.9       
## 
## loaded via a namespace (and not attached):
##  [1] matrixStats_0.61.0 fs_1.5.2           lubridate_1.8.0    RColorBrewer_1.1-2
##  [5] httr_1.4.2         tools_4.0.3        backports_1.2.1    bslib_0.3.1       
##  [9] utf8_1.2.2         R6_2.5.1           vipor_0.4.5        DBI_1.1.1         
## [13] lazyeval_0.2.2     colorspace_2.0-2   withr_2.5.0        tidyselect_1.2.0  
## [17] curl_4.3.2         compiler_4.0.3     textshaping_0.3.1  cli_3.6.0         
## [21] rvest_1.0.2        xml2_1.3.3         sandwich_3.0-1     labeling_0.4.2    
## [25] sass_0.4.0         scales_1.2.1       mvtnorm_1.1-1      systemfonts_1.0.1 
## [29] digest_0.6.29      foreign_0.8-81     rio_0.5.27         pkgconfig_2.0.3   
## [33] htmltools_0.5.2    highr_0.8          dbplyr_2.1.1       fastmap_1.1.0     
## [37] htmlwidgets_1.5.4  rlang_1.0.6        readxl_1.3.1       rstudioapi_0.13   
## [41] jquerylib_0.1.4    generics_0.1.3     farver_2.1.0       zoo_1.8-9         
## [45] jsonlite_1.7.2     zip_2.1.1          car_3.0-11         magrittr_2.0.1    
## [49] modeltools_0.2-23  Matrix_1.5-3       ggbeeswarm_0.6.0   Rcpp_1.0.7        
## [53] munsell_0.5.0      fansi_0.5.0        abind_1.4-5        lifecycle_1.0.3   
## [57] multcomp_1.4-20    stringi_1.7.6      yaml_2.2.1         carData_3.0-4     
## [61] MASS_7.3-53.1      grid_4.0.3         parallel_4.0.3     crayon_1.4.2      
## [65] lattice_0.20-41    haven_2.5.1        splines_4.0.3      hms_1.1.2         
## [69] knitr_1.31         pillar_1.8.1       ggsignif_0.6.2     codetools_0.2-18  
## [73] stats4_4.0.3       reprex_2.0.1       glue_1.6.2         evaluate_0.14     
## [77] data.table_1.14.2  modelr_0.1.8       vctrs_0.5.2        tzdb_0.3.0        
## [81] cellranger_1.1.0   gtable_0.3.0       assertthat_0.2.1   xfun_0.24         
## [85] openxlsx_4.2.5     coin_1.4-2         libcoin_1.0-9      broom_1.0.3       
## [89] ragg_1.1.1         survival_3.2-7     viridisLite_0.4.0  beeswarm_0.3.1    
## [93] TH.data_1.1-1      ellipsis_0.3.2
```
