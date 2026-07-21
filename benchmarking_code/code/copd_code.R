# library(tidyverse)
library(scater)
# library(scran)
library(Seurat)
library(miloR)
library(dawnn)
library(SingleCellExperiment)
library(patchwork)

library(ggplot2)
library(plotly)
#library(ggpubr)

source("benchmarking_utilities.R")

# This code is adapted from Dann et al.: https://github.com/MarioniLab/milo_analysis_2020/blob/main/notebooks/Fig5_liver_cirrhosis.Rmd

options("future.globals.maxSize"=10**24)
cells <- readRDS("861b6b12-f9c9-4434-8d09-695a5156ce23_as_seurat.rds")

# If wanting to nornmalize cell numbers between case and control:
# sub_idxs <- sample(which(cells$disease == "normal"), sum(cells$disease != "normal"))
# cells <- cells[, c(which(cells$disease != "normal"), sub_idxs)]

# Dawnn
cells <- FindNeighbors(cells, return.neighbor = TRUE, k.param = 1001, reduction = "X_pca")

for (da_mode in c("ada", "pda")) {
    cells <- run_dawnn(cells, label_names = "disease", label_1 = "normal",
                       label_2 = "chronic obstructive pulmonary disease",
                       tf_conda_env = "tf_env", seed = 123, recalculate_graph = FALSE,
                       da_mode = da_mode, reduced_dim = "X_pca")
    cells@meta.data[, paste0("dawnn_", da_mode, "_lfc_rev")] <- -1 * cells@meta.data$dawnn_lfc
    cells@meta.data[, paste0("dawnn_", da_mode, "_verdict")] <- cells@meta.data$dawnn_da_verdict
}
dawnn_pda_out_df <- data.frame(lfc = cells$dawnn_pda_lfc_rev,
                               is_da = cells$dawnn_pda_verdict,
                               author_cell_type = cells$author_cell_type)
dawnn_ada_out_df <- data.frame(lfc = cells$dawnn_ada_lfc_rev,
                               is_da = cells$dawnn_ada_verdict,
                               author_cell_type = cells$author_cell_type)


# DA-seq
daseq_out <- getDAcells(cells@reductions[["X_pca"]]@cell.embeddings,
                        cell.labels = as.character(cells@meta.data$disease),
                        labels.1 = c("normal"), labels.2 = c("chronic obstructive pulmonary disease"),
                        k.vector = seq(50, 500, 50))
daseq_norm <- (daseq_out$da.pred / 2) + 0.5
daseq_out$lfc <- log2(daseq_norm / (1 - daseq_norm)) - log2(0.8355184/(1-0.8355184)) # is it correct to subtract here?
# If subsampled:
# daseq_out$lfc <- log2(daseq_norm / (1 - daseq_norm)) - log2(0.5/(1-0.5)) # is it correct to subtract here?
daseq_out$is_da <- seq(1, length(daseq_out$cell.idx)) %in% c(daseq_out$da.up, daseq_out$da.down)

daseq_out_df <- data.frame(lfc = daseq_out$lfc,
                           is_da = daseq_out$is_da,
                           author_cell_type = cells$author_cell_type)


# Milo
cells_milo <- Milo(as.SingleCellExperiment(cells))
cells_milo <- buildGraph(cells_milo, d = 11, k=30, reduced.dim = "X_PCA")
cells_milo <- makeNhoods(cells_milo, k=30, d=11, prop = 0.05, refined=TRUE, reduced_dims="X_PCA")

cells_meta <- as.tibble(colData(cells_milo)[,c("disease", "sample_id")])
cells_meta <- distinct(cells_meta) %>%
  mutate(disease=factor(disease, levels=c("normal", "chronic obstructive pulmonary disease"))) %>%
  column_to_rownames("sample_id")

cells_milo <- miloR::countCells(cells_milo, samples = "sample_id",
                         meta.data = data.frame(colData(cells_milo)[,c("sample_id","disease")]) )
cells_milo <- calcNhoodDistance(cells_milo, d = 11, reduced.dim = "X_PCA")
milo_res <- testNhoods(cells_milo, design = ~ disease,
                       design.df = cells_meta[colnames(nhoodCounts(cells_milo)), ,
                                              drop = FALSE],
                       reduced.dim = "X_PCA")

milo_res <- milo_res[,!str_detect(colnames(milo_res), "author_cell_type")]
milo_res <- annotateNhoods(cells_milo, milo_res, "author_cell_type")
anno_df <- data.frame(cells_milo@colData) %>%
  distinct(author_cell_type)
milo_res <- left_join(milo_res, anno_df, by="author_cell_type")
milo_res$author_cell_type[milo_res$author_cell_type_fraction < 0.6] <- NA

milo_out_df <- data.frame(lfc = milo_res$logFC,
                          is_da = milo_res$SpatialFDR < 0.1,
                          author_cell_type = milo_res$author_cell_type)


# Save all outputs
all_methods_df <- rbind(cbind(milo_out_df, method = "milo"),
                        cbind(daseq_out_df, method = "daseq"),
                        cbind(dawnn_pda_out_df, method = "dawnn_pda"),
                        cbind(dawnn_ada_out_df, method = "dawnn_ada"))
write.csv(all_methods_df, "cf_cells_results.csv")
