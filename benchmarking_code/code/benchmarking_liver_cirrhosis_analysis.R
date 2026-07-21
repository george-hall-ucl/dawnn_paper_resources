library(tidyverse)
library(scater)
library(scran)
library(Seurat)
library(miloR)
library(dawnn)
library(SingleCellExperiment)
library(patchwork)

library(ggplot2)
library(plotly)
#library(ggpubr)

# This code is adapted from Dann et al.: https://github.com/MarioniLab/milo_analysis_2020/blob/main/notebooks/Fig5_liver_cirrhosis.Rmd

source("benchmarking_utilities.R")

# Download from https://datashare.ed.ac.uk/bitstream/handle/10283/3433/tissue.rdata?sequence=3&isAllowed=y
# Save as "liver_cells.rdata"
load("../data/input_data/liver_cells.rdata")

liver_sce <- SingleCellExperiment(assay = list(counts = tissue@raw.data,
                                               logcounts = tissue@data),
                                  colData = tissue@meta.data)
dec_liver <- modelGeneVar(liver_sce)
fit_liver <- metadata(dec_liver)

hvgs <- getTopHVGs(dec_liver, n = 3000)
liver_sce <- runPCA(liver_sce, subset_row = hvgs, ncomponents = 11)
liver_sce <- runUMAP(liver_sce, dimred = "PCA", ncomponents = 2)
liver_seu <- as.Seurat(liver_sce)

liver_seu@meta.data$condition <- as.character(liver_seu@meta.data$condition)
liver_seu@meta.data$condition[which(liver_seu@meta.data$condition == "Uninjured")] <- "Condition1"
liver_seu@meta.data$condition[which(liver_seu@meta.data$condition == "Cirrhotic")] <- "Condition2"


# Run Dawnn
liver_seu <- FindNeighbors(liver_seu, dims = 1:11, reduction = "PCA",
                           return.neighbor = TRUE, k.param = 1001)
# PDA
liver_seu <- run_dawnn(liver_seu, label_names = "condition",
                       recalculate_graph = FALSE, reduced_dim = "PCA",
                       label_1 = "Condition1", label_2 = "Condition2",
                       nn_model = "/root/.dawnn/dawnn_nn_model.h5",
                       tf_conda_env="tf_env", da_mode="pda", reduced_dim = "foo")
# Invert Dawnn's LFCs to match milo's output
liver_seu@meta.data$dawnn_pda_lfc_rev <- -1 * liver_seu@meta.data$dawnn_lfc
liver_seu@meta.data$dawnn_pda_verdict <- liver_seu@meta.data$dawnn_da_verdict

# ADA
liver_seu <- run_dawnn(liver_seu, label_names = "condition",
                       recalculate_graph = FALSE, reduced_dim = "PCA",
                       label_1 = "Condition1", label_2 = "Condition2",
                       nn_model = "/root/.dawnn/dawnn_nn_model.h5",
                       tf_conda_env="tf_env", da_mode="ada")
# Invert Dawnn's LFCs to match milo's output
liver_seu@meta.data$dawnn_ada_lfc_rev <- -1 * liver_seu@meta.data$dawnn_lfc
liver_seu@meta.data$dawnn_ada_verdict <- liver_seu@meta.data$dawnn_da_verdict


# Run DAseq
daseq_out <- getDAcells(liver_seu@reductions[["PCA"]]@cell.embeddings,
                        cell.labels = liver_seu@meta.data$condition,
                        labels.1 = c("Condition1"), labels.2 = c("Condition2"),
                        k.vector = seq(50, 500, 50))
daseq_norm <- (daseq_out$da.pred / 2) + 0.5
daseq_out$lfc <- log2(daseq_norm / (1 - daseq_norm)) - log2(0.5762363/(1-0.5762363)) # is it correct to subtract here?
daseq_out$is_da <- seq(1, length(daseq_out$cell.idx)) %in% c(daseq_out$da.up, daseq_out$da.down)


# Run Milo
liver_milo <- Milo(liver_sce)
liver_milo <- buildGraph(liver_milo, d = 11, k=30)
liver_milo <- makeNhoods(liver_milo, k=30, d=11, prop = 0.05, refined=TRUE)

colData(liver_milo)[['sort']] <- str_remove(colData(liver_milo)[['dataset']], ".+_")
colData(liver_milo)[['sort']] <- str_remove(colData(liver_milo)[['sort']], "A|B")
liver_meta <- as.tibble(colData(liver_milo)[,c("dataset","condition", 'sort')])
liver_meta <- distinct(liver_meta) %>%
  mutate(condition=factor(condition, levels=c("Uninjured", "Cirrhotic"))) %>%
  column_to_rownames("dataset")

liver_milo <- miloR::countCells(liver_milo, samples = "dataset",
                         meta.data = data.frame(colData(liver_milo)[,c("dataset","condition",'sort')]) )
liver_milo <- calcNhoodDistance(liver_milo, d = 11)
milo_res <- testNhoods(liver_milo, design = ~ condition,
                       design.df = liver_meta[colnames(nhoodCounts(liver_milo)), ])

milo_res <- milo_res[,!str_detect(colnames(milo_res), "annotation_lineage")]
milo_res <- annotateNhoods(liver_milo, milo_res, "annotation_indepth")
anno_df <- data.frame(liver_milo@colData) %>%
  distinct(annotation_lineage, annotation_indepth)
milo_res <- left_join(milo_res, anno_df, by="annotation_indepth")
milo_res$annotation_indepth[milo_res$annotation_indepth_fraction < 0.6] <- NA
milo_res$annotation_lineage[milo_res$annotation_indepth_fraction < 0.6] <- NA

milo_out_df <- data.frame(lfc = milo_res$logFC,
                          is_da = milo_res$SpatialFDR < 0.1,
                          annotation_indepth = milo_res$annotation_indepth,
                          annotation_lineage = milo_res$annotation_lineage)

daseq_out_df <- data.frame(lfc = daseq_out$lfc,
                           is_da = daseq_out$is_da,
                           annotation_indepth = liver_seu$annotation_indepth,
                           annotation_lineage = liver_seu$annotation_lineage)

dawnn_pda_out_df <- data.frame(lfc = liver_seu$dawnn_pda_lfc_rev,
                               is_da = liver_seu$dawnn_pda_verdict,
                               annotation_indepth = liver_seu$annotation_indepth,
                               annotation_lineage = liver_seu$annotation_lineage)
dawnn_ada_out_df <- data.frame(lfc = liver_seu$dawnn_ada_lfc_rev,
                               is_da = liver_seu$dawnn_ada_verdict,
                               annotation_indepth = liver_seu$annotation_indepth,
                               annotation_lineage = liver_seu$annotation_lineage)

all_methods_df <- rbind(cbind(milo_out_df, method = "milo"),
                        cbind(daseq_out_df, method = "daseq"),
                        cbind(dawnn_pda_out_df, method = "dawnn_pda"),
                        cbind(dawnn_ada_out_df, method = "dawnn_ada"))
write.csv(all_methods_df, "liver_cirrhosis_results_rerun.csv")
