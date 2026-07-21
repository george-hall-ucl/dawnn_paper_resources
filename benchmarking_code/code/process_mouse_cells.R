# Following https://github.com/MarioniLab/milo_analysis_2020/blob/dc8a30f4025396a17e40ee1708bb023e60a89072/benchmark/Fig2_benchmark_results.Rmd

library(MouseGastrulationData)
library(SingleCellExperiment)
library(scater)
library(scran)

# Starting from https://github.com/MarioniLab/milo_analysis_2020/blob/dc8a30f4025396a17e40ee1708bb023e60a89072/benchmark/Fig2_benchmark_results.Rmd#L112

late_samples <- AtlasSampleMetadata %>%
  filter(stage %in% c("E7.75", "E8.0", "E8.25", "E8.5")) %>%
  pull(sample)

embryo_sce <- EmbryoAtlasData(type="processed", samples = late_samples)
logcounts(embryo_sce) <- log1p(counts(embryo_sce))

keep.rows <- rowSums(logcounts(embryo_sce)) != 0
embryo_sce <- embryo_sce[keep.rows, ]
dec <- modelGeneVar(embryo_sce)
hvgs <- getTopHVGs(dec, n=5000)
embryo_sce <- embryo_sce[,apply(reducedDim(embryo_sce, "pca.corrected"), 1, function(x) all(!is.na(x)))]
embryo_sce <- scater::runUMAP(embryo_sce, dimred="pca.corrected")
