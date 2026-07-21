# Data from "Integrated multi-omic characterization of congenital heart disease"
# Nature 608 pp. 181-191 (2022).

library(dplyr)
library(Seurat)

heart_expr_mat <- ReadMtx(mtx = "heart_expression_matrix.mtx.gz",
                          features = "heart_genes.tsv.gz",
                          cells = "heart_barcodes.tsv.gz",
                          feature.column = 1)
cells <- CreateSeuratObject(counts = heart_expr_mat)
cells <- cells %>%
         FindVariableFeatures() %>%
         NormalizeData() %>%
         ScaleData() %>%
         RunPCA() %>%
         RunUMAP(1:20) %>%
         FindNeighbors() %>%
         FindClusters(resolution = 1.5)

saveRDS(cells, "heart_tissue_cells.RDS")

# upregulated cells: c("1", "9", "14", "20", "22", "24", "29", "31")
# then run: generate_benchmarking_data(name = "heart", rds_file = "regen_heart_dataset/heart_tissue_cells.rds", reduced_dim = "pca", cell_labels = "seurat_clusters")
