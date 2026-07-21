# Data from "Cholangiocyte organoids can repair bile ducts after transplantation in the human liver"
# SCIENCE 19 Feb 2021 Vol 371, Issue 6531 pp. 839-846

# Download data from https://www.ebi.ac.uk/gxa/sc/experiment/E-MTAB-8495/download/zip?fileType=normalised&accessKey=
# Unzip downloaded file and save in directory "organoid_data"

d <- Matrix::readMM('organoid_data/E-MTAB-8495.aggregated_filtered_normalised_counts.mtx')

col_names <- read.csv("organoid_data/E-MTAB-8495.aggregated_filtered_normalised_counts.mtx_cols", header = F)
row_names <- read.csv("organoid_data/E-MTAB-8495.aggregated_filtered_normalised_counts.mtx_rows", header = F, sep = "\t")
rownames(d) <- row_names[,1]
colnames(d) <- col_names[,1]

library(Seurat)
library(dplyr)

set.seed(123)
cells <- CreateSeuratObject(counts = d) %>%
         FindVariableFeatures() %>%
         NormalizeData() %>%
         ScaleData() %>%
         RunPCA() %>%
         RunUMAP(1:20) %>%
         FindNeighbors(cells) %>%
         FindClusters()

saveRDS(cells, "organoid_tissue_cells_regen.RDS")


# upregulated cells: c("0", "3", "4", "6", "9", "10", "12", "17")
# then run: generate_benchmarking_data(name = "organoid", rds_file = "organoid_tissue_cells.RDS", reduced_dim = "PCA", cell_labels = "seurat_clusters")
