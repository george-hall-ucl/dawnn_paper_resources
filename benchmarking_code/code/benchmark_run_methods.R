build_milo_graph <- function(cells, d, k, prop, reduced_dim = "PCA.CORRECTED") {
    cells_sce <- Milo(as.SingleCellExperiment(cells))
    cells_sce <- buildGraph(cells_sce, reduced.dim = reduced_dim, k = k,
                            d = d)
    cells_sce <- makeNhoods(cells_sce, reduced_dims = reduced_dim, k = k,
                            d = d, prop = prop)

    return(cells_sce)
}


build_milo_graph_from_sce <- function(cells_sce, d, k, prop, reduced_dim = "PCA.CORRECTED") {
    cells_sce <- Milo(cells_sce)
    cells_sce <- buildGraph(cells_sce, reduced.dim = reduced_dim, k = k,
                            d = d)
    cells_sce <- makeNhoods(cells_sce, reduced_dims = reduced_dim, k = k,
                            d = d, prop = prop)

    return(cells_sce)
}


run_milo <- function(cells, cells_sce, d, k, prop,
                     reduced_dim = "PCA.CORRECTED") {
    cells_sce <- miloR::countCells(cells_sce, meta.data = data.frame(colData(cells_sce)), sample = "synth_samples")
    cells_sce <- miloR::calcNhoodDistance(cells_sce,
                                          reduced.dim = reduced_dim, d = d)
    metadata_df <- data.frame(colData(cells_sce))[, c("synth_samples", "synth_labels")]
    design_df <- unique(metadata_df)
    rownames(design_df) <- design_df$synth_samples
    milo_out <- testNhoods(cells_sce, design = ~ synth_labels,
                           design.df = design_df, reduced.dim = reduced_dim)
    cells@meta.data$milo_verdict <- (milo2output(cells_sce, milo_out, out_type = "label") != "NotDA")


    return(cells)
}


run_daseq <- function(cells, reduced_dim = "pca.corrected",
                      two_tailed = TRUE, k_vect = NULL) {
    # labels.1 corresponds to Condition2 to match how the other methods treat
    # up-/down-regulation.
    daseq_verdict <- getDAcells(cells@reductions[[reduced_dim]]@cell.embeddings,
                                cell.labels = cells@meta.data$synth_labels,
                                labels.1 = c("Condition2"),
                                labels.2 = c("Condition1"), k.vector = k_vect)

    cells@meta.data$daseq_pc1 <- as.numeric(locfit::expit(daseq_verdict$da.pred))

    cells@meta.data$daseq_verdict <- FALSE
    cells@meta.data$daseq_verdict[daseq_verdict$da.up] <- TRUE
    if (two_tailed) {
        cells@meta.data$daseq_verdict[daseq_verdict$da.down] <- TRUE
    }

    return(cells)
}


run_cna <- function(cells, n_neighbors, cpu_number,
                    reduced_dim = "pca.corrected",
                    thing_to_bench = "accuracy") {

    # Prevent a bug in SaveH5Seurat if there's an SCT assay
    if ("SCT" %in% names(cells@assays)) {
        slot(cells$SCT@SCTModel.list[[1]], 'median_umi') = median(cells$SCT@SCTModel.list[[1]]@cell.attributes$umi)
    }

    h5_filename_base <- paste0("cna_input_", cpu_number)
    # Save cells to anndata object

    # First, fix bug if trying to save a Seurat v5 object (as per https://github.com/mojaveazure/seurat-disk/issues/30#issuecomment-1887069144)
    cells[["RNA3"]] <- as(object = cells[["RNA"]], Class = "Assay")
    DefaultAssay(cells) <- "RNA3"
    cells[["RNA"]] <- NULL
    cells[["RNA"]] <- cells[["RNA3"]]
    DefaultAssay(cells) <- "RNA"
    cells[["RNA3"]] <- NULL

    SeuratDisk::SaveH5Seurat(cells,
                             filename = paste0(h5_filename_base, ".h5Seurat"),
                             overwrite = TRUE)
    SeuratDisk::Convert(paste0(h5_filename_base, ".h5Seurat"), dest = "h5ad", overwrite = TRUE)

    cmd <- paste0('eval "$(conda shell.bash hook)"; conda activate cna; ',
                  "python3 run_cna.py ", h5_filename_base, " ", reduced_dim,
                  " ", n_neighbors, " ", thing_to_bench, "; ",
                  "conda deactivate")
    system(cmd)
    cna_out <- read.csv(paste0(h5_filename_base, "_cna_out.csv"),
                        header = FALSE)
    if (thing_to_bench %in% c("accuracy", "accuracy_verbose")) {
        cna_out <- as.logical(cna_out[, 1])
        cells@meta.data$cna_verdict <- cna_out
        return(cells)
    } else if (thing_to_bench %in% c("MSE", "time")) {
        cna_out <- as.numeric(cna_out[, 1])
        return(cna_out)
    }
}
