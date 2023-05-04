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
                      two_tailed = FALSE, k_vect = NULL) {
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

