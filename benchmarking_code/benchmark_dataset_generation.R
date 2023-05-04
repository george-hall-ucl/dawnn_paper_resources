source("milo_paper_benchmark_utils.R")

add_synth_labels <- function(cells, upreg_pop, pop_enr, reduced_dim = "PCA",
                             cell_labels = "seurat_clusters", method = "pop", seed = 42, m = 2) {
    if (method == "pop") {
        synthed <- as.Seurat(add_synthetic_labels_pop(as.SingleCellExperiment(cells),
                                                      pop = upreg_pop,
                                                      pop_enr = pop_enr,
                                                      redDim = reduced_dim,
                                                      pop_column = cell_labels, seed = seed, m = m))
        cells@meta.data <- synthed@meta.data
        cells@meta.data$is_da <- cells@meta.data$Condition1_prob > 0.55
    }
    else if (method == "cluster") {
        synthed <- as.Seurat(add_synthetic_labels_by_cluster(as.SingleCellExperiment(cells),
                                                      pop = upreg_pop,
                                                      pop_enr = pop_enr,
                                                      redDim = reduced_dim,
                                                      pop_column = cell_labels, seed = seed, m = m))
        cells@meta.data <- synthed@meta.data
        cells@meta.data$is_da <- cells@meta.data$Condition1_prob > 0.5
    }
    
    return(cells)
}


generate_benchmarking_data_traj <- function(upreg_pop_idx, data_type = "labels") {
    pops <- c("Caudal neurectoderm", "Erythroid2", "Gut", "Somitic mesoderm",
              "Pharyngeal mesoderm", "Erythroid1", "Mesenchyme", "ExE endoderm")
    upreg_pop <- pops[upreg_pop_idx]
   
    if (data_type == "labels") {
        upreg_pc1s <- seq(0.7, 0.95, 0.05)
        num_reps <- 3
    } else if (data_type == "pc1s") {
        upreg_pc1s <- c(0.7, 0.8, 0.9)
        num_reps <- 1
    } else {
        stop(paste0("Invalid data_type:", data_type))
    }

    for (pc1 in upreg_pc1s) {
        for (rep in 1:num_reps) {
            cells <- add_synth_labels(cells, upreg_pop, pc1, seed = rep,
                              reduced_dim = "PCA.CORRECTED",
                              cell_labels = "celltype")
            if (data_type == "labels") {
                out_df <- data.frame(t(c(c(upreg_pop, pc1, rep),
                                         as.numeric(cells@meta.data$synth_labels == "Condition1"),
                                         as.numeric(cells@meta.data$is_da))))
            } else if (data_type == "pc1s") {
                out_df <- data.frame(t(c(c(upreg_pop, pc1, rep),
                                         cells@meta.data$Condition1_prob,
                                         cells@meta.data$synth_labels)))
            }
            file_name <- paste0("benchmark_dataset_", upreg_pop, "_data_type_",
                                data_type, ".csv")
            write.table(out_df, file_name, append = TRUE, sep = ",",
                        col.names = FALSE, row.names = FALSE)
        }
    }
}


generate_benchmarking_data <- function(name, rds_file, reduced_dim, cell_labels, data_type = "labels") {
    cells <- readRDS(rds_file)

    if (name == "skin") {
        cells <- FindNeighbors(cells, dims = 1:50, reduction = "harmony") %>% FindClusters()
        upreg_pops <- c("1", "3", "4", "5", "7", "8", "10", "11")
    } else if (name == "organoid") {
        upreg_pops <- c("0", "3", "4", "6", "9", "10", "12", "17")
    } else if (name == "heart") {
        upreg_pops <- c("1", "9", "14", "20", "22", "24", "29", "31")
    } else if (name == "mouse_gastrulation") {
        upreg_pops <- c("Caudal neurectoderm", "Erythroid2", "Gut",
                        "Somitic mesoderm", "Pharyngeal mesoderm", "Erythroid1",
                        "Mesenchyme", "ExE endoderm")
        cells <- as.Seurat(cells)
    }

    if (data_type == "labels") {
        upreg_pc1s <- seq(0.7, 0.95, 0.05)
        num_reps <- 3
    } else if (data_type == "pc1s") {
        upreg_pc1s <- c(0.7, 0.8, 0.9)
        num_reps <- 1
    } else {
        stop(paste0("Invalid data_type:", data_type))
    }

    for (upreg_pop in upreg_pops) {
        for (pc1 in upreg_pc1s) {
            for (rep in 1:num_reps) {
                cells <- add_synth_labels(cells, upreg_pop, pc1, seed = rep,
                                  reduced_dim = reduced_dim,
                                  cell_labels = cell_labels)
                if (data_type == "labels") {
                    out_df <- data.frame(t(c(c(rds_file, upreg_pop, pc1, rep),
                                             as.numeric(cells@meta.data$synth_labels == "Condition1"),
                                             as.numeric(cells@meta.data$is_da))))
                } else if (data_type == "pc1s") {
                    out_df <- data.frame(t(c(c(upreg_pop, pc1, rep),
                                             cells@meta.data$Condition1_prob,
                                             cells@meta.data$synth_labels)))
                }
                file_name <- paste0("benchmark_dataset_", name, "_", upreg_pop, "_data_type_",
                                    data_type, ".csv")
                write.table(out_df, file_name, append = TRUE, sep = ",",
                            col.names = FALSE, row.names = FALSE)
            }
        }
    }
}

