#!/usr/bin/Rscript

library(foreach)
library(doParallel)

library(Seurat)
library(miloR)
library(scater)
library(DAseq)
library(dawnn)

source("milo_paper_benchmark_utils.R")
source("benchmark_run_methods.R")

evaluate_results <- function(verdicts, cells = NULL, ground_truth = NULL,
                             da_threshold = 0.55, two_tailed = FALSE) {
    if (is.null(ground_truth)) {
        if (two_tailed == FALSE) {
            ground_truth <- cells@meta.data$Condition1_prob > da_threshold
        } else {
            ground_truth <- (cells@meta.data$Condition1_prob > da_threshold) |
                            (cells@meta.data$Condition1_prob < 1 - da_threshold)
        }
    }
    tpr <- sum((ground_truth == TRUE) & (verdicts == TRUE)) / sum((ground_truth == TRUE))
    if (tpr > 0) {
        fpr <- sum((ground_truth == FALSE) & (verdicts == TRUE)) / sum((ground_truth == FALSE))
        fdr <- sum((ground_truth == FALSE) & (verdicts == TRUE)) / sum((verdicts == TRUE))
    } else {
        fpr <- 0
        fdr <- 0
    }

    return(c(tpr, fpr, fdr))
}


add_to_results_df <- function(results_df, cells, method_name, verdicts,
                              time_taken, upreg_pc1, upreg_pop, rep, family,
                              da_threshold = 0.55, verbose = TRUE) {

    eval <- evaluate_results(cells, verdicts, da_threshold)
    this_result <- data.frame(method = method_name, tpr = eval[1],
                              fpr = eval[2], fdr = eval[3],
                              time = time_taken, upreg_pc1 = upreg_pc1,
                              upreg_pop = upreg_pop, rep = rep, family = family,
                              da_threshold = da_threshold)
    if (verbose) {print(this_result)}
    results_df <- rbind(results_df, this_result)

    return(results_df)
}


benchmark_dawnn <- function(cells, thing_to_bench, reduced_dim) {
        message("Benchmarking Dawnn")

        if (thing_to_bench == "accuracy") {
            start_time <- Sys.time()
            cells <- run_dawnn(cells, label_names = "synth_labels",
                               label_1 = "Condition1", label_2 = "Condition2",
                               reduced_dim = reduced_dim,
                               recalculate_graph = FALSE,
                               nn_model = "final_model_dawnn_rerun.h5")
            end_time <- Sys.time()
            time_taken <- as.numeric(end_time - start_time)
            tpr_fpr_fdr <- evaluate_results(verdicts = cells@meta.data$dawnn_da_verdict,
                                            ground_truth = cells@meta.data$is_da)
            result <- c("dawnn", "dawnn", tpr_fpr_fdr, time_taken)
        } else if (thing_to_bench == "MSE") {
            cells <- run_dawnn(cells, label_names = "synth_labels",
                               label_1 = "Condition1", label_2 = "Condition2",
                               reduced_dim = "pca.corrected",
                               recalculate_graph = FALSE,
                               nn_model = "final_model_dawnn_rerun.h5")
            result <- cells$dawnn_scores
        }


        return(result)
}


benchmark_milo <- function(cells, milo_graph, d, k, prop, reduced_dim) {
        message("Benchmarking Milo")

        start_time <- Sys.time()
        cells <- run_milo(cells, milo_graph, d, k, prop, reduced_dim)
        end_time <- Sys.time()
        time_taken <- as.numeric(end_time - start_time)

        tpr_fpr_fdr <- evaluate_results(verdicts = cells@meta.data$milo_verdict,
                                        ground_truth = cells@meta.data$is_da)

        result <- c("milo", "milo", tpr_fpr_fdr, time_taken)

        return(result)
}


benchmark_daseq <- function(cells, thing_to_bench, reduced_dim, two_tailed,
                            k_vect) {
        message("Benchmarking DAseq")

        start_time <- Sys.time()
        cells <- run_daseq(cells, reduced_dim = reduced_dim,
                           two_tailed = two_tailed, k_vect = k_vect)
        if (thing_to_bench == "accuracy") {
            end_time <- Sys.time()
            time_taken <- as.numeric(end_time - start_time)

            tpr_fpr_fdr <- evaluate_results(verdicts = cells@meta.data$daseq_verdict,
                                            ground_truth = cells@meta.data$is_da)

            result <- c("daseq", "daseq", tpr_fpr_fdr, time_taken)
        } else if (thing_to_bench == "MSE") {
            result <- cells@meta.data$daseq_pc1
        }

        return(result)
}


benchmark_meld <- function(cells, thing_to_bench) {
        message("Benchmarking MELD")

        neighbor_labels <- generate_neighbor_labels(cells, "pca.corrected", label_names = "synth_labels", recalculate_graph = FALSE, verbose = TRUE)
        scores <- nn_model %>% predict(neighbor_labels, verbose = 2)
        mse <- mean((scores[,1] - cells@meta.data$Condition1_prob)^2)
        result <- scores[,1]


        return(result)
}


benchmark_accuracy <- function(cells, instance, method, model_path,
                               upreg_pc1, upreg_pop, rep, milo_graph,
                               out_file_name, methods_params) {

    reduced_dim = methods_params[["all_reduced_dim"]]

    if (method == "dawnn") {
        result <- c(benchmark_dawnn(cells, "accuracy", reduced_dim),
                    upreg_pc1, upreg_pop, rep)
    }

    if (method == "milo") {
        replicates <- paste0("R", 1:3)
        cells@meta.data$synth_samples <- factor(paste0(cells@meta.data$synth_labels,
                                                "_", replicates))

        colData(milo_graph)$synth_labels <- cells@meta.data$synth_labels
        colData(milo_graph)$synth_samples <- cells@meta.data$synth_samples
        milo_out <- benchmark_milo(cells, milo_graph,
                                   methods_params[["milo_d"]],
                                   methods_params[["milo_k"]],
                                   methods_params[["milo_prop"]],
                                   reduced_dim = toupper(reduced_dim))
        result <- c(milo_out, upreg_pc1, upreg_pop, rep)
    }

    if (method == "daseq") {
        result <- c(benchmark_daseq(cells, "accuracy", reduced_dim = reduced_dim,
                                    two_tailed = FALSE,
                                    k_vect = methods_params[["daseq_k_vect"]]),
                    upreg_pc1, upreg_pop, rep)
    }

    return(result)
}


benchmark_mse <- function(cells, instance, method,
                          out_file_name, upreg_pc1, upreg_pop, rep, methods_params,
                          nn_model = NULL) {

    num_cells <- dim(cells)[2]
    labels <- as.numeric(instance[5:(4 + num_cells)])
    labels[which(labels == 1)] <- "Condition1"
    labels[which(labels == 0)] <- "Condition2"
    cells@meta.data$synth_labels <- labels

    if (method == "dawnn") {
        estimated_pc1s <- benchmark_dawnn(cells, "MSE", methods_params[["all_reduced_dim"]])
    } else if (method == "daseq") {
        estimated_pc1s <- benchmark_daseq(cells, "MSE",
                                          reduced_dim = methods_params[["all_reduced_dim"]],
                                          two_tailed = FALSE,
                                          k_vect = methods_params[["daseq_k_vect"]])
    }
    result <- c(estimated_pc1s, upreg_pc1, upreg_pop)

    return(result)
}


benchmark_time <- function(cells, instance, method, nn_model,
                           reduced_dim, out_file_name, methods_params,
                           downsample_size, rep, num_cells) {

    set.seed(rep)
    cells <- cells[, sample(1:num_cells, downsample_size)]
    if (method == "dawnn") {
        f <- function(cells, reduced_dim) {
            cells <- FindNeighbors(cells, dims = 1:50, return.neighbor = TRUE,
                                   k.param = 1001, reduction = reduced_dim)
            cells <- run_dawnn(cells, label_names = "synth_labels",
                               label_1 = "Condition1", label_2 = "Condition2",
                               reduced_dim = reduced_dim,
                               recalculate_graph = FALSE,
                               nn_model = "final_model_dawnn_rerun.h5")
        }
        timed <- system.time(f(cells, reduced_dim))
        time_taken <- sum(timed[1:2])
    } else if (method == "milo") {
        cells_sce <- as.SingleCellExperiment(cells)
        f <- function(cells_sce) {
            milo_graph <- build_milo_graph_from_sce(cells_sce,
                                                    methods_params[["milo_d"]],
                                                    methods_params[["milo_k"]],
                                                    methods_params[["milo_prop"]],
                                                    reduced_dim = toupper(reduced_dim))
            replicates <- paste0("R", 1:3)
            cells@meta.data$synth_samples <- factor(paste0(cells@meta.data$synth_labels,
                                                    "_", replicates))

            colData(milo_graph)$synth_labels <- cells@meta.data$synth_labels
            colData(milo_graph)$synth_samples <- cells@meta.data$synth_samples

            run_milo(cells, milo_graph, methods_params[["milo_d"]],
                     methods_params[["milo_k"]],
                     methods_params[["milo_prop"]], toupper(reduced_dim))

        }
        timed <- system.time(f(cells_sce))
        time_taken <- sum(timed[1:2])
    } else if (method == "daseq") {
        f <- function(cells, reduced_dim) {
            run_daseq(cells, reduced_dim = reduced_dim,
                      k_vect = methods_params[["daseq_k_vect"]])
        }
        timed <- system.time(f(cells, reduced_dim))
        time_taken <- sum(timed[1:2])
    }

    return(c(method, downsample_size, time_taken, rep))
}


load_datasets <- function(dataset, cells_objs_paths, num_cells_objs) {
    if (dataset %in% c("mouse_embryo", "mouse_embryo_mse")) {
        cells <- readRDS(cells_objs_paths[1])
        cells <- as.Seurat(cells)
        cells_objs <- c(cells)
    } else if (dataset == "skin") {
        cells <- readRDS(cells_objs_paths[1])
        cells_objs <- c(cells)
    } else if (dataset == "organoid") {
        cells <- readRDS(cells_objs_paths[1])
        cells_objs <- c(cells)
    } else if (dataset == "heart") {
        cells <- readRDS(cells_objs_paths[1])
        cells_objs <- c(cells)
    } else if (dataset == "sim_discrete_clusters") {
        cells_objs <- c()
        for (i in 1:num_cells_objs) {
            cells <- readRDS(paste0("milo_sim_dat/cells_sim_discerete_clusters_gex_seed_", i, ".rds"))
            cells_objs <- c(cells_objs, cells)
        }
    } else if (dataset == "sim_linear_traj") {
        cells_objs <- c()
        for (i in 1:num_cells_objs) {
            cells <- readRDS(paste0("milo_sim_dat/cells_sim_linear_traj_gex_seed_", i, ".rds"))
            cells_objs <- c(cells_objs, cells)
        }
    } else if (dataset == "sim_branch_traj") {
        cells_objs <- c()
        for (i in 1:num_cells_objs) {
            cells <- readRDS(paste0("milo_sim_dat/cells_sim_branching_traj_gex_seed_", i, ".rds"))
            cells_objs <- c(cells_objs, cells)
        }
    } else {
        stop(paste("Unrecognized dataset:", dataset))
    }

    return(cells_objs)
}


set_methods_params <- function(dataset) {
    methods_params <- list(all_reduced_dim = NULL,
                           dawnn_alpha = 0.1,
                           milo_d = 30,
                           milo_k = NULL,
                           milo_prop = 0.1,
                           daseq_k_vect = NULL,
                           meld_k = 50,
                           meld_d = 30)

    if (dataset %in% c("mouse_embryo", "mouse_embryo_mse")) {
        methods_params[["all_reduced_dim"]] <- "pca.corrected"
        methods_params[["milo_k"]] <- 50
        methods_params[["daseq_k_vect"]] <- seq(50, 500, 50)
    } else if (dataset == "skin") {
        methods_params[["all_reduced_dim"]] <- "harmony"
        methods_params[["milo_k"]] <- 50
        methods_params[["daseq_k_vect"]] <- seq(50, 500, 50)
    } else if (dataset == "organoid") {
        methods_params[["all_reduced_dim"]] <- "pca"
        methods_params[["milo_k"]] <- 50
        methods_params[["daseq_k_vect"]] <- seq(50, 500, 50)
    } else if (dataset == "heart") {
        methods_params[["all_reduced_dim"]] <- "pca"
        methods_params[["milo_k"]] <- 50
        methods_params[["daseq_k_vect"]] <- seq(50, 500, 50)
    } else if (dataset == "sim_discrete_clusters") {
        methods_params[["all_reduced_dim"]] <- "pca"
        methods_params[["milo_k"]] <- 15
        methods_params[["daseq_k_vect"]] <- seq(15, 500, 50)
    } else if (dataset %in% c("sim_linear_traj", "sim_branch_traj")) {
        methods_params[["all_reduced_dim"]] <- "pca"
        methods_params[["milo_k"]] <- 20
        methods_params[["daseq_k_vect"]] <- seq(20, 500, 50)
    } else {
        stop(paste("Unrecognized dataset:", dataset))
    }

    return(methods_params)
}


write_to_csv <- function(result, out_file_name) {
    write.table(t(result), out_file_name, append = TRUE, sep = ",",
                col.names = FALSE, row.names = FALSE)
}


run_benchmarking <- function(methods_to_benchmark, thing_to_bench, dataset,
                             out_file_name, num_cpus, cell_type = NULL) {

    if (dataset == "mouse_embryo") {
        # Spaces in cell type names have been replaced with underscores
        test_data_path <- "regen_mouse_dataset/benchmark_dataset_mouse.csv"
        num_extra_fields <- 4
    } else if (dataset == "mouse_embryo_mse") {
        test_data_path <- "benchmark_data/benchmark_dataset_for_pc1_estimation.csv"
        num_extra_fields <- 4
    } else if (dataset == "skin") {
        test_data_path <- "regen_skin_dataset/benchmark_dataset_skin.csv"
        num_extra_fields <- 4
    } else if (dataset == "organoid") {
        test_data_path <- "regen_organoid_dataset/benchmark_dataset_organoid_labels.csv"
        num_extra_fields <- 4
    } else if (dataset == "heart") {
        test_data_path <- "regen_heart_dataset/benchmark_dataset_heart_data_type_labels.csv"
        num_extra_fields <- 4
    } else if (dataset == "sim_discrete_clusters") {
        # test_data_path <- "milo_sim_dat/benchmark_dataset_sim_discrete_clusters.csv"
        test_data_path <- "milo_sim_dat/benchmark_dataset_sim_discrete_clusters.csv"
        num_extra_fields <- 4
    } else if (dataset == "sim_linear_traj") {
        test_data_path <- "milo_sim_dat/benchmark_dataset_sim_linear_traj.csv"
        num_extra_fields <- 5
    } else if (dataset == "sim_branch_traj") {
        test_data_path <- "milo_sim_dat/benchmark_dataset_sim_branching_traj.csv"
        num_extra_fields <- 5
    }

    test_data <- read.csv(test_data_path, header = FALSE)
    if (!is.null(cell_type)) {
        test_data <- subset(test_data, V2 %in% cell_type)
    }
    cells_objs_paths <- unique(test_data[,1])
    num_cells_objs <- length(cells_objs_paths)
    num_cells <- (dim(test_data)[2] - num_extra_fields) / 2
    num_instances_per_cells_obj <- dim(test_data)[1] / length(cells_objs_paths)

    # Load benchmarking datasets
    cells_objs <- load_datasets(dataset, cells_objs_paths, num_cells_objs)

    methods_params <- set_methods_params(dataset)
    reduced_dim <- methods_params[["all_reduced_dim"]]

    results <- list()

    # Loop over loaded benchmarking datasets
    for (j in seq_along(cells_objs)) {
        cells <- cells_objs[[j]]

        if ("dawnn" %in% methods_to_benchmark & thing_to_bench != "time") {
            cells <- FindNeighbors(cells, dims = 1:50, return.neighbor = TRUE,
                                   k.param = 1001, reduction = reduced_dim)
        }

        if ("milo" %in% methods_to_benchmark & thing_to_bench != "time") {
            message("Constructing milo graph.")
            milo_graph <- build_milo_graph(cells, methods_params[["milo_d"]],
                                           methods_params[["milo_k"]],
                                           methods_params[["milo_prop"]],
                                           reduced_dim = toupper(reduced_dim))
        } else {
            milo_graph <- NULL
        }

        for (method in methods_to_benchmark) {

            registerDoParallel(num_cpus)

            # Loop over instances of loaded dataset in parallel
            instance_start_idx <- ((j - 1) * (num_instances_per_cells_obj)) + 1
            instance_end_idx <- (num_instances_per_cells_obj * j)
            foreach (l = instance_start_idx:instance_end_idx) %dopar% {
                out_file_name <- paste0("thread_", (l %% num_cpus), "_", out_file_name)
                instance <- test_data[l, ]
                upreg_pop <- as.character(instance[2])
                upreg_pc1 <- as.numeric(instance[3])
                rep <- as.numeric(instance[4])

                labels <- as.numeric(instance[(num_extra_fields + 1):(num_extra_fields + num_cells)])
                labels[which(labels == 1)] <- "Condition1"
                labels[which(labels == 0)] <- "Condition2"
                cells@meta.data$synth_labels <- labels

                ground_truth <- (as.numeric(instance[(num_extra_fields + num_cells + 1):length(instance)]) == 1)
                cells@meta.data$is_da <- ground_truth

                if (thing_to_bench == "accuracy") {
                    res <- benchmark_accuracy(cells, instance, method,
                                              nn_model, upreg_pc1, upreg_pop, rep,
                                              milo_graph, out_file_name, methods_params)
                    write_to_csv(res, out_file_name)
                    results <- append(results, list(res))
                } else if (thing_to_bench == "MSE") {
                    res <- benchmark_mse(cells, instance, method,
                                         out_file_name, upreg_pc1, upreg_pop, rep,
                                         methods_params, nn_model = nn_model)
                    write_to_csv(res, out_file_name)
                    results <- append(results, list(res))
                } else if (thing_to_bench == "time") {
                    for (downsample_size in seq(5000, 60000, 5000)) {
                        for (rep in 1:3) {
                            res <- benchmark_time(cells, instance, method,
                                                  nn_model, reduced_dim,
                                                  out_file_name,
                                                  methods_params,
                                                  downsample_size, rep,
                                                  num_cells)
                            write_to_csv(res, out_file_name)
                            results <- append(results, list(res))
                        }
                    }
                }
            }
            stopImplicitCluster()
        }
    }

    # Combine output files
    system2("cat", args = paste0("thread_*_", out_file_name),
            stdout = out_file_name)
}
