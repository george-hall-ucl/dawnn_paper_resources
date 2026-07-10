library(foreach)
library(doParallel)

library(Seurat)
library(miloR)

library(DAseq)
library(dawnn)

source("milo_paper_benchmark_utils.R")
source("benchmark_run_methods.R")

evaluate_results <- function(cells, verdicts) {
    results_list <- list()
    for (da_mode in c("pda", "ada")) {
        ground_truth <- cells@meta.data[, paste0("is_", da_mode)]
        tpr <- sum((ground_truth == TRUE) & (verdicts == TRUE)) / sum((ground_truth == TRUE))
        if ((tpr > 0) & (any(ground_truth == FALSE))) {
            fpr <- sum((ground_truth == FALSE) & (verdicts == TRUE)) / sum((ground_truth == FALSE))
            fdr <- sum((ground_truth == FALSE) & (verdicts == TRUE)) / sum((verdicts == TRUE))
        } else {
            fpr <- 0
            fdr <- 0
        }
        results_list[[da_mode]] <- sprintf("%.2f", c(tpr, fpr, fdr))
    }

    return(unlist(results_list[c("pda", "ada")]))
}


benchmark_dawnn <- function(cells, thing_to_bench, reduced_dim,
                            method = NULL) {
        message("Benchmarking Dawnn")

        if (thing_to_bench %in% c("accuracy", "accuracy_verbose")) {
            da_mode <- ifelse(method == "dawnn_pda", "pda", "ada")
            cells <- run_dawnn(cells, label_names = "synth_labels",
                               label_1 = "Condition1", label_2 = "Condition2",
                               reduced_dim = reduced_dim,
                               recalculate_graph = FALSE,
                               da_mode = da_mode,
                               tf_conda_env = "tf_env")

           tpr_fpr_fdr <- evaluate_results(cells, cells$dawnn_da_verdict)
           if (thing_to_bench == "accuracy") {
               result <- c(method, method, tpr_fpr_fdr)
           } else if (thing_to_bench == "accuracy_verbose") {
               result <- cells@meta.data$dawnn_da_verdict
           }
        } else if (thing_to_bench == "MSE") {
            cells <- run_dawnn(cells, label_names = "synth_labels",
                               label_1 = "Condition1", label_2 = "Condition2",
                               reduced_dim = "pca.corrected",
                               recalculate_graph = FALSE,
                               tf_conda_env = "tf_env")
            result <- cells$dawnn_scores
        }


        return(result)
}


benchmark_milo <- function(cells, thing_to_bench, milo_graph, d, k, prop, reduced_dim) {
        message("Benchmarking Milo")

        cells <- run_milo(cells, milo_graph, d, k, prop, reduced_dim)

        tpr_fpr_fdr <- evaluate_results(cells, cells$milo_verdict)

        if (thing_to_bench == "accuracy") {
            result <- c("milo", "milo", tpr_fpr_fdr)
        } else if (thing_to_bench == "accuracy_verbose") {
            result <- cells@meta.data$milo_verdict
        }

        return(result)
}


benchmark_daseq <- function(cells, thing_to_bench, reduced_dim, two_tailed,
                            k_vect) {
        message("Benchmarking DAseq")

        cells <- run_daseq(cells, reduced_dim = reduced_dim,
                           two_tailed = two_tailed, k_vect = k_vect)
        if (thing_to_bench %in% c("accuracy", "accuracy_verbose")) {

            tpr_fpr_fdr <- evaluate_results(cells, cells$daseq_verdict)

            if (thing_to_bench == "accuracy") {
                result <- c("daseq", "daseq", tpr_fpr_fdr)
            } else if (thing_to_bench == "accuracy_verbose") {
                result <- cells@meta.data$daseq_verdict
            }

        } else if (thing_to_bench == "MSE") {
            result <- cells@meta.data$daseq_pc1
        }

        return(result)
}


benchmark_cna <- function(cells, n_neighbors, cpu_number, thing_to_bench, reduced_dim) {
        message("Benchmarking CNA")

        if (thing_to_bench %in% c("accuracy", "accuracy_verbose")) {
            cna_out <- run_cna(cells, n_neighbors, cpu_number,
                               reduced_dim = reduced_dim,
                               thing_to_bench = thing_to_bench)

            if (thing_to_bench == "accuracy") {
                tpr_fpr_fdr <- evaluate_results(cells, cna_out$cna_verdict)

                result <- c("cna", "cna", tpr_fpr_fdr)
            } else {
                result <- cna_out$cna_verdict
            }
        } else if (thing_to_bench == "MSE") {
            result <- run_cna(cells, n_neighbors, cpu_number,
                         reduced_dim = reduced_dim,
                         thing_to_bench = "MSE")
        } else if (thing_to_bench == "time") {
            result <- run_cna(cells, n_neighbors, cpu_number,
                         reduced_dim = reduced_dim,
                         thing_to_bench = "time")
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
                               out_file_name, methods_params, cpu_number,
                               thing_to_bench, downsampling) {

    reduced_dim <- methods_params[["all_reduced_dim"]]

    if (downsampling) {
        freq_tab <- table(cells$synth_labels)
        rarer_cluster <- names(freq_tab)[which.min(freq_tab)]
        common_cluster <- names(freq_tab)[which.max(freq_tab)]
        if (rarer_cluster == common_cluster) {
            stop("Error in downsampling code: rare and common clusters are the same.")
        }
        common_idxs <- which(cells$synth_labels == common_cluster)
        rarer_idxs <- which(cells$synth_labels == rarer_cluster)
        downsampled_idxs <- sample(common_idxs, freq_tab[rarer_cluster])
        new_idxs <- c(rarer_idxs, downsampled_idxs)
        cells <- cells[, new_idxs]
        gc()
    }

    if (method %in% c("dawnn_pda", "dawnn_ada")) {
        if (downsampling) {
            cells <- FindNeighbors(cells, dims = 1:50, return.neighbor = TRUE,
                                   k.param = 1001, reduction = reduced_dim)
        }
        result <- c(benchmark_dawnn(cells, thing_to_bench, reduced_dim,
                                    method),
                    upreg_pc1, upreg_pop, rep)
    } else if (method == "milo") {
        replicates <- paste0("R", 1:3)
        cells@meta.data$synth_samples <- factor(paste0(cells@meta.data$synth_labels,
                                                "_", replicates))

        colData(milo_graph)$synth_labels <- cells@meta.data$synth_labels
        colData(milo_graph)$synth_samples <- cells@meta.data$synth_samples
        milo_out <- benchmark_milo(cells, "accuracy", milo_graph,
                                   methods_params[["milo_d"]],
                                   methods_params[["milo_k"]],
                                   methods_params[["milo_prop"]],
                                   reduced_dim = toupper(reduced_dim))
        result <- c(milo_out, upreg_pc1, upreg_pop, rep)

    } else if (method == "daseq") {
        result <- c(benchmark_daseq(cells, thing_to_bench, reduced_dim = reduced_dim,
                                    two_tailed = TRUE,
                                    k_vect = methods_params[["daseq_k_vect"]]),
                    upreg_pc1, upreg_pop, rep)
    } else if (method == "cna") {
        result <- c(benchmark_cna(cells, methods_params[["milo_k"]],
                                  cpu_number, thing_to_bench,
                                  reduced_dim = reduced_dim),
                    upreg_pc1, upreg_pop, rep)
    } else {
        stop(paste("Unknown method:", method))
    }

    return(result)
}


benchmark_mse <- function(cells, instance, method,
                          out_file_name, upreg_pc1, upreg_pop, rep, methods_params, cpu_number,
                          nn_model = NULL) {

    num_cells <- dim(cells)[2]
    labels <- as.numeric(instance[5:(4 + num_cells)])
    labels[which(labels == 1)] <- "Condition1"
    labels[which(labels == 0)] <- "Condition2"
    cells@meta.data$synth_labels <- labels

    if (method %in% c("dawnn_ada", "dawnn_pda")) {
        estimated_pc1s <- benchmark_dawnn(cells, "MSE", methods_params[["all_reduced_dim"]])
    } else if (method == "daseq") {
        estimated_pc1s <- benchmark_daseq(cells, "MSE",
                                          reduced_dim = methods_params[["all_reduced_dim"]],
                                          two_tailed = FALSE,
                                          k_vect = methods_params[["daseq_k_vect"]])
    } else if (method == "cna") {
        estimated_pc1s <- benchmark_cna(cells, methods_params[["milo_k"]], cpu_number, "MSE", reduced_dim = methods_params[["all_reduced_dim"]])
    } else if (method == "meld") {
        estimated_pc1s <- benchmark_meld(cells, "MSE")
    }

    result <- c(estimated_pc1s, upreg_pc1, upreg_pop)

    return(result)
}


benchmark_time <- function(cells, instance, method, nn_model,
                           reduced_dim, out_file_name, methods_params,
                           downsample_size, rep, num_cells) {

    set.seed(rep)
    cells <- cells[, sample(1:num_cells, downsample_size)]
    if (method %in% c("dawnn_pda", "dawnn_ada")) {
        f <- function(cells, reduced_dim) {
            cells <- FindNeighbors(cells, dims = 1:50, return.neighbor = TRUE,
                                   k.param = 1001, reduction = reduced_dim)
            da_mode <- substr(method, 7, 10)
            cells <- run_dawnn(cells, label_names = "synth_labels",
                               label_1 = "Condition1", label_2 = "Condition2",
                               reduced_dim = reduced_dim,
                               recalculate_graph = FALSE,
                               da_mode = da_mode,
                               tf_conda_env = "tf_env")
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
    } else if (method == "cna") {
        time_taken <- run_cna(cells, n_neighbors = methods_params[["milo_k"]],
                           cpu_number = 1, reduced_dim = reduced_dim,
                           thing_to_bench = "time")
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
            cells <- readRDS(paste0("../data/outputs/discrete_clusters/cells_sim_discrete_clusters_gex_seed_", i, ".rds"))
            cells_objs <- c(cells_objs, cells)
        }
    } else if (dataset == "sim_linear_traj") {
        cells_objs <- c()
        for (i in 1:num_cells_objs) {
            cells <- readRDS(paste0("../data/outputs/linear_trajs/cells_sim_linear_trajs_gex_seed_", i, ".rds"))
            cells_objs <- c(cells_objs, cells)
        }
    } else if (dataset == "sim_branch_traj") {
        cells_objs <- c()
        for (i in 1:num_cells_objs) {
            cells <- readRDS(paste0("../data/outputs/branching_trajs/cells_sim_branching_trajs_gex_seed_", i, ".rds"))
            cells_objs <- c(cells_objs, cells)
        }
    } else if (dataset == "fully_sim") {
        cells_objs <- c()
        for (i in 1:num_cells_objs) {
            cells <- readRDS(paste0("../data/outputs/fully_sim/fully_sim_500k.rds"))
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
    } else if (dataset %in% c("sim_linear_traj", "sim_branch_traj", "fully_sim")) {
        methods_params[["all_reduced_dim"]] <- "pca"
        methods_params[["milo_k"]] <- 20
        methods_params[["daseq_k_vect"]] <- seq(20, 500, 50)
    } else {
        stop(paste("Unrecognized dataset:", dataset))
    }

    return(methods_params)
}


write_to_file <- function(result, out_file_name) {
    result_pretty <- paste(t(result), collapse = "\t")
    write(result_pretty, out_file_name, append = TRUE)
}


run_benchmarking <- function(method, thing_to_bench, dataset, out_file_name,
                             num_cpus, downsampling) {

    if (dataset == "mouse_embryo") {
        num_extra_fields <- 4
        in_file_name <- "../data/outputs/mouse/benchmark_dataset_mouse.csv"
    } else if (dataset == "mouse_embryo_mse") {
        num_extra_fields <- 4
        in_file_name <- "../data/outputs/mouse_mse/benchmark_dataset_mouse_mse.csv"
    } else if (dataset == "skin") {
        num_extra_fields <- 4
        in_file_name <- "../data/outputs/skin/benchmark_dataset_skin.csv"
    } else if (dataset == "organoid") {
        num_extra_fields <- 4
        in_file_name <- "../data/outputs/organoid/benchmark_dataset_organoid.csv"
    } else if (dataset == "heart") {
        num_extra_fields <- 4
        in_file_name <- "../data/outputs/heart/benchmark_dataset_heart.csv"
    } else if (dataset == "sim_discrete_clusters") {
        num_extra_fields <- 5
        in_file_name <- "../data/outputs/discrete_clusters/benchmark_dataset_sim_discrete_clusters.csv"
    } else if (dataset == "sim_linear_traj") {
        num_extra_fields <- 5
        in_file_name <- "../data/outputs/linear_trajs/benchmark_dataset_sim_linear_trajs.csv"
    } else if (dataset == "sim_branch_traj") {
        num_extra_fields <- 5
        in_file_name <- "../data/outputs/branching_trajs/benchmark_dataset_sim_branching_trajs.csv"
    } else if (dataset == "fully_sim") {
        num_extra_fields <- 5
        in_file_name <- "../data/outputs/fully_sim/benchmark_dataset_fully_sim.csv"
    }

    if (dataset == "fully_sim") {
        test_data <- as.data.frame(data.table::fread(in_file_name))
    } else {
        test_data <- read.csv(in_file_name, header = FALSE)
    }
    cells_objs_paths <- unique(test_data[,1])
    num_cells_objs <- length(cells_objs_paths)
    num_cells <- (dim(test_data)[2] - num_extra_fields) / 3
    num_instances_per_cells_obj <- dim(test_data)[1] / length(cells_objs_paths)

    # Load benchmarking datasets
    cells_objs <- load_datasets(dataset, cells_objs_paths, num_cells_objs)

    methods_params <- set_methods_params(dataset)
    reduced_dim <- methods_params[["all_reduced_dim"]]

    # Fix memory issues with large objects when running in parallel
    options(future.globals.maxSize = 1024**4)

    # Loop over loaded benchmarking datasets
    for (j in seq_along(cells_objs)) {
        cells <- cells_objs[[j]]

        dawnn_method <- method %in% c("dawnn", "dawnn_pda", "dawnn_ada")
        if (dawnn_method & (thing_to_bench != "time")) {
            cells <- FindNeighbors(cells, dims = 1:50, return.neighbor = TRUE,
                                   k.param = 1001, reduction = reduced_dim)
        }

        if ((method == "milo") & (thing_to_bench != "time")) {
            message("Constructing milo graph.")
            milo_graph <- build_milo_graph(cells, methods_params[["milo_d"]],
                                           methods_params[["milo_k"]],
                                           methods_params[["milo_prop"]],
                                           reduced_dim = toupper(reduced_dim))
        } else {
            milo_graph <- NULL
        }

        registerDoParallel(num_cpus)

        # Loop over instances of loaded dataset in parallel
        instance_start_idx <- ((j - 1) * (num_instances_per_cells_obj)) + 1
        instance_end_idx <- (num_instances_per_cells_obj * j)
        foreach (l = instance_start_idx:instance_end_idx) %dopar% {
            cpu_number <- (l %% num_cpus)
            prepended_fname <- paste0("thread_", cpu_number, "_",
                                      basename(out_file_name))
            out_file_name <- file.path(dirname(out_file_name), prepended_fname)
            instance <- test_data[l, ]
            if (dataset %in% c("sim_discrete_clusters", "sim_linear_traj",
                               "sim_branch_traj")) {
                upreg_pop <- paste0(as.character(instance[1]), "_pop_", instance[4])
                upreg_pc1 <- as.numeric(instance[2])
                rep <- as.numeric(instance[3])
            } else {
                upreg_pop <- as.character(instance[2])
                upreg_pc1 <- as.numeric(instance[3])
                rep <- as.numeric(instance[4])
            }

            labels <- as.numeric(instance[(num_extra_fields + 1):(num_extra_fields + num_cells)])
            labels[which(labels == 1)] <- "Condition1"
            labels[which(labels == 0)] <- "Condition2"
            cells@meta.data$synth_labels <- labels

            ada <- (as.numeric(instance[(num_extra_fields + num_cells + 1):(num_extra_fields + (2 * num_cells))]) == 1)
            pda <- (as.numeric(instance[(num_extra_fields + (2 * num_cells) + 1):length(instance)]) == 1)
            cells@meta.data$is_ada <- ada
            cells@meta.data$is_pda <- pda

            if (thing_to_bench %in% c("accuracy", "accuracy_verbose")) {
                res <- benchmark_accuracy(cells, instance, method,
                                          nn_model, upreg_pc1, upreg_pop,
                                          rep, milo_graph, out_file_name,
                                          methods_params, cpu_number,
                                          thing_to_bench, downsampling)
                write_to_file(res, out_file_name)
            } else if (thing_to_bench == "MSE") {
                res <- benchmark_mse(cells, instance, method,
                                     out_file_name, upreg_pc1, upreg_pop, rep,
                                     methods_params, cpu_number, nn_model = nn_model)
                write_to_file(res, out_file_name)
            } else if (thing_to_bench == "time") {
                for (downsample_size in seq(50000, 250000, 50000)) {
                    message(paste("Benchmarking runtime for downsample_size =",
                                  downsample_size))
                    for (rep in 1:3) {
                        res <- benchmark_time(cells, instance, method,
                                              nn_model, reduced_dim,
                                              out_file_name,
                                              methods_params,
                                              downsample_size, rep,
                                              num_cells)
                        write_to_file(res, out_file_name)
                    }
                }
            }
        }
        stopImplicitCluster()
    }

    # Combine output files
    prepended_fname <- paste0("thread_*_", basename(out_file_name))
    out_wildcarded <- file.path(dirname(out_file_name), prepended_fname)

    system2("cat", args = out_wildcarded, stdout = out_file_name)
}
