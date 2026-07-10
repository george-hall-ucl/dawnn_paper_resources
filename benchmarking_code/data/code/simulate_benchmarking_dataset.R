suppressPackageStartupMessages({
library(SingleCellExperiment)
library(tibble)
library(dplyr)
library(igraph)
library(cydar)
library(pdist)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(ggsci)
library(umap)
library(scran)
library(scater)
library(irlba)
library(mvtnorm)
library(Rfast)
library(dyntoy)
library(optparse)
library(Seurat)
  })

.find_centroid <- function(X_emb, cluster_membership){
  cl.ixs <- split(1:nrow(X_emb), cluster_membership)  
  centroid_emb <- sapply(cl.ixs, function(x) colMeans(X_emb[x, , drop=FALSE]))
  centroid_emb
}

.scale_to_range <- function(x, min=1, max=10){
  ((x - min(x))/(max(x)-min(x)))*(max-min) + min
}


.logit <- function(x, a=1){
  1/(1+exp(- a * x))
}

# Creates random differentially expressed regions over a dataset for benchmarking.
add_synthetic_labels_pop <- function(sce, # SingleCellExperiment obj
                                     pop, pop_column="celltype",
                                     pop_enr = 0.7,
                                     redDim='pca.corrected', # embedding to use to simulate differential abundance
                                     n_conditions=2, # number of conditions to simulate
                                     n_replicates=3, # number of replicates per condition
                                     n_batches = 2, # number of technical batches per condition (at least 2 replicates per batch)
                                     condition_balance = 1, # the distribution of cells across conditions (1 = equal)
                                     m=2, # Fuzziness parameter (higher m, more fuzziness)
                                     a_logit=0.5, # logit parameter
                                     cap_enr=NULL,
                                     seed=42){
  
  # pop_sce = sce[,sce[[pop_column]]==pop]
  
  set.seed(seed)
  conditions = paste0("Condition", 1:n_conditions)
  
  X_emb = reducedDim(sce, redDim)
  
  ## Find cluster center
  cluster_membership = sce[[pop_column]]
  centroid_emb <- .find_centroid(X_emb, cluster_membership)

  ## Assign weight to each cell for each cluster center
  centroid_dist <- pdist(X_emb, t(centroid_emb))
  centroid_dist <- as.matrix(centroid_dist)
  
  w <- sapply(1:ncol(centroid_dist),  function(j) 
    sapply(1:nrow(centroid_dist), function(i) 
      1/sum(centroid_dist[i,j]/centroid_dist[i,])^(2/(m-1))
    ) 
  )
  colnames(w) <- colnames(centroid_emb)
  rownames(w) <- rownames(X_emb)
  w <- apply(scale(w), 2, .logit, a=a_logit)
  ## Normalize weights from enr_score to 0.5
  enr_scores <- rep(0.5, ncol(w)) ## Generate enrichment prob for each cluster
  # enr_scores <- runif(ncol(w)) ## Generate enrichment prob for each cluster
  names(enr_scores) <- colnames(w)
  if(length(pop_enr) == length(pop)){
    enr_scores[pop] <- pop_enr
  } else{
    # assume all pops have the same enrichment
    pop_enr <- rep(pop_enr, length(pop))
    enr_scores[pop] <- pop_enr
  }
  
  
  # altering the baseline probability can induce a skew towards a condition across _all_ cells
  enr_prob <- sapply(1:ncol(w), function(i) .scale_to_range(w[,i], min=0.5*condition_balance,
                                                            max=enr_scores[i]))
  colnames(enr_prob) <- colnames(centroid_emb)
  
  # need to integrate over these to get the condition probabilities
  # need to set relevant pops only, force the others to ~0.5
  prob_matrix <- enr_prob[,pop]
  if(is(prob_matrix, "matrix")){
    cond_probability <- rowMeans(prob_matrix)
    for(x in seq_along(pop)){
      cond_probability[sce[[pop_column]] == pop[x]] <- prob_matrix[sce[[pop_column]] == pop[x], pop[x]]
    }
  } else{
    cond_probability <- prob_matrix
  }
  
  ## Cap probabilities (to have the same number of DA cells w different maximum Fold Change)
  if (!is.null(cap_enr)) {
    cond_probability <- ifelse(cond_probability > cap_enr, cap_enr, cond_probability)
  }

  cond_probability <- sapply(cond_probability, function(x) {if (x <= 0.55) {0.5} else {x}})
  
  cond_probability = cbind(cond_probability, 1 - cond_probability)
  colnames(cond_probability) = conditions
  print(head(cond_probability, n = 20))
  print(tail(cond_probability, n = 20))

  # Generate labels for condition and replicates
  synth_labels <- sapply(1:nrow(cond_probability),  function(i) sample(colnames(cond_probability), size = 1, prob = cond_probability[i,]))
  
  replicates <- paste0("R", 1:n_replicates)
  batches <- sample(paste0("B", rep(1:n_batches, each=n_replicates)))
  
  synth_samples <- paste0(synth_labels, "_", replicates)
  if(n_batches > 1){
   names(batches) <- sort(unique(synth_samples))
  } else{
    batches <- rep("B1", length(unique(synth_samples)))
    names(batches) <- unique(synth_samples)
  }
  synth_batches <- batches[synth_samples]
  
  # Add synthetic labels and probabilities to colData
  colData(sce)[["synth_labels"]] <- synth_labels
  # colData(sce)[["synth_replicates"]] <- synth_replicates
  colData(sce)[["synth_samples"]] <- synth_samples
  colData(sce)[["synth_batches"]] <- synth_batches
  colnames(cond_probability) <- paste0(colnames(cond_probability), "_prob")
  colData(sce)[colnames(cond_probability)] <- cond_probability
  return(sce)
}


## Simple synthetic condition labelling based on cluster membership
add_synthetic_labels_by_cluster <- function(sce, # SingleCellExperiment obj
                                            pop, pop_column="celltype",
                                            pop_enr = 0.7,
                                            redDim='pca.corrected', # embedding to use to simulate differential abundance
                                            n_conditions=2, # number of conditions to simulate
                                            n_replicates=3, # number of replicates per condition
                                            n_batches = 2, # number of technical batches per condition (at least 2 replicates per batch)
                                            condition_balance = 1, # the distribution of cells across conditions (1 = equal)
                                            m=2, # Fuzziness parameter (higher m, more fuzziness)
                                            a_logit=0.5, # logit parameter
                                            cap_enr=NULL,
                                            seed=42){
  
  set.seed(seed)
  conditions = paste0("Condition", 1:n_conditions)
  
  ## Set prop != 0.5 for cells in pop
  cond_probability <- rep(0.5, ncol(sce))
  cond_probability[sce[[pop_column]] == pop] <- pop_enr
  cond_probability = cbind(cond_probability, 1 - cond_probability)
  colnames(cond_probability) = conditions
  
  # Generate labels for condition and replicates
  synth_labels <- sapply(1:nrow(cond_probability),  function(i) sample(colnames(cond_probability), size = 1, prob = cond_probability[i,]))
  
  replicates <- paste0("R", 1:n_replicates)
  batches <- sample(paste0("B", rep(1:n_batches, each=n_replicates)))
  
  synth_samples <- paste0(synth_labels, "_", replicates)
  if(n_batches > 1){
    names(batches) <- sort(unique(synth_samples))
  } else{
    batches <- rep("B1", length(unique(synth_samples)))
    names(batches) <- unique(synth_samples)
  }
  synth_batches <- batches[synth_samples]
  
  # Add synthetic labels and probabilities to colData
  colData(sce)[["synth_labels"]] <- synth_labels
  # colData(sce)[["synth_replicates"]] <- synth_replicates
  colData(sce)[["synth_samples"]] <- synth_samples
  colData(sce)[["synth_batches"]] <- synth_batches
  colnames(cond_probability) <- paste0(colnames(cond_probability), "_prob")
  colData(sce)[colnames(cond_probability)] <- cond_probability
  return(sce)
}

# Creates random differentially expressed regions over a dataset for benchmarking.
add_synthetic_labels <- function(sce, # SingleCellExperiment obj
                                 knn_graph, # for knn smoothing of probability values
                                 redDim='pca.corrected', # embedding to use to simulate differential abundance
                                 n_conditions=2, # number of conditions to simulate
                                 n_components=10, # number of components of embedding to use
                                 n_replicates=3, # number of replicates per condition
                                 n_batches = 2, # number of technical batches per condition (at least 2 replicates per batch)
                                 seed=42){
  data_embedding = reducedDim(sce, redDim)[,1:n_components]
  set.seed(seed)
  # embedding data must be mean-centered
  data_embedding = t(scale(t(data_embedding), scale=FALSE))
  
  # Randomly flip sign of each embedding dimension
  data_embedding = apply(data_embedding, 2, function(x)  x*sample(c(-1, 1), 1) )
  
  conditions = paste0("Condition", 1:n_conditions)
  cond_probability = sapply(1:(length(conditions)-1), function(x) .make_pdf(data_embedding))
  
  # KNN Smoothing to avoid regions of the graph with opposite labels
  cond_probability = .knn_smoothing(knn_graph, cond_probability, redDim=redDim, d=n_components)
  
  # Normalize to sum to 1 for each cell
  # cond_probability <- t(apply(cond_probability, 1, function(x) x/sum(abs(x))))
  cond_probability = cbind(cond_probability, 1 - rowSums(cond_probability))
  colnames(cond_probability) = conditions
  
  # Generate labels for condition and replicates
  synth_labels <- sapply(1:nrow(cond_probability),  function(i) sample(colnames(cond_probability), size = 1, prob = cond_probability[i,]))
  replicates <- paste0("R", 1:n_replicates)
  batches <- sample(paste0("B", rep(1:n_batches, each=n_replicates)))
  synth_samples <- paste0(synth_labels, "_", replicates)
  names(batches) <- sort(unique(synth_samples))
  synth_batches <- batches[synth_samples]
  
  # Add synthetic labels and probabilities to colData
  colData(sce)[["synth_labels"]] <- synth_labels
  # colData(sce)[["synth_replicates"]] <- synth_replicates
  colData(sce)[["synth_samples"]] <- synth_samples
  colData(sce)[["synth_batches"]] <- synth_batches
  colnames(cond_probability) <- paste0(colnames(cond_probability), "_prob")
  colData(sce)[colnames(cond_probability)] <- cond_probability
  return(sce)
}

simulate_linear_trajectory <- function(n.milestones = 2, num_cells = 5000,
                                       gex_seed = 42){
    set.seed(gex_seed)
    message("Running generate_dataset().")
    linear_model <- model_linear(num_milestones = n.milestones)
    dataset <- generate_dataset(model = linear_model, num_cells = num_cells,
                                num_features = 2000, add_velocity = FALSE)
    
    cnts <- t(dataset$counts)
    coldata <- data.frame(row.names = colnames(cnts),
                          dataset$prior_information$groups_id)
     
    message("Running SingleCellExperiment().")
    sce <- SingleCellExperiment(assays=list(counts = cnts,
                                            logcounts = t(dataset$expression)),
                                colData = coldata)

    message("Adding condition information.")
    branches <- dataset$prior_information$groups_id
    coldata_df <- data.frame(cell_id = colnames(sce))
    coldata_df <- left_join(coldata_df, branches)

    message("Simulating biological condition.")
    n_groups <- length(unique(branches$group_id))
    p_vec <- (1:n.milestones) / n.milestones
    if(length(p_vec) == 1){
        p_vec <- rep(p_vec, n.milestones)
    } else if (length(p_vec) != n.milestones) {
        stop("Condition proportions must == length(milestones)")
    }

    cells_label_C1 <- c()
    for (i in (1:n_groups)) {
      g <- paste0("M", i)
      p <- p_vec[i]
      cells_group_M_label_C1 <- sample(coldata_df$cell_id[coldata_df$group_id == g],
                                      size = floor(sum(coldata_df$group_id == g) * p))
      cells_label_C1 <- c(cells_label_C1, cells_group_M_label_C1)
    }

    cell_id_in_C1 <- coldata_df$cell_id %in% cells_label_C1
    coldata_df <- coldata_df %>%
                    dplyr::mutate(Condition = ifelse(cell_id_in_C1,
                                                     'Condition1',
                                                     'Condition2'))

    message("Simulating replicates.")
    coldata_df <- coldata_df %>%
                    group_by(group_id) %>%
                    dplyr::mutate(Replicate = c(rep("R1", floor(n() * 0.3)), 
                                                rep("R2", floor(n() * 0.3)), 
                                                rep("R3", n() - 2 * (floor(n() * 0.3))))
                    )
    coldata_df$Sample <- paste(coldata_df$Condition, coldata_df$Replicate,
                               sep = "_")

    sim.seu <- CreateSeuratObject(counts = sce@assays@data$logcounts) %>%
                NormalizeData() %>%
                FindVariableFeatures() %>%
                ScaleData() %>%
                RunPCA() %>%
                RunUMAP(1:20)

    sim.seu@meta.data$group_id <- coldata_df$group_id
    group_id_as_num <- as.numeric(sub(".", "", sim.seu@meta.data$group_id))
    sim.seu@meta.data$Condition1_prob <- p_vec[group_id_as_num]
    sim.seu@meta.data$synth_labels <- coldata_df$Condition
    sim.seu@meta.data$synth_samples <- coldata_df$Sample
    sim.seu@meta.data$synth_labels[which(sim.seu@meta.data$synth_labels == "A")] <- "Condition1"
    sim.seu@meta.data$synth_labels[which(sim.seu@meta.data$synth_labels == "B")] <- "Condition2"

    return(sim.seu)
}


simulate_branching_trajectory <- function(n.milestones = 1, num_cells = 7500,
                                          gex_seed = 42){
    dataset <- generate_dataset(model = model_tree(num_branchpoints = 1,
                                                   max_degree = 10),
                                num_cells = 7500, num_features = 2000,
                                add_prior_information = FALSE,
                                add_velocity = FALSE)
    cells <- CreateSeuratObject(counts=t(dataset$counts))
    cells <- cells %>%
                NormalizeData() %>%
                FindVariableFeatures() %>%
                ScaleData() %>%
                RunPCA() %>%
                RunUMAP(1:20)

    return(cells)
}


sim_discrete_clusters_gex <- function(rand_seed = 2) {
    n.clusters = 3
    total.size = 2700
    cells.per.cluster=c(900, 900, 900)

    set.seed(rand_seed)

    r.n <- 1000
    n.dim <- 50
    if(n.clusters != length(cells.per.cluster)){
        stop("Number of clusters must equal the length of cells.per.cluster")
    }

    if(total.size != sum(cells.per.cluster)){
        stop("Total sample size must match the sum of cells.per.cluster")
    }

    gex.list <- list()
    for(x in seq_along(1:n.clusters)){
        block.cells <- cells.per.cluster[x]
        # select a set of eigen values for the covariance matrix of each block, say 50 eigenvalues?
        # randomly sample a mean value
        block.mean <- runif(n=1, min=2, max=30)
        block.eigens <- sapply(1:n.dim, FUN=function(X) rexp(n=1, rate=abs(runif(n=1, min=0, max=50))))
        block.eigens <- block.eigens[order(block.eigens)]
        block.p <- qr.Q(qr(matrix(rnorm(block.cells^2, mean=4, sd=0.01), block.cells)))
        block.sigma <- crossprod(block.p*block.eigens, block.p*block.eigens)

        # The below lines of prevent a weird bug from appearing when run in
        # docker, possibly by loading some things into memory earlier than the
        # otherwise would be? Who knows! But in any case, it works and
        # presumably gives the same result as if they weren't here.
        for (nnew in c(130, 150, 200, 300, 500, 900)) {
            block.sigma <- crossprod((block.p[1:nnew, 1:nnew]) * block.eigens,
                                     (block.p[1:nnew, 1:nnew]) * block.eigens)
        }

        block.gex <- abs(Rfast::rmvnorm(n = r.n,
                                        mu = rnorm(n = block.cells,
                                                   mean = block.mean,
                                                   sd = 0.01),
                                        sigma = block.sigma))
        gex.list[[paste0("Block", x)]] <- block.gex
    }

    sim.gex <- do.call(cbind, gex.list)
    colnames(sim.gex) <- paste0("Cell", 1:ncol(sim.gex))
    rownames(sim.gex) <- paste0("Gene", 1:nrow(sim.gex))
    sim.pca <- prcomp_irlba(t(sim.gex), n = 50, scale. = TRUE, center = TRUE)

    sim.sce <- SingleCellExperiment(assays = list(logcounts = sim.gex),
                                    reducedDims = list("PCA" = sim.pca$x))
    sim.seu <- CreateSeuratObject(counts = sim.sce@assays@data$logcounts, data = sim.sce@assays@data$logcounts)
    sim.seu <- sim.seu %>%
                FindVariableFeatures() %>%
                ScaleData() %>%
                RunPCA() %>%
                RunUMAP(1:20)
    sim.seu@misc$cells.per.cluster <- cells.per.cluster

    return(sim.seu)
}


sim_discrete_clusters_labels <- function(sim.seu, rand_seed, max_pc1) {

    cells.per.cluster <- sim.seu@misc$cells.per.cluster
    condition.props = c(max_pc1, 0.5, 0.5)
    if(length(condition.props) != length(cells.per.cluster)) {
        stop("The length of condition.props must be the same as the length of cells.per.cluster")
    }

    set.seed(rand_seed)
    cond.list <- list()
    cell.list <- list()
    for(i in seq_along(condition.props)){
        block.cells <- cells.per.cluster[i]
        cell.list[[paste0("Block", i)]] <- block.cells
        
        block.cond <- rep("A", block.cells)
        block.a <- sample(1:block.cells,
                          size = floor(block.cells * condition.props[i]))
        block.b <- setdiff(1:block.cells, block.a)
        block.cond[block.b] <- "B"
        cond.list[[paste0("Block", i)]] <- block.cond
    }
    blocks <- lapply(c(1:length(cell.list)),
                     FUN = function(X) rep(paste0("B", X), cell.list[[X]]))
    rep.prop <- round(1 / length(cells.per.cluster), 2)
    reps <- lapply(c(1:length(cell.list)),
                   FUN = function(X) c(rep("R1", floor(cell.list[[X]] * rep.prop)),
                                       rep("R2", floor(cell.list[[X]] * rep.prop)),
                                       rep("R3", cell.list[[X]] - (2 * floor(cell.list[[X]] * rep.prop)))))

    meta.df <- data.frame("Block" = unlist(blocks),
                          "Condition" = unlist(cond.list),
                          "Replicate" = unlist(reps))
    colnames(meta.df) <- c("Block", "Condition", "Replicate")
    rownames(meta.df) <- paste0("Cell", 1:nrow(meta.df))
    # define a "sample" as the combination of condition and replicate
    meta.df$Sample <- paste(meta.df$Condition, meta.df$Replicate, sep = "_")
    meta.df$Vertex <- c(1:nrow(meta.df))

    sim.seu@meta.data$synth_labels <- meta.df$Condition
    sim.seu@meta.data$synth_samples <- meta.df$Sample
    sim.seu@meta.data$synth_labels[which(sim.seu@meta.data$synth_labels == "A")] <- "Condition1"
    sim.seu@meta.data$synth_labels[which(sim.seu@meta.data$synth_labels == "B")] <- "Condition2"
    l <- c()
    for (i in 1:length(cells.per.cluster)) {
        l <- c(l, rep(condition.props[i], cells.per.cluster[i]))
    }
    sim.seu@meta.data$Condition1_prob <- l

    return(sim.seu)
}

add_synth_labels <- function(cells, upreg_pop, pop_enr, reduced_dim = "PCA",
                             cell_labels = "seurat_clusters", method = "pop",
                             seed = 42, m = 2) {
    if (method == "pop") {
        synthed <- as.Seurat(add_synthetic_labels_pop(as.SingleCellExperiment(cells),
                                                      pop = upreg_pop,
                                                      pop_enr = pop_enr,
                                                      redDim = reduced_dim,
                                                      pop_column = cell_labels,
                                                      seed = seed, m = m))
        cells@meta.data <- synthed@meta.data
    }
    else if (method == "cluster") {
        synthed <- as.Seurat(add_synthetic_labels_by_cluster(as.SingleCellExperiment(cells),
                                                      pop = upreg_pop,
                                                      pop_enr = pop_enr,
                                                      redDim = reduced_dim,
                                                      pop_column = cell_labels, seed = seed, m = m))
        cells@meta.data <- synthed@meta.data
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
            out_dfs <- list()
            out_dfs[["labels"]] <- data.frame(t(c(c(upreg_pop, pc1, rep),
                                                   as.numeric(cells@meta.data$synth_labels == "Condition1"),
                                                   as.numeric(cells@meta.data$is_da))))
            if (data_type == "pc1s") {
                out_dfs[["pc1s"]] <- data.frame(t(c(c(upreg_pop, pc1, rep),
                                                    cells@meta.data$Condition1_prob,
                                                    cells@meta.data$synth_labels)))
            }
            for (data_type in names(out_dfs)) {
                file_name <- paste0("benchmark_dataset_", upreg_pop, "_data_type_",
                                    data_type, ".csv")
                write.table(out_dfs[[data_type]], file_name, append = TRUE, sep = ",",
                            col.names = FALSE, row.names = FALSE)
            }
        }
    }
}


generate_benchmarking_data <- function(name, rds_file, reduced_dim, cell_labels, epsilon, data_type = "labels") {
    cells <- readRDS(rds_file)

    if (name == "skin") {
        cells <- FindNeighbors(cells, dims = 1:50, reduction = "harmony") %>% FindClusters()
        upreg_pops <- c("1", "3", "4", "5", "7", "8", "10", "11")
    } else if (name == "organoid") {
        upreg_pops <- c("0", "3", "4", "6", "9", "10", "12", "17")
    } else if (name == "heart") {
        upreg_pops <- c("1", "9", "14", "20", "22", "24", "29", "31")
    } else if (name %in% c("mouse", "mouse_mse")) {
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
                cells <- add_ada_pda_labels(cells, epsilon)
                out_dfs <- list()
                out_dfs[["labels"]] <- data.frame(t(c(c(rds_file, upreg_pop, pc1, rep),
                                                      as.numeric(cells@meta.data$synth_labels == "Condition1"),
                                                      as.numeric(cells@meta.data$is_ada),
                                                      as.numeric(cells@meta.data$is_pda))))
                if (data_type == "pc1s") {
                    out_dfs[["pc1s"]] <- data.frame(t(c(c(upreg_pop, pc1, rep),
                                                        cells@meta.data$Condition1_prob,
                                                        cells@meta.data$synth_labels)))
                }
                for (data_type in names(out_dfs)) {
                    file_name <- paste0("outputs/", name, "/benchmark_dataset_",
                                        name, "_", upreg_pop, "_data_type_",
                                        data_type, ".csv")
                    write.table(out_dfs[[data_type]], file_name, append = TRUE, sep = ",",
                                col.names = FALSE, row.names = FALSE)
                }
            }
        }
    }
}

add_ada_pda_labels <- function(cells, epsilon) {
    cells@meta.data$is_ada <- (cells@meta.data$Condition1_prob > 0.5 + epsilon) |
                              (cells@meta.data$Condition1_prob < 0.5 - epsilon)
    c1_prop <- sum(cells@meta.data$synth_labels == "Condition1") / ncol(cells)
    cells@meta.data$is_pda <- (cells@meta.data$Condition1_prob > c1_prop + epsilon) |
                              (cells@meta.data$Condition1_prob < c1_prop - epsilon)

    return(cells)
}



# Main code:

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
simulation_type <- args[1]
epsilon <- as.numeric(args[2])

if (simulation_type %in% c("discrete_clusters", "linear_trajs", "branching_trajs")) {
    output_dir <- paste0("outputs/", simulation_type, "/")
    for (gex_seed in 1:6) {
        message(paste("Simulating", simulation_type, "with gex_seed:", gex_seed))
        if (simulation_type == "discrete_clusters") {
            cells <- sim_discrete_clusters_gex(rand_seed = gex_seed)
        } else if (simulation_type == "linear_trajs") {
            cells <- simulate_linear_trajectory(n.milestones = 3,
                                                num_cells = 10000,
                                                gex_seed = gex_seed)
        } else if (simulation_type == "branching_trajs") {
            cells <- simulate_branching_trajectory(n.milestones = 1,
                                                   num_cells = 7500,
                                                   gex_seed = gex_seed)
        } else {
            stop(paste("Unknown simulation_type", simulation_type))
        }

        rds_file_name <- paste0(output_dir, "cells_sim_", simulation_type,
                                "_gex_seed_", gex_seed, ".rds")
        saveRDS(cells, rds_file_name)

        if (simulation_type != "discrete_clusters") {
            cells <- FindNeighbors(cells, reduction = "pca", dims = 1:20)
            cells <- FindClusters(cells, resolution = 0.2)
            group_ids <- unique(cells@meta.data$seurat_clusters)
        }

        for (pc1 in seq(0.7, 0.95, 0.05)) {
            for (rep in 1:3) {
                # Create unique number with each combination of (pc1, gex_seed, rep)
                seed <- (pc1 * 10000) + (10 * gex_seed) + rep
                set.seed(seed)
                if (simulation_type == "discrete_clusters") {
                    cells <- sim_discrete_clusters_labels(cells, rand_seed = seed,
                                                          max_pc1 = pc1)
                    upreg_pop <- 0
                    upreg_prop <- 1/3
                } else {
                    upreg_pop <- sample(group_ids, 1)
                    cell_labels <- "seurat_clusters"
                    cells <- add_synth_labels(cells, upreg_pop = upreg_pop,
                                              pop_enr = pc1, cell_labels = cell_labels,
                                              seed = seed)
                    upreg_prop <- prop.table(table(cells[[cell_labels]]))[upreg_pop]
                }
                cells <- add_ada_pda_labels(cells, epsilon = epsilon)

                run_info <- c(paste0("gex_seed_", gex_seed), pc1, rep, upreg_pop, upreg_prop)
                c1_binary <- as.numeric(cells@meta.data$synth_labels == "Condition1")
                ada_binary <- as.numeric(cells@meta.data$is_ada)
                pda_binary <- as.numeric(cells@meta.data$is_pda)

                labels_df <- data.frame(t(c(run_info, c1_binary, ada_binary, pda_binary)))

                write.table(labels_df,
                            paste0(output_dir, "benchmark_dataset_sim_",
                                   simulation_type, ".csv"),
                            append = TRUE, sep = ",", col.names = FALSE,
                            row.names = FALSE)
            }
        }
    }
} else if (simulation_type %in% c("heart", "mouse", "mouse_mse", "organoid", "skin")) {
    input_data_path <- "input_data/"

    dataset_info <- list()
    dataset_info[["skin"]] <- list(rds_file = paste0(input_data_path,
                                                     "skin_data_end_pipeline_1458110522.rds"),
                                   reduced_dim = "HARMONY",
                                   cell_labels = "seurat_clusters")
    dataset_info[["mouse"]] <- list(rds_file = paste0(input_data_path,
                                                      "embryo_data_bm_1604210222.RDS"),
                                   reduced_dim = "PCA.CORRECTED",
                                   cell_labels = "celltype")
    dataset_info[["mouse_mse"]] <- dataset_info[["mouse"]]
    dataset_info[["organoid"]] <- list(rds_file = paste0(input_data_path,
                                                         "organoid_tissue_cells.RDS"),
                                   reduced_dim = "PCA",
                                   cell_labels = "seurat_clusters")
    dataset_info[["heart"]] <- list(rds_file = paste0(input_data_path, "heart_tissue_cells.RDS"),
                                   reduced_dim = "PCA",
                                   cell_labels = "seurat_clusters")

    ds <- dataset_info[[simulation_type]]
    generate_benchmarking_data(name = simulation_type,
                               rds_file = ds$rds_file,
                               reduced_dim = ds$reduced_dim,
                               cell_labels = ds$cell_labels,
                               epsilon = as.numeric(epsilon),
                               data_type = ifelse(simulation_type != "mouse_mse", "labels", "pc1s"))

} else {
    stop(paste("Unknown simulation_type:", simulation_type))
}
