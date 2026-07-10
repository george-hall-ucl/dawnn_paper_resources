source("benchmarking_utilities.R")

options(error = function() {
  calls <- sys.calls()
  if (length(calls) >= 2L) {
    sink(stderr())
    on.exit(sink(NULL))
    cat("Backtrace:\n")
    calls <- rev(calls[-length(calls)])
    for (i in seq_along(calls)) {
      cat(i, ": ", deparse(calls[[i]], nlines = 1L), "\n", sep = "")
    }
  }
  if (!interactive()) {
    q(status = 1)
  }
})

check_args <- function(thing_to_bench, method, dataset, downsampling) {
    if (!(thing_to_bench %in% c("accuracy", "MSE", "time", "accuracy_verbose"))) {
        stop(paste("Invalid thing_to_bench:", thing_to_bench))
    }

    valid_methods <- c("dawnn", "dawnn_pda", "dawnn_ada", "milo", "daseq",
                       "cna", "meld")
    if (!(method %in% valid_methods)) {
        stop(paste("Invalid method:", method))
    }

    if (!(dataset %in% c("mouse_embryo", "mouse_embryo_mse", "skin", "organoid", "heart", "sim_discrete_clusters", "sim_linear_traj", "sim_branch_traj", "fully_sim"))) {
        stop(paste("Invalid dataset:", dataset))
    }

    if (!(downsampling %in% c("TRUE", "FALSE"))) {
        stop(paste("Invalid downsampling:", downsampling))
    }
}

args <- commandArgs(trailingOnly = TRUE)
if ((length(args) == 0) | ("--help" %in% args)) {
    message("usage: Rscript benchmarking.R <thing_to_bench> <method> <dataset> <in_file_name> <out_file_name> <num_cpus> <[optional] cell_type>")
    quit()
}

thing_to_bench = args[1]
method = args[2]
dataset = args[3]
out_file_name = args[4]
num_cpus = args[5]
if (length(args) == 6) {
    downsampling = args[6]
} else {
    # Default to no downsampling
    downsampling <- "FALSE"
}

check_args(thing_to_bench, method, dataset, downsampling)
num_cpus <- as.numeric(num_cpus)
downsampling <- as.logical(downsampling)

print(paste("run as:", paste(args, collapse = " ")))
print(paste("thing_to_bench:", thing_to_bench))
print(paste("method:", method))
print(paste("dataset:", dataset))
print(paste("out_file_name:", out_file_name))
print(paste("num_cpus:", num_cpus))
print(paste("downsampling:", downsampling))

run_benchmarking(method = method,
                 thing_to_bench = thing_to_bench, dataset = dataset,
                 out_file_name = out_file_name, num_cpus = num_cpus, downsampling = downsampling)
