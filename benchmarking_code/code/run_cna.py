import sys

import random
import time

import numpy as np
import scanpy as sc
import pandas as pd
import cna

import anndata as ad
from multianndata import MultiAnnData


def run_cna_py(input_file, reduced_dim, n_neighbors, thing_to_bench):

    random.seed(123)

    n_neighbors = int(n_neighbors)
    reduced_dim = "X_" + reduced_dim

    e = ad.read_h5ad(input_file)

    e.obs["synth_labels_int"] = [int(i[-1]) for i in e.obs.synth_labels]

    e.obs["id"] = 0
    num_samples_per_label = 10
    e.obs.loc[e.obs.synth_labels_int == 1, "id"] = random.choices(
        [(i + 1) for i in range(num_samples_per_label)],
        k=len(e.obs[e.obs.synth_labels_int == 1]),
    )
    e.obs.loc[e.obs.synth_labels_int == 2, "id"] = random.choices(
        [(i + 1) for i in range(num_samples_per_label, 2 * num_samples_per_label)],
        k=len(e.obs[e.obs.synth_labels_int == 2]),
    )

    e = MultiAnnData(e)
    e.obs_to_sample(["synth_labels_int"])

    start_time = time.time()

    # Find neighbors then compute UMAP. The n_neighbors parameter could
    # conceivably be optimised. It should probably be set to the value of k
    # used by Milo.
    sc.pp.neighbors(e, use_rep=reduced_dim, n_neighbors=n_neighbors)
    sc.tl.umap(e, min_dist=0.5)

    # By default, Nnull = 1000. However, we kept getting the warning:
    #   "UserWarning: global association p-value attained minimal possible \
    #       value. Consider increasing Nnull"
    # I have done so here to improve the performance of CNA.
    res = cna.tl.association(
        e, e.samplem.synth_labels_int, "id", Nnull=10000, key_added="cna_verdict"
    )

    stop_time = time.time()

    if thing_to_bench in ["accuracy", "accuracy_verbose"]:
        return e.obs.cna_verdict_fdr <= 0.1
    elif thing_to_bench == "MSE":
        return [((x + 1) / 2) for x in res.ncorrs]
    elif thing_to_bench == "time":
        return [stop_time - start_time]


h5_filename_base = sys.argv[1]
input_file = h5_filename_base + ".h5ad"
reduced_dim = sys.argv[2]
n_neighbors = int(sys.argv[3])
thing_to_bench = sys.argv[4]
res = run_cna_py(input_file, reduced_dim, n_neighbors, thing_to_bench)
output_file = h5_filename_base + "_cna_out.csv"
pd.DataFrame(res).to_csv(output_file, header=False, index=False)
