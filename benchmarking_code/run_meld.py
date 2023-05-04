import sys
import csv
import random

import numpy as np
import pandas as pd
import graphtools as gt

import meld

random.seed(42)

n_dims = 50
k = 15

data = []
with open("mouse_embryo_pca_corrected.csv") as f:
    header = f.readline().strip().split(',')
    csv_reader = csv.reader(f, delimiter = ',')
    for row in csv_reader:
        data.append(row)

data = np.array(data)
X = data[:, range(1, (n_dims + 1))].astype(float)

sample_label_files = ["benchmark_data/benchmark_dataset_Erythroid1.csv",
                      "benchmark_data/benchmark_dataset_Gut.csv",
                      "benchmark_data/benchmark_dataset_Somitic mesoderm.csv"]

for sample_label_file in sample_label_files:
    data = []
    with open(sample_label_file) as f:
        csv_reader = csv.reader(f, delimiter = ',')
        for row in csv_reader:
            data.append(row)
    data = np.array(data, dtype = '<U100')
    upreg_pc1 = data[:, 1].astype(float)
    rep = data[:, 2].astype(int)
    sample_labels = data[:, 3:64021]

    print("Constructing graph")
    graph = gt.Graph(X, knn = k)

    upreg_pc1_idxs = np.in1d(upreg_pc1, [0.7, 0.8, 0.9])
    rep_1_idxs = (rep == 1)
    sample_idxs = [(a and b) for (a,b) in zip(upreg_pc1_idxs, rep_1_idxs)]
    for sample_labelling in sample_labels[sample_idxs]:
        sample_labelling[sample_labelling == "1"] = "Condition1"
        sample_labelling[sample_labelling == "0"] = "Condition2"

        print("Running meld.MELD()")
        meld_op = meld.MELD()
        meld_op.graph = graph
                                                                              
        print("Running meld_op.transform()")
        meld_fit = meld_op.transform(sample_labels = np.array(sample_labelling))
        conditions=["Condition1", "Condition2"]
        mean_density = pd.DataFrame(
                np.zeros(shape=(meld_fit.shape[0], len(conditions))),
                index = meld_fit.index,
                columns = conditions,)

        for c in conditions:
          c_mean = meld_fit.loc[:,[c in x for x in meld_fit.columns]].mean(1)
          mean_density[c] = c_mean
                                                                              
        print("Running normalize_densities()")
        likelihoods = meld.utils.normalize_densities(mean_density)
                                                                              
        with open(f"meld_out.csv", 'a') as f:
            print(f"{sample_label_file}," +
                    ",".join(str(x) for x in likelihoods["Condition1"]),
                  file = f)

# from the output (in R):
# library(data.table)
# meld_out <- fread("meld_out.csv", header=F)
# meld_out <- meld_out[,2:64019]
# ground_truth <- fread("benchmark_dataset_Erythroid1_data_type_pc1s.csv", header=F)
# ground_truth <- ground_truth[,5:64022]
# for (i in 1:9) {print(mean(as.numeric(ground_truth[ceiling(i/3),]) - as.numeric(meld_out[i,])))}
