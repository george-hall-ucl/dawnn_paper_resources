This repository aims to ensure reproducibility of the results in the paper
_Dawnn: single-cell differential abundance with neural networks_. The notebook
to produce the figures in the paper can be found
[here](github.com/george-hall-ucl/dawnn_paper_results/blob/main/benchmarking_results.md).

The table below shows where to locate the data and code used to generate the
files read in the above notebook. Taken together, this information should allow
all results to be reproduced. If you are unable to reproduce any result, please
contact me at `george.hall@ucl.ac.uk` or post an issue on this repository.

| Directory name | Description | DOI |
| -------------- | ----------- | --- |
| heart\_dataset | Benchmarking dataset based on heart samples. <br> • process\_heart\_cells.R - Code to generate benchmarking dataset <br> • heart\_tissue\_cells.RDS - Generated benchmarking dataset <br> • heart\_barcodes.tsv.gz - Barcode list for raw data <br> • heart\_genes.tsv.gz - Barcode list for raw data <br> • heart\_expression\_matrix.mtx.gz - Expression matrix for raw data <br> • benchmark\_dataset\_heart\_data\_type\_labels.csv - Generated benchmarking dataset | 10.5522/04/22601260 |
| skin\_dataset | Benchmarking dataset based on skin cells. <br> • benchmark\_dataset\_skin.csv - Resulting benchmarking dataset <br> • simulate\_skin\_labels\_Rscript.R - R code to generate benchmarking dataset <br> • simulate\_skin\_labels\_bash.sh - Bash script to generate benchmarking dataset <br> • skin\_data\_end\_pipeline\_1458110522.rds - Input dataset | 10.5522/04/22607236 |
| organoid\_dataset | Benchmarking dataset based on bile duct organoid cells. | 10.5522/04/22612576 |
| mouse\_dataset | Benchmarking dataset based on mouse cells. | 10.5522/04/22614004 |
| discrete\_clusters\_dataset | Benchmarking dataset based on simulated discrete clusters. | 10.5522/04/22616590 |
| linear\_trajectory\_dataset | Benchmarking dataset based on simulated linear trajectories. | 10.5522/04/22616611 |
| branching\_trajectory\_dataset | Benchmarking dataset based on simulated branching trajectories. | 10.5522/04/22619851 |
| dawnn\_trained\_model | Trained neural network model needed to run Dawnn. <br> • final\_model\_dawnn\_rerun.h5 - Final trained Dawnn model | 10.5522/04/22241017 |
dawnn\_model\_training | train\_final\_model\_regen\_seed\_123\_job\_sub.sh - Job submission script to train final selected model <br> • train\_nn\_regen\_seed\_123.py - Python script to train final selected model | 10.5522/04/22633606
training\_set\_simulation | Code and resulting data when simulating Dawnn's training set. <br> • autogen4\_code.R - Code to generate training set <br> • labels\_df.csv - Resulting generated training set | 10.5522/04/22634200
model\_evaluations | Code and results from evaluating different models. <br> • nn\_model\_choice.py - Code for hyperparameter optimization <br> • model\_evaluations\_structure\_all\_nn\_results.txt - Results from neural network hyperparameter optimization <br> • eval\_rf\_svm.py - Code to evaluate random forests and support vector machines <br> • svm\_model\_evaluations.txt - Results from evaluation of support vector machines <br> • rf\_model\_evaluations.txt - Results from evaluation of random forests | 10.5522/04/22634416
benchmarking\_code\_and\_results | • collect\_results\_all\_sim\_dat.R - Code to run and collect benchmarking for simulated datasets <br> • tpr\_fdr\_results\_discrete\_clusters\_rerun.csv - Results from benchmarking on discrete clusters dataset <br> • tpr\_fdr\_results\_linear\_traj\_rerun.csv - Results from benchmarking on linear trajectory dataset <br> • tpr\_fdr\_results\_branch\_traj\_rerun.csv - Results from benchmarking on branching trajectory dataset <br> • collecting\_results\_mouse.sh - Code to run and collect benchmarking for mouse dataset <br> • tpr\_fdr\_results\_mouse\_regen.csv - Results from benchmarking on mouse dataset <br> • collecting\_results\_skin.sh - Code to run and collect benchmarking for skin dataset <br> • tpr\_fdr\_results\_skin\_regen.csv - Results from benchmarking on skin dataset <br> • collecting\_results\_organoid.sh - Code to run and collect benchmarking for organoid dataset <br> • tpr\_fdr\_results\_organoid\_regen.csv - Results from benchmarking on organoid dataset <br> • collecting\_results\_heart.sh - Code to run and collect benchmarking for heart dataset <br> • tpr\_fdr\_results\_heart\_regen.csv - Results from benchmarking on heart dataset <br> • benchmarking\_liver\_cirrhosis\_analysis.R - Code to run methods on cirrhotic liver dataset <br> • liver\_cirrhosis\_results\_rerun.csv - Results from running on cirrhotic liver dataset | 10.5522/04/22634470
