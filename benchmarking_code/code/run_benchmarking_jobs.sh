# run inside docker with something like:
# docker run --interactive --tty --rm            --volume "$(pwd):/tmp/in_mnt"            --workdir /tmp/in_mnt/code            --oom-score-adj -1000            ${DOCKER_IMAGE_NAME} bash run_benchmarking_jobs.sh <DATASET NAME>

DATASET=$1

eval "$(conda shell.bash hook)"
conda activate r_env

# for METHOD in "dawnn_ada" "daseq" "meld"; do
for METHOD in "meld"; do

    if [[ ${METHOD} == dawnn* ]]; then
        NUM_CPU=1
    else
        NUM_CPU=2
    fi

    echo "Running benchmarking for ${METHOD}"
    Rscript benchmarking_utilities_as_script.R MSE ${METHOD} ${DATASET} ../results/${DATASET}/${METHOD}.csv ${NUM_CPU}
done

conda deactivate
