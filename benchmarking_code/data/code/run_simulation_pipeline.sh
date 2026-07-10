# Run from top directory (i.e. containing `code` and `outputs` directories etc)

if [ $# != 2 ]; then
    echo "usage: $0 <SIMULATION_TYPE> <EPSILON>"
    exit 1
fi

SIMULATION_TYPE=$1
EPSILON=$2

docker run \
    --volume ${PWD}:/tmp/in_mnt \
    --volume /var/run/docker.sock:/var/run/docker.sock \
    --workdir /tmp/in_mnt \
    --rm \
    --name ${SIMULATION_TYPE} \
    georgehallucl/dawnn_dataset_generation \
    Rscript code/simulate_benchmarking_dataset.R ${SIMULATION_TYPE} ${EPSILON}

if [ "$SIMULATION_TYPE" != "discrete_clusters" ] && \
   [ "$SIMULATION_TYPE" != "linear_trajs" ] && \
   [ "$SIMULATION_TYPE" != "branching_trajs" ]; then
    OUTPUT_PATH="outputs/${SIMULATION_TYPE}"
    cat ${OUTPUT_PATH}/benchmark_dataset_${SIMULATION_TYPE}_*_data_type_labels.csv > ${OUTPUT_PATH}/benchmark_dataset_${SIMULATION_TYPE}.csv
fi
