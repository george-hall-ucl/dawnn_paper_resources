# for SIMULATION_TYPE in "discrete_clusters" "linear_trajs" "branching_trajs" "skin" "mouse" "organoid" "heart"; do
# for SIMULATION_TYPE in "discrete_clusters" "linear_trajs" "branching_trajs"; do
for SIMULATION_TYPE in "mouse_mse"; do
    echo "Running pipeline for ${SIMULATION_TYPE}"
    sh code/run_simulation_pipeline.sh ${SIMULATION_TYPE} 0.05 > stdlogs/${SIMULATION_TYPE}.log 2>&1 &
done
