#!/bin/bash
set -e
set -x

## ******************************************************************************
## This source code is licensed under the MIT license found in the
## LICENSE file in the root directory of this source tree.
##
## Copyright (c) 2024 Georgia Institute of Technology
## ******************************************************************************

# find the absolute path to this script
SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_DIR="${SCRIPT_DIR:?}/../.."
EXAMPLE_DIR="${PROJECT_DIR:?}/examples/custom_collective"

# paths
ASTRA_SIM="${PROJECT_DIR:?}/build/astra_analytical/build/bin/AstraSim_Analytical_Congestion_Aware"
WORKLOAD="${EXAMPLE_DIR:?}/../network_analytical/workload/AllReduce_1MB"
SYSTEM="${EXAMPLE_DIR:?}/custom_collective.json"
NETWORK="${EXAMPLE_DIR:?}/../network_analytical/network.yml"
REMOTE_MEMORY="${EXAMPLE_DIR:?}/../network_analytical/remote_memory.json"

# start
echo "[ASTRA-sim] Compiling ASTRA-sim with the Analytical Network Backend..."
echo ""

# Compile
"${PROJECT_DIR:?}"/build/astra_analytical/build.sh

echo "[ASTRA-sim] Compilation finished."
echo "[ASTRA-sim] Running ASTRA-sim Example with Analytical Network Backend..."

# run ASTRA-sim
"${ASTRA_SIM:?}" \
    --workload-configuration="${WORKLOAD}" \
    --system-configuration="${SYSTEM:?}" \
    --remote-memory-configuration="${REMOTE_MEMORY:?}" \
    --network-configuration="${NETWORK:?}"

# finalize
echo "[ASTRA-sim] Finished the execution."
