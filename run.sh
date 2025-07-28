#!/bin/bash

# compile astra sim: https://astra-sim.github.io/astra-sim-docs/getting-started/build.html
# git clone git@github.com:astra-sim/astra-sim.git
# ASTRA_SIM=$(realpath ./astra-sim)
# cd ${ASTRA_SIM}
# git submodule update --init --recursive
# ./build/astra_analytical/build.sh


# compile chakra: 
# cd astra-sim/extern/graph_frontend/chakra
# pip3 install .
# or maybe:
# echo "[ASTRA-sim] Installing Chakra..."
# echo ""
# "${PROJECT_DIR:?}"/utils/install_chakra.sh
# remember to add the install path e.g.,:    export PATH=$HOME/.local/bin:$PATH

# a built-in example at examples/network_analytical/run_network_analytical.sh
set -e

## ******************************************************************************
## This source code is licensed under the MIT license found in the
## LICENSE file in the root directory of this source tree.
##
## Copyright (c) 2024 Georgia Institute of Technology
## ******************************************************************************

# find the absolute path to this script
SCRIPT_DIR=$(dirname "$(realpath "$0")")
PROJECT_DIR=$SCRIPT_DIR
MOE_DIR="${PROJECT_DIR:?}/moe"
WORKLOAD_DIR="${MOE_DIR:?}/workload"
CHAKRA_DIR="${PROJECT_DIR:?}/extern/graph_frontend/chakra"

# paths
ASTRA_SIM="${PROJECT_DIR:?}/build/astra_analytical/build/bin/AstraSim_Analytical_Congestion_Aware"
# CHAKRA_CONVERTER="${CHAKRA_DIR:?}/src/converter/converter.py"
CHAKRA_CONVERTER="chakra_converter"
MODEL="deepseek"
WORKLOAD_TXT="${MOE_DIR:?}/${MODEL:?}.txt"
WORKLOAD_PREFIX="${WORKLOAD_DIR:?}/${MODEL:?}"
SYSTEM="${MOE_DIR:?}/system.json"
NETWORK="${MOE_DIR:?}/network.yml"
REMOTE_MEMORY="${MOE_DIR:?}/remote_memory.json"

echo "[ASTRA-sim] Generating Chakra trace..."
echo ""

${CHAKRA_CONVERTER} Text \
    --input="${WORKLOAD_TXT:?}" \
    --output="${WORKLOAD_PREFIX:?}"\
    --num-npus=8 \
    --num-passes=1

echo "[ASTRA-sim] Running ASTRA-sim Example with Analytical Network Backend..."
echo ""

# run ASTRA-sim
${ASTRA_SIM} \
    --workload-configuration="${WORKLOAD}" \
    --system-configuration="${SYSTEM:?}" \
    --remote-memory-configuration="${REMOTE_MEMORY:?}" \
    --network-configuration="${NETWORK:?}"

# finalize
echo ""
echo "[ASTRA-sim] Finished the execution."
