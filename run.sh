#!/usr/bin/env bash
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
CHAKRA_CONVERTER="chakra_converter"  # must be on PATH
MODEL="deepseek"
WORKLOAD_TXT="${MOE_DIR:?}/${MODEL:?}.txt"
WORKLOAD_CSV="${MOE_DIR:?}/${MODEL:?}.csv"
WORKLOAD_PREFIX="${WORKLOAD_DIR:?}/${MODEL:?}"
SYSTEM="${MOE_DIR:?}/system.json"
NETWORK="${MOE_DIR:?}/network.yml"
REMOTE_MEMORY="${MOE_DIR:?}/remote_memory.json"
LOG_OUTPUT_PREFIX="${MOE_DIR:?}/results/log"

NUM_NPUS=1024

usage() {
  cat <<EOF
Usage:
  $0 <action> [options]

Actions (exactly one):
  build_astra        [-l | -d]   Build Astra-Sim (debug if -l or -d present)
  build_chakra                   Build & install Chakra (pip install .)
  run_chakra         -b <N>      Launch N chakra runs with --index=0..N-1
  run_astra          -b <N>      Launch N astra runs with --index=0..N-1

Examples:
  $0 build_astra
  $0 build_astra -d
  $0 build_chakra
  $0 run_chakra -b 4
  $0 run_astra -b 8
EOF
}

[ $# -ge 1 ] || { usage; exit 1; }
ACTION="$1"; shift || true

case "$ACTION" in
  build_astra)
    # Options: -l or -d (pass through to build.sh)
    DEBUG_FLAG=""
    OPTIND=1
    while getopts "ld" opt; do
      case "$opt" in
        l) DEBUG_FLAG="-l" ;;
        d) DEBUG_FLAG="-d" ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage; exit 1 ;;
      esac
    done
    shift $((OPTIND - 1))

    echo "[ASTRA-sim] Building Astra-Sim (${DEBUG_FLAG:-release})..."
    cd "${PROJECT_DIR}"
    git submodule update --init --recursive
    ./build/astra_analytical/build.sh ${DEBUG_FLAG}
    echo "[ASTRA-sim] Build complete."
    ;;

  build_chakra)
    echo "[ASTRA-sim] Building & installing Chakra..."
    cd "${CHAKRA_DIR}"
    # Ensure ~/.local/bin is on PATH for user installs
    export PATH="$HOME/.local/bin:$PATH"
    pip3 install .
    echo "[ASTRA-sim] Chakra installed. (Tip: ensure \$HOME/.local/bin is in PATH)"
    ;;

  run_chakra)
    # Option: -b <start,end> (or start:end), inclusive
    B_RANGE=
    OPTIND=1
    while getopts ":b:" opt; do
      case "$opt" in
        b) B_RANGE="$OPTARG" ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage; exit 1 ;;
        :)  echo "Option -$OPTARG requires an argument." >&2; usage; exit 1 ;;
      esac
    done
    [ -n "$B_RANGE" ] || { echo "Missing -b <start,end> for run_chakra"; usage; exit 1; }

    IFS=',:' read -r B_START B_END <<< "$B_RANGE"
    [[ $B_START =~ ^[0-9]+$ && $B_END =~ ^[0-9]+$ ]] || { echo "Invalid -b format. Use -b start,end"; exit 1; }
    (( B_END >= B_START )) || { echo "Invalid -b: end < start"; exit 1; }

    : "${NUM_NPUS:=1024}"
    NUM_RUNS=$(( B_END - B_START + 1 ))
    echo "[Chakra] Launching batches ${B_START}..${B_END} (${NUM_RUNS} run(s))..."
    for (( i=B_START; i<=B_END; ++i )); do
      echo "[RUN][chakra][$i] ${CHAKRA_CONVERTER} Yaml --input=${WORKLOAD_CSV} --output=${WORKLOAD_PREFIX}_${i} --num-npus=${NUM_NPUS} --batch-id=${i}"
      "${CHAKRA_CONVERTER}" Yaml \
        --input="${WORKLOAD_CSV}" \
        --output="${WORKLOAD_PREFIX}_${i}" \
        --num-npus="${NUM_NPUS}" \
        --batch-id="${i}" &
    done
    wait
    echo "[Chakra] All runs finished."
    ;;

  run_astra)
    # Option: -b <start,end> (or start:end), inclusive
    B_RANGE=
    OPTIND=1
    while getopts ":b:" opt; do
      case "$opt" in
        b) B_RANGE="$OPTARG" ;;
        \?) echo "Unknown option: -$OPTARG" >&2; usage; exit 1 ;;
        :)  echo "Option -$OPTARG requires an argument." >&2; usage; exit 1 ;;
      esac
    done
    [ -n "$B_RANGE" ] || { echo "Missing -b <start,end> for run_astra"; usage; exit 1; }

    IFS=',:' read -r B_START B_END <<< "$B_RANGE"
    [[ $B_START =~ ^[0-9]+$ && $B_END =~ ^[0-9]+$ ]] || { echo "Invalid -b format. Use -b start,end"; exit 1; }
    (( B_END >= B_START )) || { echo "Invalid -b: end < start"; exit 1; }

    NUM_RUNS=$(( B_END - B_START + 1 ))
    echo "[ASTRA-sim] Launching batches ${B_START}..${B_END} (${NUM_RUNS} run(s))..."
    for (( i=B_START; i<=B_END; ++i )); do
      echo "[RUN][astra][$i] ${ASTRA_SIM} --workload-configuration=${WORKLOAD_PREFIX}_${i} --system-configuration=${SYSTEM} --remote-memory-configuration=${REMOTE_MEMORY} --network-configuration=${NETWORK} --log-output-path=${LOG_OUTPUT_PREFIX}_${i}.txt --index=${i}"
      "${ASTRA_SIM}" \
        --workload-configuration="${WORKLOAD_PREFIX}_${i}" \
        --system-configuration="${SYSTEM}" \
        --remote-memory-configuration="${REMOTE_MEMORY}" \
        --network-configuration="${NETWORK}" \
        --log-output-path="${LOG_OUTPUT_PREFIX}_${i}.txt" \
        --index="${i}" &
    done
    wait
    echo "[ASTRA-sim] All runs finished."
    ;;

  *)
    echo "Unknown action: ${ACTION}" >&2
    usage
    exit 1
    ;;
esac
