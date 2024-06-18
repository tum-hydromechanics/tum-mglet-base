#!/bin/bash

set -o errexit
set -o pipefail

ACTION=$1
CASE=$2

DEFAULT_NUM_PROCESSES=1
DEFAULT_MGLET_BIN_PATH="../../build/src/mglet"

USAGE_EXAMPLE="./run.sh <bench|clean> <case 8|64> [processes=$DEFAULT_NUM_PROCESSES] [mglet-binary-path=$DEFAULT_MGLET_BIN_PATH]"

# Determine case (8 or 64 grids)
if [[ $CASE != 8 ]] && [[ $CASE != 64 ]]; then
    echo "Invalid argument: No benchmark for $CASE grids. [use 8 or 64]"
    echo "Usage example:"
    echo $USAGE_EXAMPLE
    exit 1
fi

# Set directory
BENCHMARK="$CASE-grids-32-cells"
cd $BENCHMARK

if [[ "$ACTION" == "bench" ]]; then
    NUM_PROC=${3-$DEFAULT_NUM_PROCESSES}
    MGLET_BIN="${4:-$DEFAULT_MGLET_BIN_PATH}"
    echo "Running benchmark $BENCHMARK"
    echo "Benchmark directory: $( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
    echo "MGLET binary: $MGLET_BIN"
    mpirun -N $NUM_PROC $MGLET_BIN 2>&1 | tee mglet.OUT
elif [[ "$ACTION" == "clean" ]]; then
    echo "Cleaning up files for benchmark $BENCHMARK"
    rm -rf LOGS fields.h5 mglet-perf-report.txt *.OUT
else
    echo "Invalid action: $ACTION"
    exit 1
fi
