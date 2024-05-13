#!/bin/bash

set -o errexit
set -o pipefail

ACTION=$1
CASE=$2

USAGE_EXAMPLE="./run.sh <bench|clean> <case 8|64> [processes=1] [mglet-binary-path=../build/src/mglet]"

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
    NUM_PROC=${3-1}
    MGLET_BIN="../${4:-"../build/src/mglet"}"
    echo "Running benchmark $BENCHMARK"
    echo "Using mglet binary $MGLET_BIN"
    mpirun -N $NUM_PROC $MGLET_BIN 2>&1 | tee mglet.OUT
elif [[ "$ACTION" == "clean" ]]; then
    echo "Cleaning up files for benchmark $BENCHMARK"
    rm -rf LOGS fields.h5 mglet-perf-report.txt *.OUT
else
    echo "Invalid action: $ACTION"
    exit 1
fi
