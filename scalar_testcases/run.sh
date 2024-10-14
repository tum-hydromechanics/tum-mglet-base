#!/bin/bash

set -o errexit
set -o pipefail

ACTION=$1
NUM_GRIDS=$2
NUM_CELLS=$3

DEFAULT_NUM_MPI_PROCS=1
DEFAULT_MGLET_BINARY="../../build/src/mglet"

# Check if ACTION is valid
if [[ "$ACTION" != "help" && "$ACTION" != "test" && "$ACTION" != "clean" ]]; then
    echo "Error: ACTION must be one of 'help', 'test', or 'clean'."
    echo "Check Usage example with $0 help"
    exit 1
fi

# Show example usage of this script
if [[ "$ACTION" == "help" ]]; then
    echo "Usage: $0 {test|clean} {NUM_GRIDS} {NUM_CELLS} [MPI_PROCS=$DEFAULT_NUM_MPI_PROCS] [MGLET_BINARY=$DEFAULT_MGLET_BINARY]"
    echo
    echo "    NUM_GRIDS          8 | 64 | 256 | 512 | 4096"
    echo "    NUM_CELLS          32| 64"
    echo "                       Test case may not exist for each combination"
    exit 1
fi

# Check if NUM_GRIDS or NUM_CELLS is not configured
if [[ -z "$NUM_GRIDS" || -z "$NUM_CELLS" ]]; then
    echo "Error: NUM_GRIDS and NUM_CELLS must be specified."
    echo "Check Usage example with $0 help"
    exit 1
fi

# cd into test directory
TEST="$NUM_GRIDS-grids-$NUM_CELLS-cells"
cd $TEST

# Run test or cleanup results in test directory
if [[ "$ACTION" == "test" ]]; then
    MPI_PROC=${4-$DEFAULT_NUM_MPI_PROCS}
    MGLET_BINARY="${5:-$DEFAULT_MGLET_BINARY}"
    echo "Running test $TEST"
    echo "Test directory: $( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
    echo "MGLET binary: $MGLET_BINARY"
    mpirun -N $MPI_PROC $MGLET_BINARY 2>&1 | tee mglet.OUT
elif [[ "$ACTION" == "clean" ]]; then
    echo "Cleaning up files for test $TEST"
    rm -rf LOGS fields.h5 mglet-perf-report.txt *.OUT
else
    echo "Invalid action $ACTION for test $TEST"
    exit 1
fi
