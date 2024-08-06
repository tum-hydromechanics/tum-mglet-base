#!/bin/bash

set -o errexit
set -o pipefail

ACTION=$1

if [[ "$ACTION" == "test" ]]; then
    MGLET_BIN=$2
    mpirun -n 3 $MGLET_BIN 2>&1 | tee mglet.OUT
elif [[ "$ACTION" == "newtest" ]]; then
    rm -rf LOGS fields.h5 ib_stencils.h5 mglet-perf-report.txt *.OUT Particle_Snapshots
    MGLET_BIN=$2
    mpirun -n 3 $MGLET_BIN 2>&1 | tee mglet.OUT
elif [[ "$ACTION" == "clean" ]]; then
    rm -rf LOGS fields.h5 ib_stencils.h5 mglet-perf-report.txt *.OUT Particle_Snapshots
else
    echo "Invalid action: $ACTION"
    exit 1
fi
