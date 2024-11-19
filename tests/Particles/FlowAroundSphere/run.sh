#!/bin/bash

set -o errexit
set -o pipefail

ACTION=$1

if [[ "$ACTION" == "test" ]]; then
    MGLET_BIN=$2
    mpirun -n 1 $MGLET_BIN 2>&1 | tee mglet.OUT
elif [[ "$ACTION" == "newtest" ]]; then
    rm -rf LOGS fields.h5 ib_stencils.h5 mglet-perf-report.txt *.OUT Particle_Snapshots Particle_Statistics
    MGLET_BIN=$2
    mpirun -n 1 $MGLET_BIN 2>&1 | tee mglet.OUT
elif [[ "$ACTION" == "clean" ]]; then
    rm -rf LOGS fields.h5 ib_stencils.h5 probes.h5 snapshots.h5 mglet-perf-report.txt *.OUT Particle_Snapshots Particle_Statistics
else
    echo "Invalid action: $ACTION"
    exit 1
fi
