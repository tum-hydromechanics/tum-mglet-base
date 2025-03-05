#!/bin/bash

set -o errexit
set -o pipefail

ACTION=$1

if [[ "$ACTION" == "test" ]]; then
    MGLET_BIN=$2
    cp parameters_start.json parameters.json
    mpirun -n 32 $MGLET_BIN 2>&1 | tee mglet_01.OUT
    mv fields.h5 fields_started.h5
    cp parameters_continue.json parameters.json
    mpirun -n 32 $MGLET_BIN 2>&1 | tee mglet_02.OUT
    rm parameters.json
elif [[ "$ACTION" == "clean" ]]; then
    rm -rf LOGS fields.h5 mglet-perf-report.txt *.OUT *.stl
else
    echo "Invalid action: $ACTION"
    exit 1
fi
