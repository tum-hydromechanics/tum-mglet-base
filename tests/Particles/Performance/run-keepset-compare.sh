#!/usr/bin/env bash
# Optional clean remeasure of the kept optimization set:
#   - physics guards (already in HEAD)
#   - RNG rewrite (working tree)
#   - lincon split reverted
#
# Compares against authoritative medium baselines. User-run timed tiers only.

source /etc/profile
module load mglet/gcc13 || {
    echo "Failed to load mglet/gcc13" >&2
    exit 1
}

set -Eeuo pipefail

ROOT=/home/yaydin/particle-performance
REPO="$ROOT/tum-mglet-base"
BUILD="$ROOT/build-particle-release"
RESULTS="$ROOT/benchmark-results"
BINARY="$BUILD/src/mglet"
HARNESS=tests/Particles/Performance/particle_bench.py
DRIVER_LOG="$RESULTS/keepset-compare-driver.log"

# Allow launch from repo root or from this directory.
cd "$REPO"
[[ -f "$HARNESS" ]] || {
    echo "Missing harness at $REPO/$HARNESS" >&2
    exit 1
}

mkdir -p "$RESULTS"
: > "$DRIVER_LOG"

log() {
    echo "$@" | tee -a "$DRIVER_LOG"
}

run_logged() {
    "$@" > >(tee -a "$DRIVER_LOG") 2>&1
}

log "Commit: $(git rev-parse HEAD)"
log "Dirty: $(git status --porcelain | wc -l) paths"
log "Binary: $BINARY"
log "Started: $(date -Is)"
log "Keep-set: guards + RNG; lincon split reverted"

if [[ "${SKIP_REBUILD:-0}" != "1" ]]; then
    log "Rebuilding Release mglet..."
    run_logged cmake --build "$BUILD" --target mglet -j 4
fi

[[ -x "$BINARY" ]] || {
    echo "Missing executable: $BINARY" >&2
    exit 1
}

run_logged env PYTHONDONTWRITEBYTECODE=1 python3 tests/Particles/Performance/test_particle_bench.py
run_logged env PYTHONDONTWRITEBYTECODE=1 python3 tests/Particles/Performance/test_gaussian_rng.py

compare_session() {
    local name=$1
    local baseline=$2
    local session=$3
    local comparison="$RESULTS/${name}-keepset-comparison.json"

    [[ -f "$session/summary.json" ]] || {
        log "Missing summary for $name: $session/summary.json"
        return 1
    }

    run_logged python3 "$HARNESS" compare \
        --result "$session/summary.json" \
        --baseline "$baseline" \
        --output "$comparison"

    log "$name session: $session"
    log "$name comparison: $comparison"
    if [[ -f "$session/summary.txt" ]]; then
        log "----- $name summary.txt -----"
        run_logged cat "$session/summary.txt"
    fi
}

run_and_compare() {
    local name=$1
    local baseline=$2
    shift 2

    local run_log="$RESULTS/${name}-keepset.log"
    local session rc

    set +e
    python3 "$HARNESS" run \
        --binary "$BINARY" \
        "$@" \
        --mpirun-arg=--bind-to \
        --mpirun-arg=core > "$run_log" 2>&1
    rc=$?
    set -e
    cat "$run_log" >> "$DRIVER_LOG"
    cat "$run_log"

    if (( rc != 0 )); then
        log "$name failed; see $run_log"
        return "$rc"
    fi

    session=$(
        awk -F'Results: ' '/^Results: / { result=$2 } END { print result }' "$run_log"
    )
    [[ -n "$session" && -f "$session/summary.json" ]] || {
        log "Could not locate summary for $name"
        return 1
    }

    compare_session "$name" "$baseline" "$session"
}

if [[ "${SKIP_SELF_CONTAINED:-0}" == "1" ]]; then
    log "Skipping self-contained medium suite (SKIP_SELF_CONTAINED=1)"
elif [[ -n "${SELF_CONTAINED_SESSION:-}" ]]; then
    log "Reusing self-contained session: $SELF_CONTAINED_SESSION"
    compare_session \
        self-contained-medium \
        "$RESULTS/self-contained-medium-baseline.json" \
        "$SELF_CONTAINED_SESSION"
else
    run_and_compare \
        self-contained-medium \
        "$RESULTS/self-contained-medium-baseline.json" \
        --tier medium
fi

run_and_compare \
    bcc-compute-medium \
    "$RESULTS/bcc-compute-medium-baseline.json" \
    --manifest tests/Particles/Performance/bcc-benchmarks.json \
    --tier medium \
    --case bcc-compute-medium

log "Finished: $(date -Is)"
log "Keep-set comparison completed."
log "Return these files:"
log "  $DRIVER_LOG"
log "  $RESULTS/self-contained-medium-keepset-comparison.json"
log "  $RESULTS/bcc-compute-medium-keepset-comparison.json"
