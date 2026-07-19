# Particle performance suite

This directory contains a local benchmark harness for measuring particle
optimizations without using machine-dependent CI timing thresholds.

The harness always runs correctness checks before timed cases unless
`--skip-correctness` is given. It stages inputs and writes all generated files
under `/home/yaydin/particle-performance/benchmark-results`; checked-in test
cases are never modified.

## Prerequisites

Load the server toolchain before configuring, building, or running:

```bash
source /etc/profile
module load mglet/gcc13
```

The correctness checks require the Python packages `h5py` and `numpy`, both of
which are available in the current server Python environment.

For this server, select the GCC 13 executables explicitly. The current
non-OpenMP source also calls a few OpenMP runtime query functions, so the
executable needs the OpenMP runtime at link time:

```bash
CC=gcc-13 CXX=g++-13 FC=gfortran-13 cmake \
  -S /home/yaydin/particle-performance/tum-mglet-base \
  -B /home/yaydin/particle-performance/build-particle-release \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_EXE_LINKER_FLAGS=-fopenmp

cmake --build /home/yaydin/particle-performance/build-particle-release \
  --target mglet -j 4
```

Do not compare Debug and Release results.

## Workload tiers

`benchmarks.json` is the self-contained workload definition.

- `quick`: 10,000 particles and 10 timesteps for both diffusion-focused and
  advection-diffusion cases. Intended for rapid checks, not final performance
  claims.
- `medium`: 100,000 particles × 125 steps for each single-rank compute case,
  plus a four-rank 50,000-particle × 500-step exchange case.
- `overnight`: the original 8,000,000-particle diffusion workload and
  1,000,000-particle advection-diffusion workload. These are not routine
  baselines.

Each benchmark performs one discarded warmup and three measured repetitions.
Before those runs, an untimed zero-step preflight confirms the initialized
global particle count. Reported values are medians. Snapshots, statistics,
frequent `itinfo`, terminal chatter, and particle checkpoint output are disabled
in staged benchmark parameters.

The CPU implementation honors `do_advection` and `do_diffusion`, matching the
OpenMP target path. Consequently, `diffusion-focused-*` measures the
Diffusion-Cube workload without running the disabled advection kernel.
The profiles provide deterministic per-rank GCC 13 `particle_seed` arrays so
repetitions use the same initialization and diffusion realization. On hyd38,
calibration pilots took about 6 minutes per measured compute case and 5 minutes
for the exchange case. The routine self-contained medium tier is therefore
expected to take about 75 minutes including warmups.

`bcc-benchmarks.json` defines optional production-faithful cases that bind
read-only inputs from `/home/yaydin/thesis/particle_run/10_BCC`. The harness
permits no other external root, hashes every bound input, and writes only below
`/home/yaydin/particle-performance/benchmark-results`.

- `bcc-compute-medium`: 200,000 particles × 1,000 steps on 16 ranks with
  frozen flow, ghost-cell IB, obstacles, advection, and Gaussian diffusion;
  snapshots and final particle output are disabled.
- `bcc-io-medium`: the same compute window with production snapshot cadence
  (`snapshot_step = 195`) and timer 960 required.
- `bcc-production-overnight`: the exact 39,075-step production window with
  final particle output, intended for occasional validation only.

Calibration took about 11 minutes for one BCC compute repetition. The BCC I/O
pilot was heavily affected by concurrent jobs and is only an upper bound; do
not treat calibration timings as a baseline.

## Running

Run the quick tier with fixed rank placement:

```bash
python3 tests/Particles/Performance/particle_bench.py run \
  --binary /home/yaydin/particle-performance/build-particle-release/src/mglet \
  --tier quick \
  --mpirun-arg=--bind-to \
  --mpirun-arg=core
```

Run the self-contained medium tier:

```bash
python3 tests/Particles/Performance/particle_bench.py run \
  --binary /home/yaydin/particle-performance/build-particle-release/src/mglet \
  --tier medium \
  --mpirun-arg=--bind-to \
  --mpirun-arg=core
```

Run production BCC compute and I/O separately:

```bash
python3 tests/Particles/Performance/particle_bench.py run \
  --binary /home/yaydin/particle-performance/build-particle-release/src/mglet \
  --manifest tests/Particles/Performance/bcc-benchmarks.json \
  --tier medium \
  --case bcc-compute-medium \
  --mpirun-arg=--bind-to \
  --mpirun-arg=core

python3 tests/Particles/Performance/particle_bench.py run \
  --binary /home/yaydin/particle-performance/build-particle-release/src/mglet \
  --manifest tests/Particles/Performance/bcc-benchmarks.json \
  --tier medium \
  --case bcc-io-medium \
  --mpirun-arg=--bind-to \
  --mpirun-arg=core
```

Add `--pilot` during calibration. Pilot mode removes warmups, runs one measured
repetition, and marks the result ineligible for baseline recording.

Select one benchmark by repeating `--case`, for example:

```bash
python3 tests/Particles/Performance/particle_bench.py run \
  --binary /home/yaydin/particle-performance/build-particle-release/src/mglet \
  --tier quick \
  --case diffusion-focused-quick
```

The results path is printed when the run completes. A session contains:

```text
benchmark-results/<timestamp>/
  environment.json
  summary.json
  summary.txt
  correctness/
  benchmarks/
```

Failed sessions retain their staged inputs, logs, and partial `summary.json`.

## Correctness gates

The default gates are:

1. deterministic one-rank advection with particle count, unique ID, and finite
   coordinate checks;
2. deterministic two-rank advection with the same checks;
3. uninterrupted versus split HDF5 restart continuation, comparing particles
   after sorting by `ipart`.

The deterministic smoke cases also check their expected final coordinates. The
restart check requires identical IDs, counts, schema, state/grid fields, any RNG
field present in the file, and a maximum coordinate difference of `1e-6`. A
failure stops the session before timed workloads. The current CPU gate has zero
diffusion and therefore does not prove stochastic RNG continuation; add that
case when the OpenMP particle build is usable on the selected compiler.

## Recording and comparing a baseline

Record a successful summary as a local baseline:

```bash
python3 tests/Particles/Performance/particle_bench.py record-baseline \
  --result /home/yaydin/particle-performance/benchmark-results/<run>/summary.json \
  --output /home/yaydin/particle-performance/benchmark-results/baseline.json
```

Compare a later run:

```bash
python3 tests/Particles/Performance/particle_bench.py compare \
  --result /home/yaydin/particle-performance/benchmark-results/<new-run>/summary.json \
  --baseline /home/yaydin/particle-performance/benchmark-results/baseline.json
```

The comparison reports percentage changes for wall time and every common timer.
It refuses cases whose effective parameters, source case, rank count, or
configured particle count differ, and the workload fingerprint includes hashes
of grid and dictionary inputs. It also refuses compiler, host, MPI, module,
placement, or affinity mismatches. For exploratory analysis only, the latter
check can be bypassed with `--allow-environment-mismatch`. The comparison does
not enforce a timing threshold.

## Recorded metrics

The harness records wall time, generated output bytes, timer totals and
instances, configured particle count, ranks, parameters, hostname, CPU,
compiler/CMake configuration, loaded modules, placement arguments, Git commit,
branch, and dirty state.

Timer data comes from both sections of `mglet-perf-report.txt`:

- 1 and 2: whole program and timeloop;
- 900: particle simulation;
- 920–925: integration and its subregions;
- 940–945: exchange and its subregions;
- 960 and 990: snapshots and particle finalization.

Use inclusive timer 920 for total particle integration and inclusive timer 940
for total exchange. CPU and OpenMP builds assign different meanings to timers
921–925, so those sub-timers must not be compared across build modes.

## Reproducibility requirements

Only compare runs with the same:

- machine and CPU allocation;
- compiler, optimization flags, precision, and OpenMP mode;
- MPI rank count, rank placement, and thread affinity;
- workload manifest and input grid;
- background-load policy.

Use multiple repetitions and medians. Keep the full raw sessions; do not rely
only on the text summary. Before authoritative runs, verify that no other
MGLET or MPI workloads are using the selected cores.

## Harness unit tests

```bash
PYTHONDONTWRITEBYTECODE=1 \
  python3 tests/Particles/Performance/test_particle_bench.py
```
