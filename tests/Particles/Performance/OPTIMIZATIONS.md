# Particle performance optimizations

This note records the two optimizations currently kept on the CPU Release
particle path (`MGLET_OPENMP=OFF`), and the one shared-kernel experiment that
was measured and reverted. Numbers below are median wall / exclusive timer
totals from the harness on `hyd39`.

Authoritative baselines (pre-optimization gold):

- `benchmark-results/self-contained-medium-baseline.json`
- `benchmark-results/bcc-compute-medium-baseline.json`

Post-guards clean reference sessions (guards only, before RNG keep-set):

- Self-contained: `benchmark-results/20260719-145914-762784/`
- BCC compute: `benchmark-results/20260719-150533-923471/`

Shared-kernel measurement session (included a later-reverted lincon split):

- Self-contained: `benchmark-results/20260719-154034-886572/`
- BCC compute: `benchmark-results/20260719-154522-760502/`
- Comparisons: `*-shared-kernels-comparison.json`

Timer IDs (Release/CPU names):

| ID  | Region              | Meaning                                      |
|-----|---------------------|----------------------------------------------|
| wall| `wall_seconds`      | Full case elapsed time                       |
| 921 | `ADV_VELOCITY`      | Velocity interpolation (lincon)              |
| 922 | `ADV_MOTION`        | Advective `move_particle`                    |
| 924 | `DIF_RN_GENERATION` | Random-walk displacement sampling            |
| 925 | `DIF_MOTION`        | Diffusive `move_particle`                    |

---

## 1. Physics guards (`do_advection` / `do_diffusion`)

**Commit:** `b2689a2f` — *Honor do_advection/do_diffusion on the CPU particle path.*

**Change:** On the CPU time-integration path, skip disabled advective or
diffusive kernels instead of always running both. Disabled timer regions are
still touched so the harness timer set stays stable.

**Where it helps:** Cases that only need diffusion (notably
`diffusion-focused-medium`) no longer pay for unused lincon + advective motion.

### Results vs authoritative baselines (post-guards clean remeasure)

| Case | Wall | Notes |
|------|------|-------|
| `diffusion-focused-medium` | **−35.5%** | 921/922 drop to ~0 (advection skipped) |
| `advection-diffusion-medium` | ~flat (±1%) | Both physics modes still on |
| `advection-diffusion-exchange-medium` | ~flat | |
| `bcc-compute-medium` | ~flat | |

Guards are correctness-preserving configuration honoring, not a change to the
math of enabled kernels.

---

## 2. Shared / faster diffusion RNG (kept)

**Status:** Kept in the working tree for commit (with lincon split reverted).

**Change (summary):**

- Integer walk-mode enum (`rw_uniform` / `rw_gaussian2` / `rw_rademacher`) set
  once at init; hot path branches on the integer, not strings.
- Pure `declare target` transforms (`sample_uniform`, `sample_rademacher`,
  Marsaglia polar pair + truncated scaling) shared by host and OpenMP device
  adapters.
- Host adapter still uses `RANDOM_NUMBER` (seeded host behavior preserved).
- Device adapter uses per-particle LCG and the same pure transforms.
- Pair-producing truncated Gaussian reuses the spare polar variate across
  x/y/z draws for one particle.
- `baseparticle_t%seed` always present (unified host/device/MPI/HDF5 layout).

**Unit test:** `tests/Particles/Performance/test_gaussian_rng.py`

### Results (RNG effect)

Against the **authoritative baselines**, the shared-kernel compare run reported
large wall wins that mix guards (diffusion-focused) with RNG:

| Case | Wall vs baseline | Exclusive 924 vs baseline |
|------|------------------|---------------------------|
| `diffusion-focused-medium` | −61.3% | −85.0% |
| `advection-diffusion-medium` | −18.3% | −85.1% |
| `advection-diffusion-exchange-medium` | −23.7% | −77.9% |
| `bcc-compute-medium` | −21.0% | −77.4% |

Against **post-guards clean** sessions (RNG attribution only; same compare
session still had the later-reverted lincon split):

| Case | Wall | Exclusive 924 |
|------|------|---------------|
| `diffusion-focused-medium` | −40.0% | −85.1% |
| `advection-diffusion-medium` | −18.8% | −85.8% |
| `advection-diffusion-exchange-medium` | −23.6% | −77.8% |
| `bcc-compute-medium` | −20.9% | −77.4% |

**Takeaway:** Timer **924** is the dominant win. Prefer BCC `924` and
advection–diffusion wall when claiming the RNG upgrade.

---

## 3. Lincon gather / shared stencil (reverted)

**Status:** Reverted. Host/device interpolation restored to the pre-split
`interpolate_lincon` / `interpolate_lincon_target` bodies.

**What was tried:** Gather into a local stencil, then call a shared
`lincon_from_local_stencil` Gobert kernel from both host and device paths.

**Why it was dropped:** Exclusive timer **921** regressed on important cases
in the same compare session (~+16% on `bcc-compute-medium` and
`advection-diffusion-exchange-medium` vs post-guards). That failed the
~1% idle-host keep bar for interpolation claims.

Motion helpers / CPU in-bbox fast path from the same plan pass were left in
place; they showed only small mixed movement on 922/925 and were not the
primary win or the primary regression.

---

## Keep-set for commit

| Piece | Keep? |
|-------|-------|
| Physics guards (`b2689a2f`) | Yes (already committed) |
| RNG walk-mode + polar truncated sampler + unified seed | Yes |
| Lincon shared-stencil split | No (reverted) |

Optional confirmation after this keep-set lands: from the repo root of
`tum-mglet-base`, run
`tests/Particles/Performance/run-keepset-compare.sh` (self-contained medium +
BCC compute vs the authoritative baselines). That is good practice before
calling the RNG numbers final on a lincon-reverted binary, but not required if
the tree matches this document and correctness gates pass. The wrapper expects
Release `mglet` at
`/home/yaydin/particle-performance/build-particle-release/src/mglet`.
