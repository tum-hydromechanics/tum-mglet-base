# Particle performance optimizations

This note records the optimizations kept on the CPU Release particle path
(`MGLET_OPENMP=OFF`), the shared-kernel / lincon experiment that was measured
and reverted, and the motion/exchange prep bundle. Numbers below are median
wall / exclusive timer totals from the harness on `hyd39` unless noted.

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

Motion/exchange bundle keep-set confirmation (guards + RNG + §4 vs baselines):

- Self-contained: `benchmark-results/20260730-165059-204638/`
- BCC compute: `benchmark-results/20260730-165338-397453/`
- Comparisons: `*-keepset-comparison.json`
- Driver log: `benchmark-results/keepset-compare-driver.log`

Timer IDs (Release/CPU names):

| ID  | Region              | Meaning                                      |
|-----|---------------------|----------------------------------------------|
| wall| `wall_seconds`      | Full case elapsed time                       |
| 921 | `ADV_VELOCITY`      | Velocity interpolation (lincon)              |
| 922 | `ADV_MOTION`        | Advective `move_particle`                    |
| 924 | `DIF_RN_GENERATION` | Random-walk displacement sampling            |
| 925 | `DIF_MOTION`        | Diffusive `move_particle`                    |
| 941 | `PREP_COMM`         | Exchange preparation (target grid / triage)  |

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

**Commit:** `66d101b7` — *Speed up particle diffusion RNG with shared polar sampling.*

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

## 4. Motion / exchange prep bundle (kept)

**Status:** Kept. Measured on `hyd39` via
`tests/Particles/Performance/run-keepset-compare.sh` (sessions listed above).

**Intent:** After guards + RNG, motion (`922` / `925`) and exchange prep
(`941`) remained major costs. This bundle removes duplicate cell-lookup
metadata work and a redundant same-grid exchange refresh, plus two small
timestep-invariant hoists. It does **not** revisit the reverted lincon
shared-stencil experiment, change RNG stream semantics, or alter MPI/HDF5
particle layout.

### Changes

| Piece | Where | What |
|-------|-------|------|
| Cached grid context for cell updates | `particle_basetype_mod.F90`, `particle_boundaries_mod.F90`, `particle_timeintegration_mod.F90` | Optional `kk/jj/ii` + `X/Y/Z` pointers into `update_particle_cell`; threaded through CPU `move_particle`. Removes unused `DX/DY/DZ` and `get_bbox` from the cell-update fallback. Falls back when `particle%igrid` no longer matches the bound grid (e.g. after replace). |
| Skip redundant exchange cell refresh | `particle_exchange_mod.F90` (CPU `#else` path) | When `destgrid` is unchanged **and** `iface == 0`, skip `update_particle_cell` (motion already refreshed `ijkcell`). Still refresh when `iface /= 0` (periodic same-grid wrap) and on all cross-grid paths. |
| Precomputed diffusion scales | `particle_diffusion_mod.F90`, `particle_timeintegration_mod.F90` | Compute `sqrt(2*D*dt)` once per timestep; pass scales into host displacement generation. Walk mode and RNG call order unchanged. |
| Reuse init-time RK coeffs | `particle_timeintegration_mod.F90` | CPU advection uses `A_offload` / `B_offload` instead of `prkscheme%get_coeffs` every particle/stage. Timer 921/922 placement inside the RK loop is unchanged. |

Also fixed a latent REAL64 kind typo in `particle_utils_mod.F90`
(`INTEGER(realk)` → `INTEGER(intk)` in `conditional_update_ri`) so double-precision
builds compile; single-precision Release was unaffected.

### Cumulative results vs authoritative baselines

Environment match: hyd39, GCC 13, Release, `MGLET_OPENMP=OFF`,
`--bind-to core`. Correctness gates passed. These wall numbers include
guards + RNG + this bundle.

| Case | Wall | Notable exclusive timers |
|------|------|--------------------------|
| `diffusion-focused-medium` | **−79.7%** (30.5 → 6.2 s) | 924 −85%; 925 −63%; 941 −92%; 921/922 ~0 (guards) |
| `advection-diffusion-medium` | **−58.3%** (47.8 → 19.9 s) | 924 −85%; 922 −74%; 925 −71%; 941 −94%; 921 −15% |
| `advection-diffusion-exchange-medium` | **−56.2%** (25.5 → 11.2 s) | 924 −78%; 922 −71%; 925 −72%; 941 −92%; 942 −57% |
| `bcc-compute-medium` | **−52.6%** (52.4 → 24.8 s) | 924 −78%; 922 −64%; 925 −59%; 941 −92%; 942 −54%; 921 −5% |

### Bundle attribution (approx.)

Against the earlier post-RNG / shared-kernel BCC session (~41 s wall, before
this bundle; still had the later-reverted lincon split on 921), the keep-set
BCC wall is ~25 s. The largest new exclusive cuts aligned with this section
are **941** (~2.9 → 0.23 s) and motion **922** / **925**. Timer **924** was
already dominated by the RNG keep-set.

Comparisons against the authoritative baselines alone are **cumulative**; an
isolated §4-only delta needs a clean post-RNG / pre-bundle reference under the
same environment.

---

## Keep-set

| Piece | Keep? |
|-------|-------|
| Physics guards (`b2689a2f`) | Yes |
| RNG walk-mode + polar truncated sampler + unified seed (`66d101b7`) | Yes |
| Lincon shared-stencil split | No (reverted) |
| Motion / exchange prep bundle (§4) | Yes |
| REAL64 `conditional_update_ri` kind fix | Yes |

Re-measure with
`tests/Particles/Performance/run-keepset-compare.sh` (expects Release `mglet` at
`/home/yaydin/particle-performance/build-particle-release/src/mglet`).
