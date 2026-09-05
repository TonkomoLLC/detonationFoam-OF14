# Stage G3 compact characterization plan

G3 is a **laptop-scale screening qualification**, not a formal detonation validation campaign.

## Cases

All cases use the same 12-mm canonical OpenFOAM 14 fast detonation segment and terminate at `2.7e-8 s`.

| Case | Cells | dx | maxCo | Purpose |
|---|---:|---:|---:|---|
| `dx10um_Co0p1` | 1200 | 10 µm | 0.10 | coarse spatial screen |
| `dx5um_Co0p1` | 2400 | 5 µm | 0.10 | reference |
| `dx2p5um_Co0p1` | 4800 | 2.5 µm | 0.10 | fine spatial screen |
| `dx5um_Co0p05` | 2400 | 5 µm | 0.05 | time-step sensitivity |

The 2400-cell initialized profiles are conservatively mapped onto the other meshes. Factor-two coarsening averages pairs of cells; factor-two refinement duplicates each parent value into two children. This preserves the initialized front location/integral content instead of stretching the profile.

## Metrics

The gate records leading threshold-front position, solver-reported shock position, peak pressure, peak temperature, an interpolated threshold-front displacement speed, all-33-species `sum(Y)` closure, and interpolated profile NRMSE for `p`, `T`, `Ux`, `NH3`, `O2`, `H2`, and `N2`.

The acceptance envelopes are intentionally screening-level. A failure indicates a numerical sensitivity large enough to investigate before release; a pass does **not** establish formal asymptotic grid convergence, a validated CJ velocity, or long-time detonation accuracy.

## Runtime policy

Default hard wall-clock limit is 900 s (15 min), below the project-wide ~30-minute laptop ceiling. Full-resolution convergence, extended propagation, CJ/experimental validation, broad parameter sweeps, and endurance runs remain post-release work.
