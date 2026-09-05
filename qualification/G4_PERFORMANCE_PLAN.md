# Stage G4 — compact serial/MPI performance and native-load-balancing smoke

## Purpose

G4 closes Stage G with a laptop-scale parallel execution check. It is intentionally a **microbenchmark/smoke test**, not a formal scalability study.

## Qualified matrix

All branches use the same 2400-cell, 12-mm initialized detonation segment, HLLCP fluxes, the D3-B exact legacy-property layer, `maxCo=0.1`, and `endTime=2.7e-8 s`.

- serial: 1 process, no distributor;
- MPI2: deterministic `simple (2 1 1)` initial decomposition for `decomposePar`; before `foamRun`, `system/decomposeParDict` is replaced by an equivalent 2-domain `method scotch` dictionary, so the native OF14 `loadBalancer` is constructed with Scotch;
- MPI4: the same two-phase handoff with `simple (4 1 1)` for the initial split and `method scotch` for runtime redistribution.

MPI cases must show `Selecting fvMeshDistributor loadBalancer`, `Selecting distributor scotch`, native load-imbalance diagnostics, and at least one actual `Redistributing mesh` event. Candidate 1 was intentionally not accepted after its returned log showed `Selecting distributor simple`; candidate 2 corrects that configuration rather than weakening the audit. Final shock position, peak pressure, and peak temperature are compared with the serial branch using loose smoke-test consistency limits.

## Timing interpretation

Wall-clock time, speedup, and parallel efficiency are recorded. **No positive-speedup requirement is imposed.** At 2400 cells, MPI startup, communication, chemistry-load instrumentation, and redistribution overhead can be comparable to or larger than useful parallel work. Requiring speedup here would be a misleading release criterion.

Formal strong/weak scaling, larger meshes, longer trajectories, and load-balancer amortization studies are post-release work.

## Runtime policy

Default hard budget: **900 s (15 minutes)**. This remains below the project-wide ~30-minute laptop ceiling.
