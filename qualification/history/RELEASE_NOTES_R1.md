# OpenFOAM 14 detonationFoam - R1 release-candidate notes

## Closed qualification before R1

Stages A-G are closed through G4. The final compact G4 run exercised serial, MPI2, and MPI4 HLLCP trajectories using the D3-B exact legacy-property layer and genuine native OpenFOAM 14 `loadBalancer` + Scotch redistribution.

G4 timing diagnostics on the 2400-cell laptop case were:

| Run | Wall time [s] | Speedup | Efficiency | Native redistributions |
|---|---:|---:|---:|---:|
| serial | 152.966 | 1.000 | 1.000 | 0 |
| MPI2 | 91.855 | 1.665 | 0.833 | 29 |
| MPI4 | 102.875 | 1.487 | 0.372 | 29 |

The timing data are diagnostic only; formal scaling is deferred.

## R1 packaging decisions

- `detonationFluid` remains a modular `foamRun` solver library.
- Exact legacy `NS_mixtureAverage` property/transport compatibility remains a separate runtime transport library rather than a duplicate legacy solver branch.
- The reusable true-2-D/wedge AMR capability is packaged as `src/planarFvMeshTopoChangers/` with its own `Allwmake` and root-level `AllwmakeAMR`.
- The full solver build does not require the optional AMR library.
- The AMR library does not link against `libdetonationFluidSolver.so`.
- Manuals are provided in Markdown, Word, and PDF, with a separate AMR manual as well as AMR coverage in the main solver manual.
- Long qualification tasks are not release-blocking under the laptop policy and are listed in `qualification/POST_RELEASE_DEFERRED_QUALIFICATION.md`.
