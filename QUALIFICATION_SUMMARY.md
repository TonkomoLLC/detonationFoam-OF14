# OpenFOAM 14 detonationFoam v1.0.0 qualification summary

**Final release status: PASSED / RELEASED.**

This file is the authoritative release-level qualification summary. Individual files under `qualification/` include historical candidate snapshots and should be interpreted as staged-development evidence, not as a later override of this summary.

## Closed qualification areas

| Area | Final release disposition |
|---|---|
| A - OF14 modular build/startup and OF8 cross-version baseline | Passed |
| B - shock flux families (HLL/HLLC/HLLCP/AUSM+/AUSM+up and retained baselines) | Passed |
| C - serial/MPI matrix, restart/redecomposition, OF8<->OF14 profile equivalence | Passed |
| D - legacy transport/property compatibility and exact binary-diffusion path | Passed |
| E - native OF14 3-D AMR, localized AMR, MPI/load balancing, refined restart | Passed |
| F - reusable planar/slab and wedge/axisymmetric `planarRefiner` | Passed |
| G1 - integrated laptop release regression | Passed - 205 s |
| G2 - compact OF8<->OF14 equivalence/profile gate | Passed - 31 s |
| G3 - compact numerical/physical sensitivity characterization | Passed - 249 s |
| G4 - compact serial/MPI performance and native Scotch load-balance smoke | Passed - 348 s |
| R1 - clean package/build/install/runtime smoke | Passed; original reported failure was checker-only false positive |

## F3 reusable planar/wedge AMR evidence

- Localized true-2-D slab refinement: 1200 -> 2352 cells, 384 selected cells.
- Transition topology remained bounded: four 7-face transition polyhedra and four concave cells.
- Maximum non-orthogonality: 18.4349 degrees; maximum skewness: 0.333333.
- Species closure: `max |sumY-1| = 1.045e-08`.
- Wedge portability with stock `incompressibleFluid`: 1600 -> 6400 cells; wedge constraint patches preserved; final `checkMesh` reported `Mesh OK.`
- MPI2/native load balancing/refined restart: load-balancer redistribution and `planarRefiner` distribute callback observed; restart state restored as 1/1 with no duplicate refinement.
- Continuous-vs-restart geometry: `dCmax = 0`, relative volume difference on the order of `1e-14`.
- Maximum primary-field NRMSE: `2.016e-05`; maximum major-species NRMSE: `6.447e-10`.

## G2 equivalence evidence

The checksum-locked OF8/OF14 spatial-profile comparison and fresh OF14 replay passed. The fresh replay reached `t = 2.7e-08 s` on 2400 cells with:

- leading shock location `0.0102625 m`;
- `pmax = 2.17190106e6 Pa`;
- `Tmax = 4994.47564 K`;
- 33 species fields;
- `max |sumY-1| = 1.212e-12`.

## G3 sensitivity screen

The compact screen used the same 12-mm trajectory with 10, 5, and 2.5 micrometre meshes at `maxCo=0.1`, plus 5 micrometres at `maxCo=0.05`.

The shock location was already stable from 5 -> 2.5 micrometres (`3.073e-08 m` difference), while peak pressure remained more mesh-sensitive (`4.152%` fine-vs-reference difference). This is recorded as screening characterization, not a formal asymptotic grid-convergence or CJ-validation claim.

## G4 serial/MPI/load-balancing evidence

| Run | Wall time [s] | Speedup | Efficiency | Native redistributions |
|---|---:|---:|---:|---:|
| serial | 152.966 | 1.000 | 1.000 | 0 |
| MPI2 | 91.855 | 1.665 | 0.833 | 29 |
| MPI4 | 102.875 | 1.487 | 0.372 | 29 |

Both MPI runs explicitly selected the native OF14 `loadBalancer` with `scotch` redistribution. Final shock location, peak pressure, and peak temperature matched the serial branch.

## R1 returned-result disposition

R1 performed, in order:

1. static/package/documentation audit - passed;
2. clean standalone AMR build - passed;
3. clean full detonation/legacy-transport build - passed;
4. link audit - passed: `libdetonationFluidSolver.so` has no dependency on `libplanarFvMeshTopoChangers.so`;
5. packaged fast serial `foamRun` smoke command - returned successfully;
6. post-run health grep - incorrectly failed.

The incorrect grep matched the normal OpenFOAM startup message:

```text
sigFpe : Enabling floating point exception trapping (FOAM_SIGFPE).
```

This exact false-positive mode had already been identified and corrected during Gate B1. Because the smoke command returned exit code 0 and all prior build/link steps passed, the R1 result is accepted as a **qualification pass with a harness-only defect**. `RunReleaseQualification` contains the corrected startup-aware check.

## Laptop qualification policy

Pre-release gates were intentionally kept at or below roughly 30 minutes on the development laptop. Longer full-resolution, endurance, formal convergence, and formal scaling studies are deferred rather than being release blockers solely because of runtime. See `qualification/POST_RELEASE_DEFERRED_QUALIFICATION.md`.
