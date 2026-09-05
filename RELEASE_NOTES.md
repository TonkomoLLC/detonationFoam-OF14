# detonationFoam for OpenFOAM Foundation 14 - v1.0.0 release notes

## Release status

**Released / qualified for OpenFOAM Foundation 14.**

The staged migration program is closed through Stage G4 and the R1 clean-build/install/runtime qualification is accepted. The R1 runner's reported failure was a known false-positive health check: it matched OpenFOAM's normal startup message `sigFpe : Enabling floating point exception trapping (FOAM_SIGFPE).` The smoke command itself returned successfully after both clean build paths passed. The final release checker now distinguishes that startup message from a real floating-point exception.

## Major delivered capability

- Modular OpenFOAM 14 `foamRun` solver module: `detonationFluid`.
- Restored flux families including Kurganov/Tadmor, HLL, HLLC, HLLCP, AUSM+, and AUSM+up paths qualified during the staged migration.
- Native OpenFOAM 14 thermo/chemistry integration.
- Exact legacy mixture-average/property compatibility layer through `libdetonationLegacyThermophysicalTransportModels.so` where stock OF14 behavior is not equivalent.
- Native OpenFOAM 14 3-D AMR/load balancing/restart integration.
- Reusable true-2-D slab and wedge/axisymmetric AMR library `libplanarFvMeshTopoChangers.so`, runtime type `planarRefiner`.
- MPI execution with native `loadBalancer` and Scotch redistribution.
- Restart/redecomposition and refined/distributed restart qualification.
- OF8-to-OF14 migration documentation and packaged reference/fast tutorials.

## Reusable AMR packaging

The 2-D/wedge AMR component is intentionally solver-independent. It may be built without detonationFoam using either:

```bash
./AllwmakeAMR
```

or:

```bash
cd src/planarFvMeshTopoChangers
./Allwmake
```

The full detonation build remains separate:

```bash
./Allwmake
```

`libdetonationFluidSolver.so` has no runtime link dependency on `libplanarFvMeshTopoChangers.so`.

## Compact release evidence

- G1 integrated laptop regression: 205 s.
- G2 compact OF8<->OF14 equivalence replay: 31 s; 33 species; `max |sumY-1| = 1.212e-12`.
- G3 numerical/physical screen: 249 s; 1200/2400/4800-cell mesh screen plus `maxCo=0.10/0.05` time-step screen.
- G4 compact performance/load-balance smoke: serial 152.966 s; MPI2 91.855 s; MPI4 102.875 s; 29 native redistributions in both MPI branches; identical final shock position, peak pressure, and peak temperature across serial/MPI2/MPI4.
- F3 reusable planar/wedge AMR: localized slab, wedge/`incompressibleFluid` portability, MPI/load balancing, and refined/distributed restart passed.
- R1: independent AMR clean build passed; full clean build passed; no AMR link dependency; packaged fast smoke command returned successfully. R1 post-run failure was solely the corrected `FOAM_SIGFPE` string-matching false positive.

## Deferred post-release work

The release intentionally does not block on long laptop-unfriendly studies. Full-resolution long-horizon equivalence, formal asymptotic grid convergence, formal CJ/experimental validation, endurance/restart accumulation, and formal strong/weak scaling are listed in `qualification/POST_RELEASE_DEFERRED_QUALIFICATION.md`.
