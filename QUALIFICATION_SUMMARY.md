# OpenFOAM 14 detonationFoam v1.1.0 qualification summary

**Final release status: PASSED / RELEASED.**

This file is the authoritative release-level qualification summary. Version
1.1.0 retains the accepted v1.0.0 A-G/R1 migration evidence and adds qualified
automatic planar unrefinement plus published-OF8 H/H2 Soret compatibility.

## Closed qualification areas

| Area | Final release disposition |
|---|---|
| A - OF14 modular build/startup and OF8 cross-version baseline | Passed |
| B - shock flux families | Passed |
| C - serial/MPI matrix, restart/redecomposition, OF8<->OF14 profile equivalence | Passed |
| D - legacy transport/property compatibility and binary diffusion | Passed |
| E - native OF14 3-D AMR, localized AMR, MPI/load balancing, refined restart | Passed |
| F - reusable planar/slab and wedge/axisymmetric `planarRefiner` | Passed |
| G1-G4 - integrated laptop regression, equivalence, sensitivity, performance/load-balance smoke | Passed |
| R1 - v1.0.0 clean package/build/install/runtime smoke | Passed |
| v1.1-U - reversible automatic planar unrefinement | **Passed** |
| v1.1-S - published OF8 H/H2 Soret compatibility | **Passed** |
| v1.1-I - coupled H2/O2 laptop integration smoke | **Passed** |

## v1.1-U automatic-unrefinement evidence

- Uniform slab: `1200 -> 4800 -> 1200 -> 4800 -> 1200`.
- Localized slab: `1200 -> 2352 -> 1200 -> 2352 -> 1200`.
- Annular wedge: `1600 -> 6400 -> 1600 -> 6400 -> 1600`.
- Volume/uniform/tracer integral relative conservation: `<= 2e-10`.
- Serial `checkMesh` passed at topology states.
- MPI2 fixed-decomposition reversible cycle passed.
- Restart with active reversible ancestry was intentionally refused.
- Restart from a fully coarsened state followed by re-refinement passed.

Runtime redistribution/load balancing with active reversible ancestry remains
intentionally unsupported in v1.1.0.

## v1.1-S published-OF8 Soret evidence

The runtime flux gate used OpenFOAM 14's public multicomponent species-flux
interface and a synthetic diagnostic coefficient set. It passed:

- H2-only polynomial drove H flux: `max|jH| = 1.690421e-03`;
- inferred direct H2 Soret residual: `1.595e-16` relative;
- published H/H2 corrected-flux identity: `1.603e-16` relative;
- zero-net-mass species-flux closure: `1.283e-16` relative;
- Soret-off control: exactly zero H/H2 Soret flux.

The mode reproduces the published OF8 H2-to-H assignment intentionally. It is
not a corrected-Soret claim, and no physical coefficient dataset was invented.

## v1.1-I integrated H2/O2 laptop smoke

The 120-cell `foamRun`/`detonationFluid` tutorial exercised both new features in
a reacting run and passed all acceptance checks. The accepted run reported:

- repeated planar refinement;
- repeated automatic unrefinement;
- no OpenFOAM fatal error/FPE;
- normal solver termination;
- final `Tmin/Tmax = 300 / 3581.900841 K`;
- final leading-shock position `0.00051 m`.

This is an integration smoke, not formal CJ/grid/experimental validation.

## Selected retained v1.0.0 evidence

- F3 localized planar slab: 1200 -> 2352 cells with bounded transition topology.
- F3 wedge portability with stock `incompressibleFluid`: 1600 -> 6400; `Mesh OK.`.
- G2 OF8/OF14 replay: 33 species; `max |sumY-1| = 1.212e-12`.
- G3 compact mesh/time-step characterization completed.
- G4 native load-balancer/Scotch serial/MPI2/MPI4 smoke completed.
- R1 independent AMR clean build, full solver clean build, link audit, and packaged `foamRun` smoke passed.

## Laptop qualification policy

Pre-release gates were intentionally kept at or below roughly 30 minutes on the
development laptop. Longer full-resolution, endurance, formal convergence,
formal scaling, corrected-Soret development, and reversible-history
serialization/distribution are deferred rather than being release blockers
solely because of runtime. See `qualification/POST_RELEASE_DEFERRED_QUALIFICATION.md`.
