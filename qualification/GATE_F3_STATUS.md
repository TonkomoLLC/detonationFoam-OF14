# Gate F3 status - PASSED

Stage F3 passed the reusable planar/wedge AMR qualification on OpenFOAM Foundation 14.

- localized true-2-D slab AMR passed: 1200 -> 2352 cells with 384 selected cells;
- bounded transition topology: four 7-face transition polyhedra / four concave transition cells;
- wedge/axisymmetric portability passed with stock `incompressibleFluid`: 1600 -> 6400 cells and final `Mesh OK.`;
- MPI2 with native `loadBalancer` passed;
- refined/distributed in-place restart restored refinement state without duplicate refinement;
- continuous-vs-restart geometry and field comparisons passed;
- no `detonationFluid` production-source patch was required.

The reusable library is `libplanarFvMeshTopoChangers.so`, runtime type `planarRefiner`.
