# Gate B3 status — HLLCP restoration

## Status: CLOSED / PASSED

The returned OpenFOAM Foundation v14 result bundle confirms Gate B3 is closed.

Observed qualification:

- `Allwmake` rebuilt `libdetonationFluidSolver.so` successfully;
- `foamRun` selected `detonationFluid fluxScheme = HLLCP`;
- no fatal, floating-point, segmentation, or NaN runtime failure occurred;
- the reduced-domain leading shock reached **0.0102675 m**, exceeding the **0.0102625 m** two-cell propagation target;
- final reported values were approximately `pmax = 2.2385621 MPa` and `Tmax = 4992.31 K`;
- the run terminated normally with `End`.

The B3 offline algebra audit remains part of every later gate: 5,000 randomized HLLCP conservative-flux comparisons, including the pressure-ratio sensor, Mach/star-pressure blending, `phip`, arbitrary tangential velocity, face orientation/area, and ALE pressure work.
