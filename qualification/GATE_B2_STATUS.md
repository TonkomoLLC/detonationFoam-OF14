# Gate B2 status — HLLC restoration

## Status: CLOSED / PASSED

The returned OpenFOAM Foundation v14 result bundle confirms that Gate B2 is closed.

Observed qualification:

- `Allwmake` rebuilt `libdetonationFluidSolver.so` successfully;
- `foamRun` selected `detonationFluid fluxScheme = HLLC`;
- no fatal, floating-point, segmentation, or NaN runtime failure occurred;
- the reduced-domain leading shock reached **0.0102675 m**, exceeding the **0.0102625 m** two-cell propagation target;
- final reported values were approximately `pmax = 2.2383867 MPa` and `Tmax = 4992.31 K`;
- the run terminated normally with `End`.

The B2 offline algebra audit also remains part of every later gate: 5,000 randomized HLLC conservative-flux comparisons, including tangential velocity, arbitrary face orientation/area, and ALE pressure work.
