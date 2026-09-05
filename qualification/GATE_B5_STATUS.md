# Gate B5 status — AUSM+up

**CLOSED / PASSED** on OpenFOAM Foundation 14.

Returned runtime evidence:

- `libdetonationFluidSolver.so` compiled and linked successfully;
- `foamRun` selected `detonationFluid` with `fluxScheme = AUSM+up`;
- final reported time: `3.78464899e-08 s`;
- final leading-shock location: `0.0102675 m`;
- propagation target: `0.0102625 m`;
- final `pmax = 2241186.53 Pa`;
- final `Tmax = 4992.26559 K`;
- normal OpenFOAM `End` marker present.

Gate B5 therefore closes the original custom flux-family restoration. The modular OF14 solver now contains Kurganov, Tadmor, HLL, HLLC, HLLCP, AUSM+, and AUSM+up.
