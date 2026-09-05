# Gate D3-B status

**CANDIDATE 2 — runtime qualification pending.**

Scope:

- register native OpenFOAM-14 `logPolynomialTransport<8>` with `coefficientWilkeMulticomponentMixture`;
- retain native Wilke mixture viscosity because it is algebraically identical to the OF8 rule;
- implement the missing OF8 arithmetic/harmonic mixture conductivity inside `legacyMixtureAverageFourier`;
- retain D2-B exact binary diffusion;
- keep Soret off;
- do not modify `detonationFluid`.

Offline static/algebra qualification: **PASSED**.

Candidate 2 adds the OpenFOAM-14 `thermo.H` prerequisite before `forThermo.H` in the custom thermo registration translation unit. No transport algebra or solver source is changed.
