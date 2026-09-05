# Stage D2 status

**Candidate 2 — READY FOR USER COMPILE/RUNTIME QUALIFICATION**

Stage D1 is closed/passed. D2 adds a runtime-selectable OpenFOAM 14 laminar thermophysical-transport model named `legacyMixtureAverageFourier` while leaving the A/B/C-qualified `detonationFluid` solver source unchanged.

Implemented in this candidate:

- legacy mixture-average binary closure
  `D_i = (1-Y_i)/sum_{j!=i}(X_j/D_ij)`;
- binary `D_ij(p,T)` inputs through native OF14 `Function2<scalar>` dictionaries;
- legacy molecular-weight-gradient flux correction;
- zero-net diffusive mass-flux correction across all species;
- corrected species flux exposed through OF14 `j()` / `divj()`;
- species sensible-enthalpy transport in `q()` / `divq()` using that same corrected flux;
- no diffusion-coefficient cache, so there is no new mesh-redistribution cache to invalidate.

Deliberately not enabled yet:

- legacy H/H2 Soret/thermal-diffusion path. The published OF8 source appears to write the H2 contribution into `TDRatio_H` rather than `TDRatio_H2`; candidate 2 accepts only `legacyThermalDiffusionMode off` so this ambiguity cannot be silently hidden.
- legacy custom viscosity/thermal-conductivity mixing. OF14 native thermo transport remains in ownership pending a D3 numerical equivalence test.
- turbulent legacy mixture-average extension. Native OF14 turbulence transport remains available separately.

The runtime smoke uses the same 33-species NH3/O2 case, generates all 528 binary pairs with a benign constant `D_ij=2e-5 m2/s`, and runs only long enough to prove runtime selection and equation coupling.

## Candidate-2 correction

Candidate 1 reached the OF14 compiler and failed only in the patch heat-flux
expression because unary minus was applied directly to `oneField`. Candidate 2
uses the native OF14 operator grouping and removes an unrelated unused-variable
warning. Solver physics and D2 transport formulas are unchanged.

## Candidate 3 correction

Candidate 2 compiled and advanced the smoke case to `End`, but OF14 selected default `unityLewisFourier`. The generated `constant/thermophysicalTransport` lacked a `FoamFile` header, so `header.headerOk()` failed and the selector deliberately took its default path. Candidate 3 adds a standard `class dictionary; object thermophysicalTransport;` header and explicitly rejects any default-model fallback in the runtime analyzer. No transport equations or `detonationFluid` source changed.
