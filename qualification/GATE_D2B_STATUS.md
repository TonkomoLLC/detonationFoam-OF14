# Stage D2-B status

**Candidate 1 — READY FOR USER COMPILE/RUNTIME QUALIFICATION**

Scope:

- carry forward the passed D2 `legacyMixtureAverageFourier` transport model;
- add exact OF8 `Diff1..Diff4` pressure/temperature binary-diffusion evaluation as `legacyBinaryDiffusionCoefficient`;
- keep `legacyThermalDiffusionMode off` for runtime qualification;
- record an explicit Soret compatibility decision;
- leave `detonationFluid` and all A/B/C-qualified flux/chemistry/load-balancing code unchanged.

Pass criteria:

1. D2 returned PASS evidence is present.
2. A/B/C solver source hashes remain unchanged.
3. D2-B algebra/formula/Soret audit passes.
4. New Function2 library compiles in OpenFOAM Foundation 14.
5. `foamRun` selects `legacyMixtureAverageFourier` and advances with 528 `legacyBinaryDiffusionCoefficient` pair entries to normal `End`.
