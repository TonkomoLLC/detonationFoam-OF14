# Gate B1 status — HLL restoration

## Status: CLOSED / PASSED

The returned OpenFOAM Foundation 14 run compiled and linked `libdetonationFluidSolver.so`, selected `fluxScheme HLL`, advanced the leading shock beyond the 0.0102625 m fast-gate target to 0.0102675 m, and terminated normally with `End`.

The candidate-1 analyzer's `FAILED (runtime error token)` message was a false positive. It lower-cased the log and interpreted OpenFOAM's normal startup line `SIGFPE : Enabling floating point exception trapping (FOAM_SIGFPE)` as an actual floating-point exception. `RuntimeLogHealth.py` fixes this by recognizing only true fatal/signal forms, and its self-test is now part of the static audits.

Gate B1 is therefore closed; no rerun is required solely for that analyzer bug.
