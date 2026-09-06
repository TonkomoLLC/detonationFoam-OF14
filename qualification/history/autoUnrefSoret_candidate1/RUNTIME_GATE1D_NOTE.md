# Candidate 1 Runtime Gate 1d

Gate 1c showed that the published-OF8 Soret transport model itself constructs,
reads `constant/thermoDiff`, and enables the intended compatibility mode.  The
qualification executable then stopped because both tiny transport-probe cases
were missing `system/fvSolution`.

Gate 1d adds a minimal valid OpenFOAM 14 `fvSolution` dictionary to both
`soretPublishedOF8` and `soretOff`.  The probe performs no PDE linear solve, so
an empty `solvers {}` dictionary is sufficient; the file exists solely because
OF14 model construction expects the solution dictionary to be present.

No production source, transport implementation, AMR implementation, physics
coefficient, or runtime-probe C++ source is changed by this overlay.
