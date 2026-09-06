# Candidate 1 Runtime Gate 1 Status

**Status:** prepared; requires execution in the user's OpenFOAM Foundation 14 environment.

The production Candidate 1 source is unchanged from the compile-passed tree.
This gate adds only qualification drivers, cases, analyzers, and runners.

## Automatic unrefinement acceptance targets

| Gate | Target |
|---|---|
| AU-R1 | Uniform slab `1200 -> 4800 -> 1200 -> 4800 -> 1200` |
| AU-R2 | Localized slab `1200 -> 2352 -> 1200 -> 2352 -> 1200` |
| AU-R3 | Volume, uniform-scalar integral, and linear-tracer integral conserved through every topology change; every serial state passes `checkMesh` |
| AU-R4 | Fixed-decomposition MPI2 uniform slab completes the reversible cycle |
| AU-R5 | Restart with active reversible ancestry is refused deliberately |
| AU-R6 | Restart from fully coarsened state succeeds and can re-refine |

Runtime redistribution/load balancing remains intentionally outside Candidate 1 reversible mode because Candidate 1 rejects `distribute()` while reversible cut ancestry is active.

## Published OF8 H/H2 Soret acceptance targets

The diagnostic uses uniform composition

```text
Y_H  = 0.005
Y_H2 = 0.100
Y_N2 = 0.895
```

with an axial temperature gradient.  Every H thermal-diffusion polynomial is zero; only

```text
H2-N2/ThermDiff_1 = 100
```

is nonzero.  This coefficient is synthetic and qualification-only.

For the exact published OF8 assignment, the H2 polynomial contributes to the H raw thermal flux and the direct H2 raw thermal flux remains zero.  With uniform composition and the legacy zero-net-mass correction,

```text
j_H  = -(1-Y_H) q
j_H2 = Y_H2 q
```

and therefore

```text
(1-Y_H) j_H2 + Y_H2 j_H = 0.
```

The gate requires:

- the real `legacyMixtureAverageFourier` runtime model to load `constant/thermoDiff` in `publishedOF8` mode;
- nonzero H flux produced even though all H polynomial coefficients are zero;
- the H/H2 corrected-flux identity above to hold face-by-face;
- an independently inferred direct H2 Soret contribution to be zero;
- the sum of all 33 corrected species diffusion fluxes to close to zero;
- the otherwise identical `legacyThermalDiffusionMode off` control to give zero uniform-composition diffusive flux.

## Run

```bash
./qualification/candidate1Runtime/RunCandidate1RuntimeQualification \
    2>&1 | tee log.Candidate1.runtimeQualification
```

If a gate fails, return `log.Candidate1.runtimeQualification` and the generated `qualification/candidate1Runtime/results*` directories.
