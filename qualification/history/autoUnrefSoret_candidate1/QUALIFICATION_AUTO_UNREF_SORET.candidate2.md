# Automatic planar unrefinement + published OF8 H/H2 Soret qualification

Status: **PASSED on OpenFOAM 14** (user runtime gate, 2026-09-06).

This package consolidates Candidate 1 production code with the corrected Runtime Gate 1 qualification harness. The production implementation is unchanged from the compile-passed Candidate 1.

## Automatic planar unrefinement

Runtime qualification passed for:

- Uniform slab reversible cycle: `1200 -> 4800 -> 1200 -> 4800 -> 1200`.
- Localized slab reversible cycle: `1200 -> 2352 -> 1200 -> 2352 -> 1200`.
- Annular wedge reversible cycle: `1600 -> 6400 -> 1600 -> 6400 -> 1600`.
- Volume/uniform/tracer integral relative conservation tolerance: `<= 2e-10`.
- Serial `checkMesh` at topology states.
- MPI2 fixed-decomposition reversible cycle.
- Intentional refusal to restart with active reversible ancestry.
- Restart after complete coarsening followed by successful re-refinement.

Runtime redistribution/load balancing while reversible ancestry is active remains intentionally refused by design.

## Published OF8 H/H2 Soret compatibility

The compatibility mode intentionally reproduces the published OF8 source behavior in which the H2 thermal-diffusion polynomial contribution is accumulated into the H ratio and the direct H2 ratio remains zero. This is a reproducibility mode, not a claim that the legacy assignment is the corrected physical formulation.

Runtime qualification passed using the public OpenFOAM 14 species-flux interface with a synthetic diagnostic coefficient set:

- `publishedOF8` mode loaded `constant/thermoDiff`.
- H2-only diagnostic polynomial drove H flux: `max|jH| = 1.690421e-03`.
- Direct H2 Soret contribution inferred zero from corrected-flux identity: relative error `1.595e-16`.
- Published H/H2 corrected-flux identity: relative error `1.603e-16`.
- Zero-net-mass species-flux closure: relative error `1.283e-16`.
- Soret-off control: `max|jH| = 0`, `max|jH2| = 0`.

The synthetic coefficient used by the qualification case is only a diagnostic signal and is not a recommended physical transport coefficient.

## Qualification commands

```bash
./qualification/candidate1Runtime/RunPlanarAutoUnrefRuntimeGate
./qualification/candidate1Runtime/RunPublishedOF8SoretRuntimeGate
```

Both new production capabilities remain opt-in/default-off for backward compatibility.
