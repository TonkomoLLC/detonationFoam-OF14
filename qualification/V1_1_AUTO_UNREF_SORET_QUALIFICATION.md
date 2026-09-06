# OpenFOAM 14 v1.1.0 automatic-unrefinement and published-OF8 Soret qualification

Status: **PASSED on OpenFOAM Foundation 14** (2026-09-06).

## Automatic planar unrefinement

Runtime qualification passed for:

- uniform slab cycle: `1200 -> 4800 -> 1200 -> 4800 -> 1200`;
- localized slab cycle: `1200 -> 2352 -> 1200 -> 2352 -> 1200`;
- annular wedge cycle: `1600 -> 6400 -> 1600 -> 6400 -> 1600`;
- volume/uniform/tracer integral relative conservation: `<= 2e-10`;
- serial `checkMesh` at refined and coarsened topology states;
- MPI2 fixed-decomposition reversible cycle;
- intentional refusal of restart with active reversible ancestry;
- successful restart from a fully coarsened state followed by re-refinement.

Runtime mesh redistribution/load balancing while reversible ancestry is active
remains intentionally refused by design.

## Published OF8 H/H2 Soret compatibility

`legacyThermalDiffusionMode publishedOF8` intentionally reproduces the source
behavior published in the supplied OF8 code: both H and H2 thermal-diffusion
polynomial sums accumulate into the H ratio, while the direct H2 ratio remains
zero. This is a reproducibility mode, not a claim that the legacy assignment is
a corrected modern Soret formulation.

The direct OpenFOAM 14 species-flux gate passed with a synthetic diagnostic
coefficient set:

- `publishedOF8` loaded the original `constant/thermoDiff` layout;
- H2-only polynomial drove H flux: `max|jH| = 1.690421e-03`;
- inferred direct H2 Soret contribution: relative residual `1.595e-16`;
- published H/H2 corrected-flux identity: relative residual `1.603e-16`;
- zero-net-mass species-flux closure: relative residual `1.283e-16`;
- Soret-off control: `max|jH| = 0`, `max|jH2| = 0`.

The synthetic coefficient is a qualification signal only and is not physical
transport data.

## Integrated reacting laptop smoke

A 120-cell H2/O2 `foamRun` / `detonationFluid` case exercised the two new
features together and passed:

- `publishedOF8` Soret mode constructed;
- `planarRefiner` automatic unrefinement constructed;
- repeated refinement events occurred;
- repeated automatic unrefinement events occurred;
- no OpenFOAM fatal error or floating-point exception;
- normal solver termination;
- final temperature range: `300 .. 3581.900841 K`;
- final leading-shock position: `0.00051 m`.

This is a compact integration smoke, not formal detonation validation.

## Reproducibility commands

The detailed isolated feature gates are retained under the historical
Candidate-1 qualification directory because those exact scripts produced the
accepted runtime evidence:

```bash
./qualification/candidate1Runtime/RunPlanarAutoUnrefRuntimeGate
./qualification/candidate1Runtime/RunPublishedOF8SoretRuntimeGate
```

The integrated laptop tutorial is:

```bash
cd tutorials/H2_O2_laptop_autoUnref_Soret_OF14
./Allrun
```
