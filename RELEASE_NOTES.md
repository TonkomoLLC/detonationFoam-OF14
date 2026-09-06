# detonationFoam for OpenFOAM Foundation 14 - v1.1.0 release notes

## Release status

**Released / qualified for OpenFOAM Foundation 14.**

Version 1.1.0 is a feature release built on the qualified v1.0.0 OF14 port. The
new production capabilities are automatic reversible unrefinement in
`planarRefiner` and explicit reproduction of the published OF8 H/H2 Soret
behavior. Both remain opt-in/default-off.

## New in v1.1.0

### Automatic `planarRefiner` unrefinement

`planarRefiner` can now reverse its own planar directional cuts through
OpenFOAM's `undoableMeshCutter` history.

Key controls are:

```text
automaticUnrefinement          true;
unrefineInterval               1;
lowerUnrefineLevel             700;
upperUnrefineLevel             5000;
maxRefinementLevel             1;
maxUnrefinementPassesPerUpdate 2;
```

Coarsening requires both siblings of a reversible split to lie on the same
outside side of the hysteresis band. One planar level is two binary directional
cuts, so the default allows two undo passes per update. Per-cell binary cut
depth allows repeated refine/unrefine cycles without exhausting the old global
refinement-event counter.

Qualified runtime cycles were:

- slab uniform: `1200 -> 4800 -> 1200 -> 4800 -> 1200`;
- slab localized: `1200 -> 2352 -> 1200 -> 2352 -> 1200`;
- annular wedge: `1600 -> 6400 -> 1600 -> 6400 -> 1600`;
- MPI2 fixed decomposition: reversible slab cycle passed;
- integral conservation through topology changes: relative error `<= 2e-10`.

Safety restrictions are explicit: active reversible ancestry is not serialized
or redistributed, so active-history restart, mesh-to-mesh mapping, and runtime
mesh redistribution/load balancing are refused while reversible mode is active.
Restart from a fully coarsened state was qualified.

### Published OF8 H/H2 Soret behavior

`legacyMixtureAverageFourier` now accepts:

```text
legacyThermalDiffusionMode publishedOF8;
```

This mode intentionally reproduces the supplied OF8 source behavior in which
the H2 thermal-diffusion polynomial contribution is accumulated into the H
ratio and the direct H2 ratio remains zero. It is a compatibility/reproducibility
mode, **not** a corrected physical formulation.

The mode accepts either the original `constant/thermoDiff` layout or inline
`thermalDiffusionCoeffs`. Missing pair data are rejected rather than silently
reconstructed from the old `trandat`/`groupSpecies` shortcut.

The direct OF14 flux gate passed with:

- H2-only diagnostic polynomial driving H flux: `max|jH| = 1.690421e-03`;
- inferred direct H2 Soret residual: `1.595e-16` relative;
- published H/H2 corrected-flux identity: `1.603e-16` relative;
- zero-net-mass species-flux closure: `1.283e-16` relative;
- Soret-off control exactly zero for H/H2.

The qualification coefficients are synthetic diagnostics only; the supplied OF8
archive contained no physical `thermoDiff` dataset.

### Integrated H2/O2 laptop smoke

A new 120-cell tutorial exercises chemistry, published-OF8 Soret compatibility,
and repeated planar refine/unrefine in a single `foamRun` calculation:

```text
tutorials/H2_O2_laptop_autoUnref_Soret_OF14/
```

The accepted run reported repeated refinement and automatic-unrefinement events,
`Tmax = 3581.900841 K`, and a final leading-shock position of `0.00051 m` with
normal termination and no OpenFOAM fatal/FPE.

## Capability retained from v1.0.0

- modular `foamRun` solver module `detonationFluid`;
- Kurganov/Tadmor, HLL, HLLC, HLLCP, AUSM+, and AUSM+up flux paths;
- native OF14 thermo/chemistry integration;
- exact legacy mixture-average/property compatibility library;
- native OF14 3-D AMR/load-balancing/restart integration;
- reusable true-2-D slab and wedge/axisymmetric `planarRefiner` library;
- MPI/restart/redecomposition qualification;
- OF8-to-OF14 migration tooling/documentation.

## Backward compatibility

Both new features default to off:

```text
automaticUnrefinement false;
legacyThermalDiffusionMode off;
```

Existing v1.0.0 cases therefore do not enable either behavior unless their
dictionaries are changed.

## Deferred post-release work

The release still does not claim:

- formal asymptotic grid/time-step or CJ/experimental validation;
- long full-resolution/endurance equivalence;
- formal strong/weak scaling;
- active reversible-history restart or runtime redistribution/load balancing;
- a corrected H/H2 Soret formulation;
- validated physical H/H2 thermal-diffusion coefficients for arbitrary mechanisms.

See `qualification/POST_RELEASE_DEFERRED_QUALIFICATION.md`.
