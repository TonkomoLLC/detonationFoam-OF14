# detonationFoam OF14 Candidate 1: automatic planar unrefinement + published OF8 Soret compatibility

This candidate adds two opt-in features without changing the default behavior of the existing OF14 port.

## 1. `planarRefiner` automatic unrefinement

Enable with:

```text
automaticUnrefinement          true;
unrefineInterval               1;
lowerUnrefineLevel             700;
upperUnrefineLevel             5000;
maxRefinementLevel             1;
maxUnrefinementPassesPerUpdate 2;
```

For a refinement band of `800..4900`, the example above refines inside that band but coarsens only when both sibling cells are below 700 or both are above 5000. This creates hysteresis.

Implementation notes:

- the existing directional mesh cutter is retained;
- a persistent `undoableMeshCutter(mesh, true)` owns the reversible split tree;
- one complete planar level is two binary cuts, so up to two undo passes are used by default;
- a runtime `planarRefinementLevel` field tracks binary cut depth and is mapped with mesh topology changes;
- `maxRefinementLevel` replaces the old global-event cap for reversible mode, allowing repeated refine/unrefine cycles;
- the default `automaticUnrefinement false` preserves the previously qualified refinement-only behavior.

### Candidate-1 safety limits

OpenFOAM's `undoableMeshCutter` history is not serialized by this component. Consequently:

1. restart is refused when the prior write recorded active reversible splits;
2. mesh-to-mesh mapping and runtime redistribution/load balancing are refused while automatic unrefinement is enabled;
3. if `removeSplitFaces()` asks to remove more faces than explicitly selected, the run is stopped instead of accepting ambiguous coarsening bookkeeping.

These restrictions are intentional until reversible ancestry can be serialized/distributed and qualified.

Enable `automaticUnrefinement` from the base/unrefined mesh. Candidate 1 cannot adopt reversible ancestry for cells that were already refined by an earlier refinement-only run or by another topology changer.

## 2. Published OF8 H/H2 Soret behavior

Enable in the `legacyMixtureAverageFourier` coefficient dictionary with:

```text
legacyThermalDiffusionMode publishedOF8;

thermalDiffusionCoeffs
{
    H-H2
    {
        ThermDiff_1  ...;
        ThermDiff_2  ...;
        ThermDiff_3  ...;
        ThermDiff_4  ...;
    }

    // Supply H-X and H2-X pairs for every other species X.
}
```

The same coefficients may instead be supplied in the original OF8 location `constant/thermoDiff`; when the inline `thermalDiffusionCoeffs` dictionary is absent, Candidate 1 reads that file automatically. Candidate 1 requires explicit H/H2 pair entries and does not infer the old `trandat`/`groupSpecies` coefficient-sharing fallback.

The compatibility mode intentionally reproduces the source-level OF8 assignment:

```text
H polynomial contributions  -> TDRatio_H
H2 polynomial contributions -> TDRatio_H
TDRatio_H2                   -> remains zero
```

This is therefore a **legacy reproduction mode**, not a corrected physical Soret model. It is named explicitly so a later corrected mode can coexist without changing historical results.

The Soret raw flux is included in the same OF8-style mixture correction used for molecular-weight diffusion, preserving the zero-net-mixture-mass-flux construction.

The supplied OF8 source archive contains the thermal-diffusion reader and equations but no tutorial `constant/thermoDiff` coefficient file. This candidate therefore does not invent coefficient values.

## 3. Qualification status

This candidate has source/static checks in the packaged tree but has **not** been compiled or run in a real OpenFOAM 14 environment in this ChatGPT sandbox. The next gate should include:

1. full `./Allwmake` and `./AllwmakeAMR` compile;
2. refinement-only regression with `automaticUnrefinement false`;
3. slab refine -> unrefine -> re-refine smoke with exact cell-count checks and final `checkMesh`;
4. wedge refine -> unrefine smoke and final `checkMesh`;
5. serial/MPI (fixed decomposition, no runtime redistribution) reversible AMR comparison;
6. Soret-off regression against the existing D3C transport gate;
7. `publishedOF8` unit/smoke case checking that `TDRatio_H` receives both H and H2 polynomial sums while the direct H2 ratio remains zero;
8. species-sum and energy-flux conservation checks with Soret enabled.
