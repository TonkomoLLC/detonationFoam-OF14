# detonationFoam for OpenFOAM Foundation 14 - v1.1.0

This package is the OpenFOAM Foundation 14 port of the OpenFOAM 8
`detonationFoam` code line. The reacting density-based solver is delivered as
the modular `foamRun` solver module `detonationFluid`.

Version 1.1.0 adds two opt-in capabilities that were requested after v1.0.0:

- **automatic unrefinement for the reusable 2-D/wedge `planarRefiner`**; and
- **explicit reproduction of the published OF8 H/H2 Soret behavior** through
  `legacyThermalDiffusionMode publishedOF8`.

Both features default to **off**, so v1.0.0-style cases retain their prior
behavior unless the new options are selected.

## Important scope note

This release was qualified primarily with compact laptop-friendly gates. It is
appropriate to perform case-specific mesh/time-step sensitivity and physical
validation before using a new mechanism, geometry, or operating condition for
production conclusions. Long full-resolution/endurance studies remain listed
as deferred qualification rather than being hidden behind the release label.

## Quick build

Source OpenFOAM Foundation 14 and build the detonation solver plus legacy
transport compatibility library:

```bash
./Allwmake
```

Run cases with:

```bash
foamRun
```

## Optional reusable 2-D/wedge AMR

The solver-independent AMR component is in:

```text
src/planarFvMeshTopoChangers/
```

Its runtime type is `planarRefiner`. It is an OpenFOAM
`fvMeshTopoChanger`, not an `fvModel`, and has no dependency on
`detonationFluid`.

Build only the AMR library with:

```bash
./AllwmakeAMR
```

This produces:

```text
$FOAM_USER_LIBBIN/libplanarFvMeshTopoChangers.so
```

### Automatic unrefinement

Automatic coarsening is opt-in. A representative `dynamicMeshDict` block is:

```text
topoChanger
{
    type                         planarRefiner;
    libs                         ("libplanarFvMeshTopoChangers.so");
    geometry                     slab;
    refineInterval               1;
    field                        T;
    lowerRefineLevel             800;
    upperRefineLevel             4900;
    maxCells                     200000;

    automaticUnrefinement        true;
    unrefineInterval             1;
    lowerUnrefineLevel           700;
    upperUnrefineLevel           5000;
    maxRefinementLevel           1;
    maxUnrefinementPassesPerUpdate 2;
}
```

Coarsening uses reversible split ancestry and a hysteresis band. With
`automaticUnrefinement true`, start from the base/unrefined mesh. Restart from
a time containing active reversible splits and runtime mesh redistribution/load
balancing are deliberately refused because reversible cutter ancestry is not
serialized/distributed by this implementation. Fixed-decomposition MPI2 was
qualified, and restart from a fully coarsened state was qualified.

## Published OF8 H/H2 Soret compatibility

The legacy mixture-average transport model defaults to:

```text
legacyThermalDiffusionMode off;
```

To reproduce the **published OF8 source behavior**, use:

```text
legacyThermalDiffusionMode publishedOF8;
```

This mode intentionally reproduces the OF8 assignment in which both the H and
H2 thermal-diffusion polynomial sums accumulate into the H thermal-diffusion
ratio while the direct H2 ratio remains zero. It is provided for
reproducibility and must not be interpreted as a corrected modern H/H2 Soret
model.

Coefficients can be supplied in the original OF8 `constant/thermoDiff` layout
or under `thermalDiffusionCoeffs`. Explicit H/H2 pair entries are required; the
legacy `trandat`/`groupSpecies` coefficient-sharing shortcut is not inferred.

The supplied OF8 source archive did not contain a physical `thermoDiff`
coefficient dataset, so this release does not invent one.

## Small integrated H2/O2 laptop smoke

A 120-cell tutorial is included specifically to exercise the two v1.1.0
features together:

```bash
cd tutorials/H2_O2_laptop_autoUnref_Soret_OF14
./Allrun 2>&1 | tee log.laptopSmoke
```

The accepted qualification run refined and automatically unrefined repeatedly,
reached a final maximum temperature of about 3582 K, and advanced the leading
shock to 0.51 mm without a fatal error/FPE. The tutorial uses one explicitly
**synthetic diagnostic Soret coefficient** because physical OF8 `thermoDiff`
data were not supplied; it is an integration smoke, not quantitative Soret
validation.

## Documentation

The main solver manual and the standalone reusable-AMR manual are supplied in
GitHub Markdown, Word, and PDF:

```text
docs/detonationFoam_OF14_Manual.md
docs/detonationFoam_OF14_Manual.docx
docs/detonationFoam_OF14_Manual.pdf

docs/planarRefiner_OF14_Manual.md
docs/planarRefiner_OF14_Manual.docx
docs/planarRefiner_OF14_Manual.pdf
```

The manuals cover theory, usage, OF8-to-OF14 case migration, capability
classification, automatic unrefinement, published-OF8 Soret compatibility,
parallel/restart constraints, qualification, limitations, and troubleshooting.

## Release qualification

The authoritative status is:

```text
QUALIFICATION_SUMMARY.md
```

The v1.1.0 feature evidence is summarized in:

```text
qualification/V1_1_AUTO_UNREF_SORET_QUALIFICATION.md
qualification/INTEGRATED_H2_LAPTOP_SMOKE.md
```

The detailed historical runtime harnesses that produced the accepted feature
evidence remain in `qualification/candidate1Runtime/` for reproducibility.

To run the compact release audit/build/runtime smoke on another OF14
installation:

```bash
./VerifyRelease
./RunReleaseQualification
```

The release qualification remains laptop-oriented; long full-resolution,
formal convergence/CJ validation, endurance, corrected-Soret development, and
reversible-history serialization/distribution remain explicitly deferred.
