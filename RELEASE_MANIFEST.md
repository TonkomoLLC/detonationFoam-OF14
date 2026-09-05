# detonationFoam OF14 v1.0.0 release manifest

## Core production components

- `applications/modules/detonationFluid/` - modular OpenFOAM 14 `foamRun` solver.
- `src/detonationLegacyThermophysicalTransportModels/legacyMixtureAverageFourier/` - legacy mixture-average/property compatibility layer.
- `src/planarFvMeshTopoChangers/` - optional reusable 2-D slab/wedge `fvMeshTopoChanger` library.

## Build entry points

Full detonation package:

```bash
./Allwmake
```

Standalone AMR only:

```bash
./AllwmakeAMR
```

or:

```bash
cd src/planarFvMeshTopoChangers
./Allwmake
```

## Manuals

- `docs/detonationFoam_OF14_Manual.md`
- `docs/detonationFoam_OF14_Manual.docx`
- `docs/detonationFoam_OF14_Manual.pdf`
- `docs/planarRefiner_OF14_Manual.md`
- `docs/planarRefiner_OF14_Manual.docx`
- `docs/planarRefiner_OF14_Manual.pdf`

## Tutorials

- `tutorials/1D_NH3_O2_cracking_0.3_detonation_OF14/`
- `tutorials/1D_NH3_O2_cracking_0.3_detonation_OF14_fast/`
- `tutorials/1D_NH3_O2_cracking_0.3_detonation_OF8_fast_reference/` - migration/equivalence reference data, not an OF14 executable path.

## Qualification

`QUALIFICATION_SUMMARY.md` is the authoritative release-level status. Files under `qualification/` preserve detailed gate evidence and historical staged-development records. Intermediate candidate status files are retained for traceability and are not the authority for final release state.
