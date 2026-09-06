# detonationFoam OF14 v1.1.0 release manifest

## Core production components

- `applications/modules/detonationFluid/` - modular OpenFOAM 14 `foamRun` solver.
- `src/detonationLegacyThermophysicalTransportModels/legacyMixtureAverageFourier/` - legacy mixture-average/property compatibility layer plus opt-in `publishedOF8` H/H2 Soret reproduction.
- `src/planarFvMeshTopoChangers/` - optional reusable 2-D slab/wedge `fvMeshTopoChanger` with opt-in reversible automatic unrefinement.

## Build entry points

Full detonation package:

```bash
./Allwmake
```

Standalone AMR only:

```bash
./AllwmakeAMR
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
- `tutorials/H2_O2_laptop_autoUnref_Soret_OF14/` - 120-cell integrated v1.1.0 smoke; synthetic Soret diagnostic coefficient.
- `tutorials/1D_NH3_O2_cracking_0.3_detonation_OF8_fast_reference/` - migration/equivalence reference data, not an OF14 executable path.

## Qualification

- `QUALIFICATION_SUMMARY.md` - authoritative release-level status.
- `qualification/V1_1_AUTO_UNREF_SORET_QUALIFICATION.md` - isolated v1.1.0 feature gates.
- `qualification/INTEGRATED_H2_LAPTOP_SMOKE.md` - accepted coupled reacting smoke.
- `qualification/candidate1Runtime/` - exact historical runtime harnesses used to produce the accepted automatic-unrefinement/Soret evidence.
- `qualification/POST_RELEASE_DEFERRED_QUALIFICATION.md` - intentionally deferred long/advanced work.

Intermediate candidate records are retained under `qualification/history/` for traceability and are not authoritative release status.
