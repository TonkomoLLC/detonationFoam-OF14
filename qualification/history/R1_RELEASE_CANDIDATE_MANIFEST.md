# R1 release-candidate manifest

Status: **candidate - runtime clean-build/install qualification required**

## Production source

- `applications/modules/detonationFluid/` - OpenFOAM 14 `foamRun` solver module.
- `src/detonationLegacyThermophysicalTransportModels/legacyMixtureAverageFourier/` - legacy mixture-average/property compatibility layer.
- `src/planarFvMeshTopoChangers/` - optional reusable true-2-D slab/wedge `fvMeshTopoChanger` library.

## Independent build paths

Full detonation package:

```bash
./Allwmake
```

AMR only, without compiling detonationFoam:

```bash
./AllwmakeAMR
```

or:

```bash
cd src/planarFvMeshTopoChangers
./Allwmake
```

## Documentation deliverables

- `docs/detonationFoam_OF14_Manual.md`
- `docs/detonationFoam_OF14_Manual.docx`
- `docs/detonationFoam_OF14_Manual.pdf`
- `docs/planarRefiner_OF14_Manual.md`
- `docs/planarRefiner_OF14_Manual.docx`
- `docs/planarRefiner_OF14_Manual.pdf`

The manuals use monospaced typography for OpenFOAM keywords, dictionaries/files, and shell commands. Markdown mathematics uses GitHub-compatible math delimiters. Word/PDF TOCs contain page numbers.

## R1 runtime qualification scope

The laptop-scale R1 gate performs:

1. package/source/documentation static audit;
2. clean independent AMR-only build;
3. clean full detonation/legacy-transport build;
4. linkage audit proving the detonation solver does not require `libplanarFvMeshTopoChangers.so`;
5. short packaged OF14 serial `foamRun` tutorial smoke;
6. runtime-log health checks and result bundling.

Default R1 hard budget: 900 s (15 minutes).

Long full-resolution equivalence, formal grid/CJ convergence, endurance, and formal scaling remain deferred post-release.
