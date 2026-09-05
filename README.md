# detonationFoam for OpenFOAM Foundation 14 - v1.0.0

This package is the released OpenFOAM Foundation 14 port of the OpenFOAM 8 `detonationFoam` code line. The detonation solver is delivered as the modular `foamRun` solver module `detonationFluid`.

# Important

This was largely a porting exercise with testing on smaller cases. Care should be taken to confirm that this works for larger cases.

The 2D AMR is available as an fvModel and can be used by other non-`detonationFoam` solvers with OpenFOAM 14. 
See below for more details.


## Quick build

Source OpenFOAM Foundation 14, then build the detonation solver and legacy transport compatibility library:

```bash
./Allwmake
```

Run cases with:

```bash
foamRun
```

## Optional reusable 2-D/wedge AMR - independent of detonationFoam

The reusable planar/axisymmetric AMR component is in:

```text
src/planarFvMeshTopoChangers/
```

Its runtime type is `planarRefiner`. It is an OpenFOAM `fvMeshTopoChanger`, not an `fvModel`, and has no dependency on `detonationFluid`.

Build only the AMR library from the release root:

```bash
./AllwmakeAMR
```

or directly:

```bash
cd src/planarFvMeshTopoChangers
./Allwmake
```

This builds:

```text
$FOAM_USER_LIBBIN/libplanarFvMeshTopoChangers.so
```

without compiling the detonation solver.

The normal `./Allwmake` path does not require or link the optional AMR library.

## Documentation

The manuals are supplied in Markdown, Word, and PDF:

```text
docs/detonationFoam_OF14_Manual.md
docs/detonationFoam_OF14_Manual.docx
docs/detonationFoam_OF14_Manual.pdf

docs/planarRefiner_OF14_Manual.md
docs/planarRefiner_OF14_Manual.docx
docs/planarRefiner_OF14_Manual.pdf
```

The main manual includes theory, usage, OF8-to-OF14 migration, a capability-disposition table, native OpenFOAM substitutions, omitted/deferred capability, AMR, MPI/load balancing, restart, qualification results, limitations, and troubleshooting. OpenFOAM keywords, dictionaries/files, runtime types, paths, and shell commands are monospaced in Word/PDF and use backticks/code blocks in Markdown. Markdown mathematics uses GitHub-compatible `$...$` and `$$...$$` syntax. Word/PDF tables of contents include page numbers.

## Release qualification

OpenFOAM 14 v1.0.0 is released after the staged A-G migration and R1 clean-build/runtime gate. See:

```text
QUALIFICATION_SUMMARY.md
qualification/POST_RELEASE_DEFERRED_QUALIFICATION.md
```

The final R1 returned run clean-built the standalone AMR library, clean-built the full detonation/legacy-transport libraries, confirmed no runtime link dependency from `libdetonationFluidSolver.so` to `libplanarFvMeshTopoChangers.so`, and successfully returned from the packaged short `foamRun` tutorial. The only R1 failure message was a post-run checker false positive caused by matching OpenFOAM's normal `FOAM_SIGFPE` startup line; that checker is corrected in `RunReleaseQualification`.

To reproduce the compact release qualification on another OpenFOAM Foundation 14 installation:

```bash
./VerifyRelease
./RunReleaseQualification
```

The default runtime budget is 900 s (15 minutes). Long full-resolution equivalence, formal grid/CJ convergence, endurance, and formal strong/weak scaling remain explicitly deferred post-release.
