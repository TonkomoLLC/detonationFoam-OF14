# Reusable OpenFOAM 14 planar / axisymmetric AMR

This directory builds `libplanarFvMeshTopoChangers.so`, providing the runtime-selectable `planarRefiner` `fvMeshTopoChanger`.

It is **solver-independent** and has no dependency on `detonationFluid`.

Build only this component:

```bash
cd src/planarFvMeshTopoChangers
./Allwmake
```

or from the release root:

```bash
./AllwmakeAMR
```

This component is an `fvMeshTopoChanger`, **not an `fvModel`**. It owns 2-D slab/empty and wedge/axisymmetric topology refinement. Current release scope is refinement-only; automatic unrefinement/coarsening is not included.

See `docs/planarRefiner_OF14_Manual.md`, `.docx`, or `.pdf`.
