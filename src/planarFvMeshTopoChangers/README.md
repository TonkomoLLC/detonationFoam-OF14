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

This component is an `fvMeshTopoChanger`, **not an `fvModel`**. It owns 2-D slab/empty and wedge/axisymmetric topology refinement.

OpenFOAM 14 v1.1.0 adds optional automatic unrefinement using OpenFOAM's reversible `undoableMeshCutter` history. The default remains refinement-only (`automaticUnrefinement false`). When reversible coarsening is enabled, runtime mesh redistribution/load balancing and restart from a time containing active reversible splits are deliberately refused because the cutter ancestry is not serialized by this component.

See `docs/planarRefiner_OF14_Manual.md`, `.docx`, or `.pdf`.
