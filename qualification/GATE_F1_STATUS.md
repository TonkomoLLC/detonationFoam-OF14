# Gate F1 status

**PASSED** — OpenFOAM 14 contains the low-level directional mesh-cutting and
2-D geometry support needed to pursue a reusable planar/axisymmetric dynamic
AMR component.

The selected architecture is a runtime-selectable `fvMeshTopoChanger`
(`planarRefiner`) backed by `multiDirRefinement`/`meshCutter`/`cellLooper`,
not an `fvModel`.  The component is intended to support both slab/empty and
wedge/axisymmetric meshes and to be reusable by compatible OF14 solvers that
participate in the normal `mesh.update()` lifecycle.

No production solver or mesh-refinement source is changed by F1.
