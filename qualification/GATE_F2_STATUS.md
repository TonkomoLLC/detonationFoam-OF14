# Gate F2 status

**PASSED** — the reusable OpenFOAM 14 `planarRefiner` fvMeshTopoChanger compiled as `libplanarFvMeshTopoChangers.so` and dynamically refined a true slab/`empty` 2-D detonation case while preserving registered fields and the qualified D3-B transport/property stack.

The F2 implementation uses OF14's directional mesh-cutting machinery with an adapted iterator that invokes `fvMesh::preChange()` and `fvMesh::topoChange(map)` after each directional split. F2 does not yet qualify wedge/axisymmetric geometry, persistent refinement history/unrefinement, or MPI restart; those remain F3.
