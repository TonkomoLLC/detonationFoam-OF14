# Candidate 1 Runtime Gate 1a harness correction

This overlay changes only the qualification harness. Production sources are unchanged.

The original Runtime Gate 1 launched the two-rank qualification driver with MPI, but then invoked
`checkMesh -parallel -time 4` directly as a single process. For a two-way decomposed case the
parallel mesh check must itself be launched with two MPI ranks.

Corrected invocation:

```bash
mpirun -np 2 checkMesh -parallel -time 4
```

The returned Runtime Gate 1 results already passed the serial uniform slab, localized slab, wedge,
and the MPI2 refinement/unrefinement cycle plus integral-conservation analyzer. Gate 1a therefore
corrects a test-harness defect; it does not change `planarRefiner` or `detonationFluid`.
