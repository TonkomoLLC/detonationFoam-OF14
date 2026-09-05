# Gate A status

## A1 — build and modular startup: PASS

`libdetonationFluidSolver.so` builds in OpenFOAM Foundation 14 and is selected
through `foamRun` as `solver detonationFluid`.

## A2-fast — PASS by superseding C1/C2 evidence

The reduced 2400-cell, 5-micron case has now been qualified more strongly than
the original A2 criterion:

- Gate C1: all seven flux schemes pass the serial propagation matrix.
- Gate C2: all seven schemes pass on four MPI ranks with repeated native
  OpenFOAM-14 load balancing/redistribution and reproduce the serial results.

The earlier 20,000-cell run also demonstrated stable full-domain propagation.

## A3-fast — pending final OF8/OF14 comparison

Gate C3 carries the original fixed-time OF8-vs-OF14 profile comparator forward.
Gate A closes when that actual cross-version comparison passes.
