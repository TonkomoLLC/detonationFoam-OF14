# Gate E4 status

**PASSED** — native OpenFOAM 14 localized AMR, MPI2 execution, Scotch
load-balancing, and in-place restart from a refined/redistributed processor
checkpoint are mutually qualified.

The checkpoint serial-vs-MPI comparison and the final continuous-vs-restart
comparison both satisfy the E4 geometry, p/T/U, peak-state, shock-location,
all-species, species-bound, and sum(Y) criteria. Species with peak mass
fraction below 1e-6 are evaluated using the strict absolute trace-species
criterion because normalized errors are ill-conditioned near zero.

No production `detonationFluid` or transport source was changed for E4.
