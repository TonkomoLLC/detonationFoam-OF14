# Gate C2 status — CLOSED / PASSED

Gate C2 is closed on returned candidate-4 evidence.

The seven restored flux schemes all passed on 4 MPI ranks with native
OpenFOAM-14 chemistry-aware load balancing and forced repeated redistribution.
Each scheme reproduced its Gate-C1 serial arrival time exactly at the reported
precision; pmax differences were only O(1e-8) to O(1e-7) relative.

The qualification exercised 28–40 real `Redistributing mesh` events per
scheme. This also qualifies the C2 redistribution fixes:

1. `multiConstraint false` for Scotch-compatible summed cell weights;
2. the zero-sized `type internal` redistribution polyPatch; and
3. clearing topology-dependent face reconstruction/flux caches whenever
   `mesh.distributing()` is active before `mesh_.update()`.
