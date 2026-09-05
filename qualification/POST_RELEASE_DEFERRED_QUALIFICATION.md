# Post-release deferred qualification

The pre-release program intentionally targets laptop-friendly gates of about 30 minutes or less. The following are useful follow-up studies but are **not release blockers solely because of runtime**:

- long/full-resolution detonation propagation and endurance trajectories;
- formal asymptotic mesh/time-step convergence and CJ-speed/thermochemical validation;
- broad mesh/time-step/chemistry parameter sweeps;
- full-resolution OF8↔OF14 trajectory/profile equivalence;
- formal strong/weak scaling over larger meshes and more MPI ranks;
- long-duration native load-balancer performance/amortization studies;
- extended AMR + load-balance + restart endurance and repeated refinement cycles.

Release documentation should state these deferred items explicitly and distinguish compact pre-release qualification from extended post-release validation.
