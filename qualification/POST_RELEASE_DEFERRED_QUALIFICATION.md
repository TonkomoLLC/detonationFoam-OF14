# Post-release deferred qualification

The release program intentionally targets laptop-friendly gates of about
30 minutes or less. The following are useful follow-up studies but are **not
release blockers solely because of runtime**:

- long/full-resolution detonation propagation and endurance trajectories;
- formal asymptotic mesh/time-step convergence and CJ-speed/thermochemical validation;
- broad mesh/time-step/chemistry parameter sweeps;
- full-resolution OF8<->OF14 trajectory/profile equivalence;
- formal strong/weak scaling over larger meshes and more MPI ranks;
- long-duration native load-balancer performance/amortization studies;
- long endurance with repeated refinement/unrefinement cycles;
- serialization and MPI redistribution of active reversible `planarRefiner`
  ancestry, so automatic unrefinement can eventually coexist with active-history
  restart and runtime redistribution/load balancing;
- a separately derived and validated **corrected** H/H2 Soret formulation,
  distinct from the intentionally bug-compatible `publishedOF8` reproduction mode;
- physical validation of real H/H2 thermal-diffusion coefficients for target
  mechanisms/cases (the packaged smoke/qualification coefficients are synthetic
  diagnostics only).

Release documentation distinguishes these deferred items from the compact
pre-release qualification that has already passed.
