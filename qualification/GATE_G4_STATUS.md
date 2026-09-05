# Stage G4 status — PASSED

Stage G4 candidate 2 passed the compact serial/MPI performance and native load-balancing smoke.

- Serial: 152.966 s
- MPI2: 91.855 s; speedup 1.665; efficiency 0.833
- MPI4: 102.875 s; speedup 1.487; efficiency 0.372
- MPI2/MPI4 each observed 29 native OpenFOAM 14 load-balancer redistributions using the Scotch distributor.
- Final shock location, peak pressure, and peak temperature were identical across serial/MPI2/MPI4: shock = 0.0102625 m, pmax = 2,158,684.42 Pa, Tmax = 4994.47644 K.

Performance is diagnostic only. Formal scaling remains deferred post-release.
