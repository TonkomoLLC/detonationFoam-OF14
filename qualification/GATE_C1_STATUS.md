# Gate C1 status — CLOSED / PASSED

Gate C1 is formally **CLOSED / PASSED** from the returned OpenFOAM Foundation 14 runtime bundle.

All seven flux schemes compiled from one solver build and passed the identical 2,400-cell / 12-mm serial detonation qualification:

- Kurganov
- Tadmor
- HLL
- HLLC
- HLLCP
- AUSM+
- AUSM+up

Every scheme reached `0.0102675 m` versus the `0.0102625 m` propagation target and terminated normally. The returned cross-scheme relative spans were approximately 1.04% in arrival time, 1.15% in `pmax`, and 1.60e-5 in `Tmax`.

`GateC1_serial_baseline.csv` is the exact returned serial baseline carried into Gate C2.
