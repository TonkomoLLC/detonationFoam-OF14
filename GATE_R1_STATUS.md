# R1 release qualification status - PASSED

R1 is accepted as passed for OpenFOAM Foundation 14.

The returned qualification clean-built the standalone planar AMR library, clean-built the detonation/legacy-transport libraries, confirmed that `libdetonationFluidSolver.so` does not depend on `libplanarFvMeshTopoChangers.so`, and successfully returned from the packaged fast serial `foamRun` smoke command.

The original final error was a harness-only false positive caused by matching OpenFOAM's normal `FOAM_SIGFPE` startup line as though it were a runtime exception. The final release checker corrects that logic. See `qualification/R1_RETURNED_RESULT_ASSESSMENT.md`.
