# R1 returned-result assessment - PASSED with harness-only false positive

The returned R1 run passed the substantive qualification steps: standalone AMR clean build, full detonation/legacy-transport clean build, link-independence audit, and the packaged fast serial `foamRun` command.

The runner then reported `ERROR: fast tutorial runtime health check failed` because its post-run grep matched OpenFOAM's normal startup message:

```text
sigFpe : Enabling floating point exception trapping (FOAM_SIGFPE).
```

This is not a floating-point exception. It announces that floating-point exception trapping is enabled. The same failure mode had previously been documented in Gate B1.

The smoke command itself returned success; otherwise the runner's `run_step` wrapper would have failed at the `packaged OF14 fast serial foamRun smoke` step before reaching the post-run grep.

Disposition: **R1 PASSED; release is not blocked.** The final `RunReleaseQualification` checker accepts the normal `FOAM_SIGFPE` startup line while still rejecting actual `FOAM FATAL`, segmentation-fault, core-dump, and true floating-point-exception lines.
