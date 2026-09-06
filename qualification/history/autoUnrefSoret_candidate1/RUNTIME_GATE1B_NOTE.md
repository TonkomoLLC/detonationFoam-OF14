# Candidate 1 Runtime Gate 1b

This overlay changes qualification code only. Production detonationFluid,
legacy transport, and planarRefiner sources are unchanged.

## Reason

OpenFOAM 14 `multicomponentThermo::specieIndex()` accepts a
`volScalarField`, not a species-name string.  The Soret qualification probe
now obtains the H and H2 fields with `thermo.Y("H")` and `thermo.Y("H2")`
and passes those fields directly to the public thermophysical-transport
`j(Yi)` API.

## Recommended rerun

The automatic-unrefinement runtime gate already passed. Run only:

```bash
./qualification/candidate1Runtime/RunPublishedOF8SoretRuntimeGate \
    2>&1 | tee log.Candidate1.Soret.Gate1b
```
