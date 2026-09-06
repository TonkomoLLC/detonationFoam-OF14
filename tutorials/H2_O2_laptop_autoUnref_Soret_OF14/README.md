# H2/O2 laptop smoke: automatic planar unrefinement + published OF8 Soret

This is a deliberately small OpenFOAM Foundation 14 integration smoke for
`foamRun` / `detonationFluid`.

It exercises, in one reacting run:

- H2/O2 finite-rate chemistry;
- `legacyMixtureAverageFourier` with `legacyThermalDiffusionMode publishedOF8`;
- `planarRefiner` refinement and automatic reversible unrefinement;
- field mapping through repeated topology change.

The base mesh contains only 120 cells. Run with:

```bash
./Allrun 2>&1 | tee log.laptopSmoke
```

`Allrun` verifies that `setFields` actually created the hot/high-pressure
ignition kernel before `foamRun` starts, then analyzes `log.foamRun` for at
least one refinement and one automatic-unrefinement event.

## Important Soret qualification note

The supplied OF8 source archive did not contain the physical `constant/thermoDiff`
coefficient dataset. This tutorial therefore uses one **synthetic diagnostic**
H2-O2 thermal-diffusion polynomial only to exercise the coupled runtime path.
It is not a quantitative Soret validation case and its coefficient must not be
interpreted as recommended transport data.

The quantitative compatibility check is packaged under
`qualification/candidate1Runtime/RunPublishedOF8SoretRuntimeGate`, which was
passed during v1.1.0 qualification and verifies the published OF8 H/H2 flux
identity and zero-net-mass-flux correction to roundoff.

## Qualified laptop run

The release qualification run reported:

```text
refinement events: 120->150, 150->153, 153->156, 156->159,
                   159->162, 162->165, 146->149, 144->147
unrefinement selected split counts: 6, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1
final T min/max: 300, 3581.900841 K
final shock position: 0.00051 m
```

The smoke is intentionally short and is not a formal CJ-speed, grid-convergence,
or long-duration detonation validation.
