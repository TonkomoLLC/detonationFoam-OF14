# Integrated H2/O2 laptop smoke - accepted result

Status: **PASSED on OpenFOAM Foundation 14** (2026-09-06).

Case: `tutorials/H2_O2_laptop_autoUnref_Soret_OF14/`

The accepted user run exercised `foamRun` / `detonationFluid`, the
`publishedOF8` H/H2 Soret compatibility path, and reversible `planarRefiner`
automatic unrefinement in one short calculation.

Accepted diagnostics:

```text
PASS: publishedOF8 Soret mode constructed
PASS: planarRefiner automatic unrefinement constructed
PASS: at least one planar refinement event occurred
PASS: at least one automatic unrefinement pass occurred
PASS: no OpenFOAM fatal/FPE in solver log
PASS: solver reached a normal stop/end condition

refinement events:
120->150, 150->153, 153->156, 156->159, 159->162, 162->165, 146->149, 144->147

unrefinement selected split counts:
6, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1

final T min/max: 300, 3581.900841 K
final shock position: 0.00051 m
```

The original r2 initialization verifier printed `120` as the minimum of both
`T` and `p` because it accidentally included the `List<scalar>` element count
in the numeric scan. This was a **post-setFields parser artifact**, not a field
value. The release tutorial fixes the parser by reading only the list body.
The solver log independently reported a final minimum temperature of 300 K.
