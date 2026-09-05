# Gate E3 status

**PASSED** — the localized native OpenFOAM 14 AMR solution agrees with a
same-finest-spacing 19200-cell fixed reference within the E3 pressure,
temperature, velocity, shock-location, species-profile, species-integral, and
local closure tolerances while using fewer cells. Localized hexRef8 transition
polyhedra are explicitly characterized and accepted only if all non-concavity
mesh-quality checks remain good.

No production `detonationFluid` or transport source was changed for E3.
