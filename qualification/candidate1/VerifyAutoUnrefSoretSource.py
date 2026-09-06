#!/usr/bin/env python3
"""Static source gate for Candidate 1 automatic unrefinement + published OF8 Soret."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]

checks = []
def check(ok, msg):
    checks.append((bool(ok), msg))
    print(("PASS" if ok else "FAIL") + ": " + msg)

h = (ROOT / "src/planarFvMeshTopoChangers/planarRefiner/planarRefiner.H").read_text()
c = (ROOT / "src/planarFvMeshTopoChangers/planarRefiner/planarRefiner.C").read_text()
mdh = (ROOT / "src/planarFvMeshTopoChangers/planarRefiner/planarMultiDirRefinement.H").read_text()
mdc = (ROOT / "src/planarFvMeshTopoChangers/planarRefiner/planarMultiDirRefinement.C").read_text()
th = (ROOT / "src/detonationLegacyThermophysicalTransportModels/legacyMixtureAverageFourier/legacyMixtureAverageFourier.H").read_text()
tc = (ROOT / "src/detonationLegacyThermophysicalTransportModels/legacyMixtureAverageFourier/legacyMixtureAverageFourier.C").read_text()
manualPath = ROOT / "CANDIDATE1_AUTO_UNREF_SORET.md"
if not manualPath.exists():
    manualPath = ROOT / "qualification/history/autoUnrefSoret_candidate1/CANDIDATE1_AUTO_UNREF_SORET.md"
manual = manualPath.read_text()

# Default-off regression protection
check('lookupOrDefault<Switch>("automaticUnrefinement",false)' in c,
      "automatic planar unrefinement defaults OFF")
check('legacyThermalDiffusionMode_("off")' in tc and '"legacyThermalDiffusionMode",\n        "off"' in tc,
      "legacy Soret defaults OFF")

# Reversible cutter path
check('undoableMeshCutter(mesh,true)' in c,
      "planarRefiner creates an undo-enabled persistent cutter")
check('reversibleCutter_->getSplitFaces()' in c,
      "unrefinement candidates come from the reversible split tree")
check('reversibleCutter_->removeSplitFaces(splitFaces,meshMod)' in c.replace(' ', ''),
      "selected split faces are removed through undoableMeshCutter")
check('maxUnrefinementPassesPerUpdate",2' in c.replace(' ', ''),
      "default allows two binary undo passes per planar level")
check('fld[own] < lowerUnrefineLevel_' in c and 'fld[nei] < lowerUnrefineLevel_' in c
      and 'fld[own] > upperUnrefineLevel_' in c and 'fld[nei] > upperUnrefineLevel_' in c,
      "coarsening requires both siblings on the same outside side of the hysteresis band")
check('refinementLevel_' in h and '2*maxRefinementLevel_' in c,
      "per-cell binary cut depth limits repeated reversible refinement")
check('reversibleCutter_(),' in c and 'undoableMeshCutter& cutter' in mdh,
      "planarMultiDirRefinement accepts the externally owned reversible cutter")
check('refineAllDirs' in mdc and 'cutter,' in mdc,
      "external reversible cutter is passed into directional refinement")

# Candidate-1 safety limits
check('cannot map' in c and 'automaticUnrefinement v1.1.0' in c,
      "mesh-to-mesh mapping is refused with active reversible mode")
check('does not' in c and 'redistribute undoableMeshCutter ancestry' in c,
      "runtime mesh redistribution/load balancing is refused in reversible mode")
check('activeReversibleSplits' in c and 'does not serialize undoableMeshCutter history' in c,
      "restart state records/blocks active reversible ancestry")

# Published OF8 Soret semantics
check('legacyThermalDiffusionMode_ != "publishedOF8" || i != HIndex_' in tc,
      "direct thermal-diffusion ratio is returned only for H in publishedOF8 mode")
check('accumulate(HIndex_, thermalDiffH_);' in tc and 'accumulate(H2Index_, thermalDiffH2_);' in tc,
      "both H and H2 polynomial sums accumulate into the H ratio")
check('TDRatio_H2 remains zero' in tc,
      "source explicitly documents the published H2-zero assignment")
check('tFlux.ref() += tThermal();' in tc,
      "Soret raw flux is included in the total raw mixture diffusion flux")
check('-mwFlux - tThermal() + fvc::interpolate(Yi)*sumRaw' in tc,
      "species equation includes direct Soret plus zero-net-mass-flux correction")
check('ThermDiff_1' in tc and 'ThermDiff_4' in tc and 'thermalDiffusionCoeffs' in tc,
      "legacy four-coefficient thermal-diffusion polynomials are read per H/H2 pair")
check('"thermoDiff"' in tc and 'mesh.time().constant()' in tc
      and 'IOobject::MUST_READ_IF_MODIFIED' in tc,
      "publishedOF8 mode accepts the original constant/thermoDiff file layout")
check('does not infer the OF8' in tc and 'trandat/groupSpecies' in tc,
      "missing pair coefficients are not silently inferred from legacy transport groups")
check('legacy reproduction mode' in manual.lower() and 'no tutorial `constant/thermoDiff` coefficient file' in manual,
      "candidate notes distinguish reproduction from corrected physics and do not invent coefficients")

failed = [m for ok, m in checks if not ok]
if failed:
    print(f"\nFAILED: {len(failed)} of {len(checks)} static checks")
    sys.exit(1)
print(f"\nPASSED: {len(checks)} static checks")
