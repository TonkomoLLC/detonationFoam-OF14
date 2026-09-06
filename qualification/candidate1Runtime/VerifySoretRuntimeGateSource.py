#!/usr/bin/env python3
from pathlib import Path
import re, sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SRC=ROOT/'src/detonationLegacyThermophysicalTransportModels/legacyMixtureAverageFourier/legacyMixtureAverageFourier.C'
APP=HERE/'applications/soretFluxQualification/soretFluxQualification.C'
ON=HERE/'cases/soretPublishedOF8'
OFF=HERE/'cases/soretOff'
checks=[]
def ck(cond,msg):
    if not cond: print('FAIL:',msg); sys.exit(1)
    print('PASS:',msg); checks.append(msg)
s=SRC.read_text(); a=APP.read_text(); td=(ON/'constant/thermoDiff').read_text(); tr=(ON/'constant/thermophysicalTransport').read_text(); troff=(OFF/'constant/thermophysicalTransport').read_text()
ck('i != HIndex_' in s, 'production direct thermal ratio remains H-only')
ck('accumulate(HIndex_, thermalDiffH_)' in s and 'accumulate(H2Index_, thermalDiffH2_)' in s, 'production H ratio receives both H and H2 polynomial sums')
ck('thermo.Y("H")' in a and 'thermo.Y("H2")' in a and 'thermophysicalTransport->j(YH)' in a and 'thermophysicalTransport->j(YH2)' in a, 'runtime probe uses OF14 named-species fields and public species-flux API')
ck('(1.0-yh)*h2 + yh2*h' in a, 'runtime probe checks published H/H2 corrected-flux identity')
ck('inferredRawH2' in a, 'runtime probe independently infers direct H2 contribution')
ck('legacyThermalDiffusionMode publishedOF8;' in tr, 'published control case enables compatibility mode')
ck('legacyThermalDiffusionMode off;' in troff, 'control case disables Soret')
ck((ON/'system/fvSolution').is_file() and (OFF/'system/fvSolution').is_file(), 'both runtime cases provide the OF14 fvSolution dictionary required during model construction')
ck(td.count('ThermDiff_1')==63, 'all explicit H/H2 pair dictionaries are present')
vals=[]
for m in re.finditer(r'([A-Za-z0-9_()+.-]+-[A-Za-z0-9_()+.-]+)\s*\{\s*ThermDiff_1\s+([^;]+);',td,re.S):
    vals.append((m.group(1),float(m.group(2))))
nonzero=[x for x in vals if abs(x[1])>0]
ck(nonzero==[('H2-N2',100.0)], 'only the H2-N2 polynomial is nonzero in the diagnostic case')
ck(tr.count('type legacyBinaryDiffusionCoefficient;')==528, 'complete 33-species binary-diffusion table is present')
print(f'PASS: {len(checks)}/{len(checks)} Soret runtime-gate source checks')
