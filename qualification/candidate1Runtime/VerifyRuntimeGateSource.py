#!/usr/bin/env python3
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
q=Path(__file__).resolve().parent
checks=[]
def ck(cond,msg): checks.append((bool(cond),msg))
prod=(ROOT/'src/planarFvMeshTopoChangers/planarRefiner/planarRefiner.C').read_text()
h=(q/'applications/planarRefinerQualification/planarRefinerQualification.C').read_text()
ck('automaticUnrefinement_' in prod,'production candidate contains automatic unrefinement')
ck('QUAL_METRIC' in h and 'mesh.update()' in h,'qualification driver reports metrics around mesh.update')
ck('uniformTracer' in h and 'tracerIntegral' in h,'field-mapping conservation probes present')
ck((q/'cases/slabUniform/system/blockMeshDict').exists(),'uniform slab case present')
ck((q/'cases/slabLocalized/system/blockMeshDict').exists(),'localized slab case present')
ck((q/'cases/wedgeUniform/system/blockMeshDict').exists(),'annular wedge case present')
ck('automaticUnrefinement true;' in (q/'cases/slabUniform/constant/dynamicMeshDict').read_text(),'automatic unrefinement enabled in runtime case')
ck('maxUnrefinementPassesPerUpdate 2;' in (q/'cases/slabUniform/constant/dynamicMeshDict').read_text(),'two undo passes configured')
ck('(50 24 1)' in (q/'cases/slabUniform/system/blockMeshDict').read_text(),'slab base mesh is 1200 cells')
ck('(80 20 1)' in (q/'cases/wedgeUniform/system/blockMeshDict').read_text(),'wedge base mesh is 1600 cells')
for ok,msg in checks: print(('PASS: ' if ok else 'FAIL: ')+msg)
if not all(x for x,_ in checks): sys.exit(1)
print(f'PASS: {sum(x for x,_ in checks)}/{len(checks)} runtime-gate source checks')
