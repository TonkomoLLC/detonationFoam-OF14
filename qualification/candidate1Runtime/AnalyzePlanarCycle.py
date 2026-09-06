#!/usr/bin/env python3
from __future__ import annotations
import argparse, math, re, sys
from pathlib import Path

ap=argparse.ArgumentParser()
ap.add_argument('log',type=Path)
ap.add_argument('--expected',required=True,help='comma-separated post-update cell counts')
ap.add_argument('--rel-tol',type=float,default=2e-10)
a=ap.parse_args()
text=a.log.read_text(errors='replace')
pat=re.compile(r'QUAL_METRIC tag=(afterRefine|afterUnrefine) step=(\d+) cells=(\d+) volume=([^ ]+) uniformIntegral=([^ ]+) tracerIntegral=([^ ]+) uniformMin=([^ ]+) uniformMax=([^\s]+)')
rows=[]
for m in pat.finditer(text):
    rows.append((m.group(1),int(m.group(2)),int(m.group(3)),*[float(m.group(i)) for i in range(4,9)]))
exp=[int(x) for x in a.expected.split(',')]
if len(rows)!=len(exp):
    print(f'FAIL: expected {len(exp)} post-update metric rows, found {len(rows)}')
    sys.exit(1)
counts=[r[2] for r in rows]
if counts!=exp:
    print(f'FAIL: cell sequence {counts} != expected {exp}')
    sys.exit(1)
base_match=re.search(r'QUAL_METRIC tag=base step=\d+ cells=(\d+) volume=([^ ]+) uniformIntegral=([^ ]+) tracerIntegral=([^ ]+) uniformMin=([^ ]+) uniformMax=([^\s]+)',text)
if not base_match:
    print('FAIL: missing base metrics'); sys.exit(1)
base=[float(base_match.group(i)) for i in range(2,7)]
baseV,baseUI,baseTI,_,_=base

def relerr(x,y): return abs(x-y)/max(abs(y),1e-300)
for r in rows:
    tag,step,cells,V,UI,TI,uMin,uMax=r
    if relerr(V,baseV)>a.rel_tol:
        print(f'FAIL: volume drift at step {step}: {V} vs {baseV}'); sys.exit(1)
    if relerr(UI,baseUI)>a.rel_tol:
        print(f'FAIL: uniform integral drift at step {step}: {UI} vs {baseUI}'); sys.exit(1)
    if relerr(TI,baseTI)>a.rel_tol:
        print(f'FAIL: tracer integral drift at step {step}: {TI} vs {baseTI}'); sys.exit(1)
    if abs(uMin-2.5)>1e-12 or abs(uMax-2.5)>1e-12:
        print(f'FAIL: uniform tracer mapping at step {step}: min/max={uMin}/{uMax}'); sys.exit(1)
if 'QUAL_END' not in text:
    print('FAIL: qualification driver did not reach QUAL_END'); sys.exit(1)
print('PASS: cell cycle =', ' -> '.join(map(str,[int(base_match.group(1))]+counts)))
print(f'PASS: volume/uniform/tracer integral relative tolerance <= {a.rel_tol:g}')
