#!/usr/bin/env python3
from pathlib import Path
import re, sys

log = Path('log.foamRun')
if not log.exists():
    print('FAIL: log.foamRun not found')
    sys.exit(2)
s = log.read_text(errors='replace')

checks = []
def check(ok, msg):
    checks.append(bool(ok))
    print(('PASS: ' if ok else 'FAIL: ') + msg)

check('published OF8 H/H2 Soret compatibility enabled' in s,
      'publishedOF8 Soret mode constructed')
check('automaticUnrefinement=true' in s,
      'planarRefiner automatic unrefinement constructed')

refines = re.findall(r'Planar refined from\s+(\d+)\s+to\s+(\d+)', s)
unrefs = re.findall(r'planarRefiner: unrefined\s+(\d+)', s)
check(bool(refines), 'at least one planar refinement event occurred')
check(bool(unrefs), 'at least one automatic unrefinement pass occurred')

fatal = ('FOAM FATAL' in s) or ('Floating point exception' in s)
check(not fatal, 'no OpenFOAM fatal/FPE in solver log')

ended = ('End\n' in s) or ('Shock-position limit reached' in s and 'ExecutionTime' in s)
check(ended, 'solver reached a normal stop/end condition')

if refines:
    seq = ', '.join(f'{a}->{b}' for a,b in refines[:8])
    print('INFO: refinement events:', seq)
if unrefs:
    print('INFO: unrefinement selected split counts:', ', '.join(x for x in unrefs[:12]))

# Report temperature extrema from the last available time if OpenFOAM prints them.
for pat, label in [
    (r'min/max\(T\)\s*=\s*([^\n]+)', 'T min/max'),
    (r'Leading shock location\s*=\s*([^\n]+)', 'shock position')
]:
    vals = re.findall(pat, s, re.I)
    if vals:
        print(f'INFO: final {label}: {vals[-1]}')

if all(checks):
    print('='*70)
    print('H2/O2 laptop integrated auto-unrefinement + Soret smoke: PASSED')
    print('='*70)
    sys.exit(0)

print('='*70)
print('H2/O2 laptop smoke did not satisfy every acceptance check.')
print('Return log.foamRun; the case is intentionally diagnostic.')
print('='*70)
sys.exit(1)
