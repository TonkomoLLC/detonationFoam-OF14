#!/usr/bin/env python3
from __future__ import annotations
import argparse, re
from pathlib import Path

METRIC_RE = re.compile(
    r"SORET_METRIC\s+internalFaces=(?P<faces>\d+)\s+"
    r"maxAbsJH=(?P<h>[-+0-9.eE]+)\s+"
    r"maxAbsJH2=(?P<h2>[-+0-9.eE]+)\s+"
    r"maxAbsSumJ=(?P<sum>[-+0-9.eE]+)\s+"
    r"sumRel=(?P<sumrel>[-+0-9.eE]+)\s+"
    r"publishedRelationRel=(?P<rel>[-+0-9.eE]+)\s+"
    r"inferredDirectH2Rel=(?P<h2direct>[-+0-9.eE]+)"
)

def parse(path: Path):
    text=path.read_text(errors='replace')
    m=METRIC_RE.search(text)
    if not m:
        raise SystemExit(f"FAIL: no SORET_METRIC in {path}")
    vals={k:(int(v) if k=='faces' else float(v)) for k,v in m.groupdict().items()}
    if 'SORET_END' not in text:
        raise SystemExit(f"FAIL: no SORET_END in {path}")
    return text, vals

ap=argparse.ArgumentParser()
ap.add_argument('published')
ap.add_argument('off')
a=ap.parse_args()
ptext,p=parse(Path(a.published)); otext,o=parse(Path(a.off))

assert p['faces']==99, f"expected 99 internal faces, got {p['faces']}"
assert o['faces']==99, f"off case expected 99 internal faces, got {o['faces']}"
assert 'published OF8 H/H2 Soret compatibility enabled' in ptext
assert 'H2 -> TDRatio_H assignment; TDRatio_H2 remains zero' in ptext
assert 'reading published OF8 Soret coefficients from constant/thermoDiff' in ptext
assert p['h'] > 1e-8, f"published H flux not activated: {p['h']:.3e}"
assert p['h2'] > 1e-9, f"published H2 corrected flux not activated: {p['h2']:.3e}"
assert p['sumrel'] < 2e-10, f"sum species flux not closed: {p['sumrel']:.3e}"
assert p['rel'] < 2e-10, f"published H/H2 flux identity failed: {p['rel']:.3e}"
assert p['h2direct'] < 2e-10, f"inferred direct H2 Soret is not zero: {p['h2direct']:.3e}"
assert o['h'] < 1e-12, f"Soret-off H flux unexpectedly nonzero: {o['h']:.3e}"
assert o['h2'] < 1e-12, f"Soret-off H2 flux unexpectedly nonzero: {o['h2']:.3e}"
assert 'published OF8 H/H2 Soret compatibility enabled' not in otext

print('PASS: publishedOF8 mode loaded constant/thermoDiff')
print(f"PASS: H2-only polynomial drives H flux: max|jH|={p['h']:.6e}")
print(f"PASS: H2 direct Soret remains zero by corrected-flux identity: rel={p['h2direct']:.3e}")
print(f"PASS: published H/H2 corrected-flux identity: rel={p['rel']:.3e}")
print(f"PASS: zero-net-mass species-flux closure: rel={p['sumrel']:.3e}")
print(f"PASS: off control remains zero: max|jH|={o['h']:.3e}, max|jH2|={o['h2']:.3e}")
