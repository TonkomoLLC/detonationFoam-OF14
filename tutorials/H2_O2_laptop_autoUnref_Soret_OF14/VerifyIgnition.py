#!/usr/bin/env python3
from pathlib import Path
import re
import sys

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


def internal_values(path):
    text = Path(path).read_text(errors="replace")
    m = re.search(r"\binternalField\b(.*?)(?=\bboundaryField\b)", text, re.S)
    if not m:
        raise RuntimeError(f"{path}: cannot locate internalField")
    chunk = m.group(1).strip()

    # setFields may write either a uniform value or a nonuniform List<scalar>.
    mu = re.match(rf"uniform\s+({NUMBER})\s*;", chunk)
    if mu:
        return [float(mu.group(1))]

    mn = re.search(r"nonuniform\s+List<scalar>\s+\d+\s*\((.*?)\)\s*;", chunk, re.S)
    if mn:
        vals = [float(x) for x in re.findall(NUMBER, mn.group(1))]
        if vals:
            return vals

    raise RuntimeError(f"{path}: unsupported internalField representation")


try:
    T = internal_values("0/T")
    p = internal_values("0/p")
except Exception as e:
    print(f"FAIL: post-setFields verification: {e}")
    sys.exit(2)

Tmin, Tmax = min(T), max(T)
pmin, pmax = min(p), max(p)

print(f"INFO: post-setFields T range = {Tmin:g} .. {Tmax:g} K")
print(f"INFO: post-setFields p range = {pmin:g} .. {pmax:g} Pa")

ok = True
if Tmax < 1900:
    print("FAIL: ignition temperature was not written to 0/T")
    ok = False
else:
    print("PASS: ignition temperature present")

if pmax < 1.9e6:
    print("FAIL: ignition pressure was not written to 0/p")
    ok = False
else:
    print("PASS: ignition pressure present")

if Tmin > 400 or pmin > 2e5:
    print("FAIL: cold/unburned background state was not retained")
    ok = False
else:
    print("PASS: cold/unburned background retained")

sys.exit(0 if ok else 1)
