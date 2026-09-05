#!/usr/bin/env python3
"""Utilities for detonationFoam/reactingDNS constant/binaryDiff dictionaries."""
from __future__ import annotations
from pathlib import Path
import re

SPECIES_RE = re.compile(r"species\s+\d+\s*\((.*?)\)\s*;", re.S)


def strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    text = re.sub(r"//.*?$", "", text, flags=re.M)
    return text


def read_species(path: Path) -> list[str]:
    text = strip_comments(path.read_text())
    m = SPECIES_RE.search(text)
    if not m:
        raise ValueError(f"cannot parse species list from {path}")
    species = re.findall(r'"[^"]+"|[^\s()]+', m.group(1))
    return [s.strip('"') for s in species if s.strip()]


def _find_matching(text: str, start: int) -> int:
    depth = 0
    quote = None
    i = start
    while i < len(text):
        c = text[i]
        if quote:
            if c == quote and (i == 0 or text[i-1] != "\\"):
                quote = None
        elif c in "'\"":
            quote = c
        elif c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise ValueError("unmatched '{' in OpenFOAM dictionary")


def top_level_blocks(text: str):
    """Yield (key, body) for top-level dictionary blocks."""
    text = strip_comments(text)
    i = 0
    n = len(text)
    token = re.compile(r'"[^"]+"|[^\s{};]+')
    while i < n:
        m = token.search(text, i)
        if not m:
            break
        key = m.group(0).strip('"')
        j = m.end()
        while j < n and text[j].isspace():
            j += 1
        if j < n and text[j] == '{':
            end = _find_matching(text, j)
            yield key, text[j+1:end]
            i = end + 1
        else:
            semi = text.find(';', j)
            i = (semi + 1) if semi >= 0 else j + 1


def parse_legacy_binary_diff(path: Path) -> dict[str, tuple[float,float,float,float]]:
    coeffs = {}
    for key, body in top_level_blocks(path.read_text()):
        vals = []
        ok = True
        for name in ('Diff1','Diff2','Diff3','Diff4'):
            m = re.search(rf'\b{name}\s+([^;]+);', body)
            if not m:
                ok = False
                break
            vals.append(float(m.group(1).strip()))
        if ok:
            coeffs[key] = tuple(vals)  # type: ignore[assignment]
    if not coeffs:
        raise ValueError(f"no Diff1..Diff4 pair blocks found in {path}")
    return coeffs


def expected_pairs(species: list[str]) -> list[tuple[str,str]]:
    return [(a,b) for i,a in enumerate(species) for b in species[i+1:]]


def resolve_pair(coeffs: dict[str, tuple[float,float,float,float]], a: str, b: str):
    if f'{a}-{b}' in coeffs:
        return coeffs[f'{a}-{b}']
    if f'{b}-{a}' in coeffs:
        return coeffs[f'{b}-{a}']
    raise KeyError(f"missing binary diffusion pair {a}-{b} (or reversed)")


def render_thermophysical_transport(species: list[str], coeffs: dict[str, tuple[float,float,float,float]]) -> str:
    lines = [
        '/*--------------------------------*- C++ -*----------------------------------*\\',
        '| =========                 |                                                 |',
        '| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |',
        '|  \\    /   O peration     | Version:  14                                    |',
        '|   \\  /    A nd           |                                                 |',
        '|    \\/     M anipulation  |                                                 |',
        '\\*---------------------------------------------------------------------------*/',
        'FoamFile','{','    format      ascii;','    class       dictionary;',
        '    location    "constant";','    object      thermophysicalTransport;','}','',
        '// Converted from legacy constant/binaryDiff (Diff1..Diff4).',
        'laminar','{','    model legacyMixtureAverageFourier;',
        '    legacyThermalDiffusionMode off;','', '    D','    {'
    ]
    for a,b in expected_pairs(species):
        c = resolve_pair(coeffs,a,b)
        lines += [
            f'        {a}-{b}', '        {',
            '            type legacyBinaryDiffusionCoefficient;',
            '            coeffs (' + ' '.join(f'{x:.17g}' for x in c) + ');',
            '        }'
        ]
    lines += ['    }','}','', '// ************************************************************************* //','']
    return '\n'.join(lines)


def render_legacy_binary_diff(species: list[str], pair_coeffs) -> str:
    lines = [
        '/*--------------------------------*- C++ -*----------------------------------*\\',
        '| =========                 |                                                 |',
        '| \\      /  F ield         | OpenFOAM legacy transport table                 |',
        '\\*---------------------------------------------------------------------------*/',
        'FoamFile','{','    format ascii;','    class dictionary;',
        '    location "constant";','    object binaryDiff;','}','',
        '// Synthetic Stage-D3 legacy Diff1..Diff4 data', ''
    ]
    for i,(a,b) in enumerate(expected_pairs(species)):
        c = pair_coeffs(i,a,b)
        lines += [f'{a}-{b}','{',
                  f'    Diff1 {c[0]:.17g};',f'    Diff2 {c[1]:.17g};',
                  f'    Diff3 {c[2]:.17g};',f'    Diff4 {c[3]:.17g};','}','']
    return '\n'.join(lines)
