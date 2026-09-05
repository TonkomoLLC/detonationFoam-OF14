#!/usr/bin/env python3
from pathlib import Path
import argparse, json
from LegacyBinaryDiffTools import read_species, parse_legacy_binary_diff, expected_pairs, resolve_pair, render_thermophysical_transport

def main():
    ap=argparse.ArgumentParser(description='Convert detonationFoam/reactingDNS constant/binaryDiff to OF14 thermophysicalTransport D entries.')
    ap.add_argument('--binary-diff',required=True)
    ap.add_argument('--species-file',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--summary-json')
    a=ap.parse_args()
    sp=read_species(Path(a.species_file)); coeff=parse_legacy_binary_diff(Path(a.binary_diff))
    missing=[]
    for x,y in expected_pairs(sp):
        try: resolve_pair(coeff,x,y)
        except KeyError: missing.append(f'{x}-{y}')
    if missing:
        raise SystemExit('missing binaryDiff pairs: '+', '.join(missing[:20]) + (' ...' if len(missing)>20 else ''))
    Path(a.output).write_text(render_thermophysical_transport(sp,coeff))
    info={'species':len(sp),'requiredPairs':len(expected_pairs(sp)),'inputPairBlocks':len(coeff),'output':str(Path(a.output))}
    if a.summary_json: Path(a.summary_json).write_text(json.dumps(info,indent=2,sort_keys=True)+'\n')
    print(f"Converted legacy binaryDiff: species={info['species']}, required-pairs={info['requiredPairs']}, input-pair-blocks={info['inputPairBlocks']}")
    print(f"Output: {a.output}")
if __name__=='__main__': main()
