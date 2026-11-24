"""Lift words in the t sections (A3h4_frame n rowidx word) and frame234 sections: admissible?"""
from common import *
bad = 0; seen = set()
for name in TARGETS:
    for o in read_jsonl_xz(f'{STORE}/{name}.core.jsonl.xz'):
        if o['section'] in ('t', 'frame'):
            for t in o['terms']:
                if t[0] == 'A3h4_frame': w = tuple(t[3]); n = t[1]; k = 4 - n
                elif t[0] == 'frame234': w = tuple(t[4]); n = t[2]; k = 5 - n
                else: continue
                seen.add((t[0], n, w))
                if len(w) != k or (w and (w[0] == 0 or w[-1] == -1)): bad += 1
log(f'lift words in t and frame sections of all 74 core files: {len(seen)} distinct (kind, n, word); non-admissible or wrong length: {bad}', 'check_twords.log')
