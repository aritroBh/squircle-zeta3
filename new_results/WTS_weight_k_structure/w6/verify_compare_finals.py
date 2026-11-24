"""verify step: exact comparison of the SLOW sympy reducer output (verify_slow_reduce.py) with the fast reducer
finals_K.json: every key, every pole, every Q(sqrt2) coefficient must agree as exact rationals.
usage: verify_compare_finals.py SLOWDIR K [K ...]
"""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, W)
from fractions import Fraction as Fr
import verify_lib as VL

SLOW = sys.argv[1]
for k in [int(x) for x in sys.argv[2:]]:
    fast = VL.load_finals(os.path.join(W, f"finals_{k}.json"))
    slow = json.load(open(os.path.join(SLOW, f"slow_finals_{k}.json")))["slow"]
    ok = True
    keys_fast = {kk for kk, d in fast.items() if d}
    keys_slow = {tuple(map(int, s.split(","))) for s, d in slow.items() if any(not l.startswith("_") for l in d)}
    if keys_fast != keys_slow:
        ok = False
        print(f"k={k}: KEY SETS DIFFER fast-only {sorted(keys_fast-keys_slow)} slow-only {sorted(keys_slow-keys_fast)}")
    nchk = 0
    for key in sorted(keys_fast & keys_slow):
        sd = slow[f"{key[0]},{key[1]}"]
        assert sd["_leftover_is_zero"], (k, key)
        spoles = {l: (Fr(v[0]), Fr(v[1])) for l, v in sd.items() if not l.startswith("_")}
        fpoles = fast[key]
        if set(spoles) != set(fpoles):
            ok = False
            print(f"k={k} key {key}: pole sets differ fast={sorted(fpoles)} slow={sorted(spoles)}")
        for l in set(spoles) & set(fpoles):
            nchk += 1
            if spoles[l] != fpoles[l]:
                ok = False
                print(f"k={k} key {key} pole {l}: fast {fpoles[l]} slow {spoles[l]}")
    print(f"k={k}: slow == fast exactly for all {len(keys_fast & keys_slow)} keys / {nchk} residues: {ok}")
