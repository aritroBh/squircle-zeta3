import os, json
from fractions import Fraction as Fr
cwd = "<project>/research/08_theory/G2_M4alpha"
out = open(cwd + "/flat_exact_N640_K4.out").read()
print(out[-200:])
if os.path.exists(cwd + "/flat_exact_N640_K4.json"):
    ref = json.load(open(cwd + "/flat_exact_N640_K4.json"))['4']
    frm = json.load(open(cwd + "/flat_exact_N640_K4_formula.json"))['4']
    ok = True
    for nv in range(521, 641):
        for k in ['(0, 0, 0)', '(1, 0, 0)', '(0, 1, 0)', '(0, 0, 1)', '(2, 0, 0)']:
            if Fr(ref[str(nv)].get(k, '0')) != Fr(frm[str(nv)].get(k, '0')):
                ok = False; print("MISMATCH", nv, k)
    print("formula table == flat_exact reference also for 521 <= n <= 640:", ok)
