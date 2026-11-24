import os, json, time, sys
cwd = os.getcwd()
sys.path.insert(0, cwd)
import w5blib as W
from w5blib import QD
from fractions import Fraction as Ff

# Write initial summary.md
summary = """# ZQ5c: weight-5 formulas for L(chi_D,5), capped pools: D=29, D=13

## Status board
| D | pool | kernel | LLL dims/precisions | relation |
|---|------|--------|--------------------|----------|
| 29 | out/d29/candidates.json (77, B=8, S={29,2,5}) | not yet | not yet | not yet |
| 13 | out/d13/candidates.json (487, B=8, S={13,2,3}) -> will cap ~150 | not yet (ZQ5b: 263 full) | not yet (ZQ5b: 285@400/600/900 spurious) | not yet |
| 21 | out/d21/candidates.json (1210) | - | - | only if 29 or 13 succeeds |

## Weight-3 JSONs (step 0): pending

## Setup
- libs copied: w5blib.py, w5pipe.py (from ZQ5b), planner_li5.py (from ZQ5, validated Li5)
- d17_final.json copied as format reference
- workdir: ZQ5c_zagier_w5_capped
"""
open('summary.md','w').write(summary)

# Verify d29 candidates: re-check factorizations on 20 candidates (step 1 requirement)
J = json.load(open('out/d29/candidates.json'))
print("d29 names:", J['names'], "split:", J['split'], "inert:", J['inert'], "B:", J['B'])
cands = [dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),29), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in J['cands']]
F29 = W.make_field(29, split=[5], inert=[2])
print("F29 letters:", F29['names'])
import random
random.seed(29)
sample = random.sample(cands, 20)
ok = 0
for c in sample:
    # verify x = product letters^vx and 1-x = product letters^v1m (up to sign)
    x = QD(1,0,29)
    for g,e in zip(F29['letters'], c['vx']):
        if e: x = x*(g**int(e))
    okx = (x == c['x'])
    f1 = W.factor_sunit(1-c['x'], F29)
    ok1 = (f1 == c['v1m'])
    # also |x|<1
    okabs = abs(c['x'].emb(60)) < 1
    if okx and ok1 and okabs: ok += 1
    else:
        print("  MISMATCH:", c['x'], okx, ok1, okabs, f1, c['v1m'])
print(f"d29 factorization re-check: {ok}/20 OK")
with open('progress.log','a') as f:
    f.write(f"step_00b_verify: d29 candidates 20/20 factorizations rechecked OK={ok}; summary.md written\n")
