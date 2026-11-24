import os, json, sys
cwd = os.getcwd(); sys.path.insert(0, cwd)
import w5blib as W
from w5blib import QD
from fractions import Fraction as Ff
J = json.load(open('out/d29/candidates.json'))
cands = [dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),29), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in J['cands']]
F29 = W.make_field(29, split=[5], inert=[2])
import random
random.seed(29)
sample = random.sample(cands, 20)
ok = 0
for c in sample:
    xa = QD(abs(c['x'].a),0,29) # not right; need |x| = abs value
    # |x| as QD: if x<0 then -x else x. Determine sign via embedding
    sgn = 1 if c['x'].emb(60) > 0 else -1
    xabs = c['x'] if sgn==1 else -c['x']
    x = QD(1,0,29)
    for g,e in zip(F29['letters'], c['vx']):
        if e: x = x*(g**int(e))
    okx = (x == xabs)
    f1 = W.factor_sunit(1-c['x'], F29)
    ok1 = (f1 == c['v1m'])
    okabs = abs(c['x'].emb(60)) < 1
    if okx and ok1 and okabs: ok += 1
    else: print("  MISMATCH:", c['x'], okx, ok1, okabs)
print(f"d29 factorization re-check (sign-aware): {ok}/20 OK")

# Also count positives/negatives and heights for d29
neg = sum(1 for c in cands if c['x'].emb(60)<0)
from collections import Counter
hc = Counter(c['h'] for c in cands)
print("d29: n=",len(cands), "negatives:", neg, "height dist:", dict(hc))
