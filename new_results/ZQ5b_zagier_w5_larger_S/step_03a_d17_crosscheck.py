
import sys, json
from fractions import Fraction as Ff
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
ZQ='../ZQ5_zagier_weight5'
J = json.load(open(f'{ZQ}/out/d17/candidates.json'))
print("their letters:", J['names'], J['disp'])
# my field for D=17 with S={17,2}: letters eps, sqrt17, pi2, pi2'  (same set)
F = W.make_field(17, split=[2], inert=[])
print("my letters:", F['names'])
def their_x(c):
    return QD(Ff(c['x'][0]), Ff(c['x'][1]), 17)
n=len(J['cands']); badx=0; bad1=0
import random
random.seed(1)
idxs = random.sample(range(n), 40)
for i in idxs:
    c=J['cands'][i]
    x=their_x(c)
    fx = W.factor_sunit(x,F) if not (x.a==0 and x.b==0) else None
    f1 = W.factor_sunit(1-x,F)
    if fx != c['vx']: badx+=1; print("vx mismatch", i, c['vx'], fx)
    if f1 != c['v1m']: bad1+=1; print("v1m mismatch", i, c['v1m'], f1)
print(f"checked {len(idxs)} candidates: vx mismatches={badx}, v1m mismatches={bad1}, total cands={n}")
