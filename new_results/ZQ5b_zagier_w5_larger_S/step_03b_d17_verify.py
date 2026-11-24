
import sys, json, itertools, time
from fractions import Fraction as Ff
from collections import defaultdict
import mpmath as mp
sys.path.insert(0,'.')
import w5blib as W
from w5blib import QD
ZQ='../ZQ5_zagier_weight5'
D=17
J = json.load(open(f'{ZQ}/out/d17/candidates.json'))
R = json.load(open(f'{ZQ}/planner/d17_relation_planner.json'))
cands=[dict(x=QD(Ff(c['x'][0]),Ff(c['x'][1]),17), vx=c['vx'], v1m=c['v1m']) for c in J['cands']]
names=J['names']; d=len(names)
T_coef=Ff(R['T_coef']); z5_coef=Ff(R['z5_coef'])
# map li5 args (strings "a+b*sqrt17") to coefficients
li5map = {k: Ff(v) for k,v in R['li5'].items()}
def xkey(x): return f"{x.a}+{x.b}*sqrt17"
prods = {k: Ff(v) for k,v in R['products'].items()}
# match li5 args to candidates
cvec = {}   # candidate index -> coeff
missing=[]
for ks, cf in li5map.items():
    found=None
    for i,c in enumerate(cands):
        if xkey(c['x'])==ks: found=i; break
    if found is None: missing.append(ks)
    else: cvec[found]=cf
print("li5 terms:", len(li5map), "matched:", len(cvec), "missing:", missing)

# (b) exact full symbol (Li5 + log5 products; pi2/z5/T zero symbol)
T = defaultdict(Ff)
for i, cf in cvec.items():
    s = W.symbol_liw_dense(cands[i]['vx'], cands[i]['v1m'], d, 5)
    for idx, v in enumerate(s):
        if v: T[idx]+= cf*v
idx_of = {n:k for k,n in enumerate(names)}
for nm, cf in prods.items():
    toks = nm.split('*')
    if 'pi2' in toks: continue
    ms = [idx_of[t] for t in toks]
    assert len(ms)==5
    v = W.prod_symbol_dense(ms, d, 5)
    for i2, vv in enumerate(v):
        if vv: T[i2]+= cf*vv
Tz={k:v for k,v in T.items() if v}
print("(b) exact full symbol:", "ZERO" if not Tz else f"NONZERO ({len(Tz)} entries)")
# (c) cobracket
dic = {W.idx_to_tup(k,d,5): v for k,v in T.items()}
for mode in ("full","cut"):
    r=W.cobracket32(dic,d,mode)
    print(f"(c) delta_{{3,2}} total [{mode}]: {len(r)} nonzero")
# Li5 part alone
Tl = defaultdict(Ff)
for i, cf in cvec.items():
    s = W.symbol_liw_dense(cands[i]['vx'], cands[i]['v1m'], d, 5)
    for idx, v in enumerate(s):
        if v: Tl[idx]+= cf*v
dicl = {W.idx_to_tup(k,d,5): v for k,v in Tl.items()}
for mode in ("full","cut"):
    r=W.cobracket32(dicl,d,mode)
    print(f"(c) delta_{{3,2}} Li5-part [{mode}]: {len(r)} nonzero")
# (e)
used = [(k,v) for k,v in list(li5map.items())+list(prods.items()) if v]
dens=[v.denominator for _,v in used]; nums=[abs(v.numerator) for _,v in used]
import math, sympy as sp
print("(e) terms:", len(used)+1, "(+z5)", " max|num|=", max(nums), " lcm den=", math.lcm(*dens), "=", sp.factorint(math.lcm(*dens)))
print("    T_coef:", T_coef, "= 17^5:", abs(T_coef)==17**5)
# (a) fresh value at 1300 digits
P=1300
mp.mp.dps=P+30
t0=time.time()
s=mp.sqrt(17)
Lv=[4+s, s, (3+s)/2, abs((3-s)/2)]
lg=[mp.log(l) for l in Lv]
Lchi=W.Lchi_hurwitz(17,5,P)
print("L(chi17,5) =", mp.nstr(Lchi,30))
tot = T_coef*Lchi/s + z5_coef*mp.zeta(5)
li5cache={}
for i, cf in cvec.items():
    x=cands[i]['x'].emb(P+30)
    v=W.li5(x,P+10)
    tot += cf*v
for nm, cf in prods.items():
    v=mp.mpf(cf.numerator)/cf.denominator
    for t in nm.split('*'): v*= (mp.pi**2 if t=='pi2' else lg[idx_of[t]])
    tot += v
print(f"(a) fresh residual @1300d: {mp.nstr(tot,4)}  ({time.time()-t0:.0f}s)")
