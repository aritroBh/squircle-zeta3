import sys, json, re
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import cor8br_sym as E
M=(1,0,0,0); C=(0,1,0,0); P=(0,0,1,0); Q=(0,0,0,1)
def g(a):
    if a==0: return {C:Fr(1)}
    if a==1: return {P:Fr(-1)}
    if a==-1: return {Q:Fr(1)}
def S_H(w):
    s={():Fr(1)}
    for a in reversed(w):
        s=E.tensor_right(s, g(a))
    return s
d = json.load(open('<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
def parse(term):
    m = re.match(r'c\^(-?\d+)\s*(.*)', term[0]); e=int(m.group(1)); rest=m.group(2).strip()
    w=[int(x) for x in rest[2:-1].split(',')] if rest.startswith('H') else None
    return e,w
sym={}
for t in d:
    e,w=parse(t)
    if w is None or len(w)!=4: continue
    cf=Fr(t[1])
    for k,v in S_H(w).items():
        sym[(e,k)]=sym.get((e,k),Fr(0))+cf*v
        if sym[(e,k)]==0: del sym[(e,k)]
def norm_key(e,k,v):
    # normalize sign of each letter: flip letter to make first nonzero component positive, flipping coeff
    kk=[]; vv=v
    for letter in k:
        l=list(letter); sgn=1
        for x in l:
            if x!=0:
                if x<0: sgn=-1
                break
        if sgn<0: l=[-x for x in l]; vv=-vv
        kk.append(tuple(l))
    return (e,tuple(kk)), vv
norm={}
for (e,k),v in sym.items():
    key,vv = norm_key(e,k,v)
    norm[key]=norm.get(key,Fr(0))+vv
    if norm[key]==0: del norm[key]
def keystr(e,k):
    return f"c^{e} " + "|".join("".join(str(x) for x in v) for v in k)
mine = {keystr(e,k): str(v) for (e,k),v in norm.items()}
fs = json.load(open('<project>/research/08_theory/COR4_M3alpha_classical/F_symbol.json'))
only_mine = {k:v for k,v in mine.items() if fs.get(k)!=v}
only_fs = {k:v for k,v in fs.items() if mine.get(k)!=v}
print("normalized: my terms:", len(mine), " COR4 terms:", len(fs))
print("differing:", len(only_mine), len(only_fs))
for k in list(only_mine)[:6]: print("mine:", k, only_mine[k], " COR4:", fs.get(k))
for k in list(only_fs)[:6]: print("COR4:", k, only_fs[k])
