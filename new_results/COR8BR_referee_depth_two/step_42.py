# Compare my S(Omega) with COR4's F_symbol.json exactly
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

def keystr(e,k):
    def vecstr(v):
        return "".join(str(x) for x in v)
    return f"c^{e} " + "|".join(vecstr(v) for v in k)
mine = {keystr(e,k): str(v) for (e,k),v in sym.items()}
fs = json.load(open('<project>/research/08_theory/COR4_M3alpha_classical/F_symbol.json'))
only_mine = {k:v for k,v in mine.items() if fs.get(k)!=v}
only_fs = {k:v for k,v in fs.items() if mine.get(k)!=v}
print("my terms:", len(mine), " COR4 terms:", len(fs))
print("differing (mine\\COR4):", len(only_mine), " (COR4\\mine):", len(only_fs))
for k in list(only_mine)[:10]: print("mine only:", k, only_mine[k], " COR4:", fs.get(k))
for k in list(only_fs)[:10]: print("COR4 only:", k, only_fs[k])
open('progress.log','a').write("step_38: S(Omega) byte-compared with COR4 F_symbol.json.\n")
