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

def normalize(symdict):
    out={}
    for (e,k),v in symdict.items():
        kk=[]; vv=v
        for letter in k:
            l=list(letter)
            for x in l:
                if x!=0:
                    if x<0: l=[-y for y in l]; vv=-vv
                    break
            kk.append(tuple(l))
        key=(e,tuple(kk))
        out[key]=out.get(key,Fr(0))+vv
    return {k:v for k,v in out.items() if v}
mine_n = normalize(sym)
fs = json.load(open('<project>/research/08_theory/COR4_M3alpha_classical/F_symbol.json'))
fs_raw={}
for k0,v0 in fs.items():
    m=re.match(r'c\^(-?\d+)\s+(.*)', k0)
    e=int(m.group(1))
    letters=tuple(tuple(int(x) for x in re.findall(r'-?\d+', part)) for part in m.group(2).split('|'))
    fs_raw[(e,letters)]=Fr(v0)
fs_n = normalize(fs_raw)
diffs1 = {k:(mine_n[k],fs_n.get(k)) for k in mine_n if fs_n.get(k)!=mine_n[k]}
diffs2 = {k:(fs_n[k],mine_n.get(k)) for k in fs_n if mine_n.get(k)!=fs_n[k]}
print("mine", len(mine_n), "fs", len(fs_n), "diffs:", len(diffs1), len(diffs2))
print("IDENTICAL:", len(diffs1)==0 and len(diffs2)==0)
for k,v in list(diffs1.items())[:5]: print("mine:",k,v)
for k,v in list(diffs2.items())[:5]: print("fs:",k,v)
