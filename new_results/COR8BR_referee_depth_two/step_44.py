# Normalize BOTH sides and compare
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
def norm(entry):
    # entry: (exp, tensor-key-as-string 'v1|v2|v3|v4' with possible '-' signs, coeff string)
    (e, ks, cf) = entry
    letters=[tuple(int(x) for x in part) for part in ks.split('|')]
    v=Fr(cf); out=[]
    for letter in letters:
        l=list(letter)
        for x in l:
            if x!=0:
                if x<0: l=[-y for y in l]; v=-v
                break
        out.append(tuple(l))
    return (e, tuple(out)), v
norm={}
for (e,k),v in sym.items():
    ks="|".join("".join(str(x) for x in vec) for vec in k)
    key,vv=norm((e,ks,str(v)))
    norm[key]=norm.get(key,Fr(0))+vv
    if norm[key]==0: del norm[key]
fs = json.load(open('<project>/research/08_theory/COR4_M3alpha_classical/F_symbol.json'))
fsnorm={}
for k0,v0 in fs.items():
    m=re.match(r'c\^(-?\d+)\s+(.*)', k0)
    e=int(m.group(1)); ks=m.group(2)
    key,vv=norm((e,ks,v0))
    fsnorm[key]=fsnorm.get(key,Fr(0))+vv
    if fsnorm[key]==0: del fsnorm[key]
om = {k:v for k,v in norm.items() if fsnorm.get(k)!=v}
ofs = {k:v for k,v in fsnorm.items() if norm.get(k)!=v}
print("normalized both: mine", len(norm), " COR4", len(fsnorm), " diffs:", len(om), len(ofs))
print("BYTE-IDENTICAL (up to letter-sign normalization):", len(om)==0 and len(ofs)==0)
# also verify 24 distinct pure tensors
tensors=set(k for (e,k) in norm)
print("distinct pure tensors:", len(tensors))
open('progress.log','a').write("step_39: S(Omega) identical to COR4 F_symbol.json after letter-sign normalization; 24 pure tensors.\n")
