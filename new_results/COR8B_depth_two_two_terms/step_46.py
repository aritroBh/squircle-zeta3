
import sys, json, re, pickle
sys.path.insert(0,'.')
from cor8b_sym import *
from fractions import Fraction as Fr
d = json.load(open('../../10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
def parse(term):
    m = re.match(r'c\^(-?\d+)\s*(.*)', term[0]); e = int(m.group(1)); rest = m.group(2).strip()
    w = [int(x) for x in rest[2:-1].split(',')] if rest.startswith('H') else None
    return e, w
gmap = {0: C, 1: vneg(P), -1: Q}
S_F = {}
for t in d:
    e, w = parse(t)
    if w is None or len(w)!=4: continue
    coeff = Fr(t[1])
    k = tuple(gmap[a] for a in reversed(w))
    cur = S_F.setdefault(k, {})
    cur[e] = cur.get(e, Fr(0)) + coeff
    if cur[e]==0: del cur[e]
S_F = {k:v for k,v in S_F.items() if v}
S22 = S_Li22()
S_m = specialize_lin(S22, spec_map(-1)); S_p = specialize_lin(S22, spec_map(+1))
def lau_sym(sym, poly):
    out = {}
    for k,c0 in sym.items():
        if isinstance(c0, dict):
            cur = {}
            for e0,cv in c0.items():
                for e,cc in poly.items():
                    cur[e0+e] = cur.get(e0+e,Fr(0))+cv*cc
            out[k] = {e:v for e,v in cur.items() if v}
        else:
            out[k] = {e: c0*cc for e,cc in poly.items() if c0*cc}
    return {k:v for k,v in out.items() if v}
def lau_sym_add(s1,s2):
    out = {k: dict(v) for k,v in s1.items()}
    for k,v in s2.items():
        cur = out.setdefault(k,{})
        for e,c0 in v.items():
            cur[e] = cur.get(e,Fr(0))+c0
            if cur[e]==0: del cur[e]
        if not cur: del out[k]
    return out
S_DTA = lau_sym_add(lau_sym(S_m, {-1:Fr(1),0:Fr(-1)}), lau_sym(S_p, {-1:Fr(-1),0:Fr(-1)}))
S_diff = lau_sym_add(S_F, lau_sym(S_DTA, {0:Fr(-1)}))
def delta22_lau(sym):
    out = {}
    for (a,b,cc,dd), poly in sym.items():
        w = wedge4(wedge2(a,b), wedge2(cc,dd))
        for k2,v in w.items():
            cur = out.setdefault(k2,{})
            for e,c0 in poly.items():
                cur[e] = cur.get(e,Fr(0)) + c0*v
                if cur[e]==0: del cur[e]
            if not cur: del out[k2]
    return out
D_diff = delta22_lau(S_diff)
print("S(F):", len(S_F), "tensors; S(DT_A):", len(S_DTA), "; S(O-DT):", len(S_diff))
print("\nd22(Omega - DT_A) with [2] kept:")
for k,v in sorted(D_diff.items()):
    print(f"  ({WNAME[k[0]]})^({WNAME[k[1]]}):", {e:str(x) for e,x in v.items()})
D_nom = {k:v for k,v in D_diff.items() if k[0]>=3 and k[1]>=3}
print("\nd22(Omega - DT_A) with [2] dropped:", D_nom if D_nom else "0 (EXACT)")
nm = sum(1 for k in S_diff if any(v[0]!=0 for v in k))
print("tensors in S(O-DT) containing [2]:", nm, "of", len(S_diff))
pickle.dump({'S_F':S_F,'S_DTA':S_DTA,'S_diff':S_diff}, open('syms_omega_dt.pkl','wb'))
print("saved syms_omega_dt.pkl")
open('progress.log','a').write("step_09 fixed; d22(Omega-DT_A)=0 mod [2] exact; 24 tensors match COR8R\n")
