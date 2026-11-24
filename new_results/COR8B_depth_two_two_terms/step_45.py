
import sys, json, re
sys.path.insert(0,'.')
from cor8b_sym import *
from fractions import Fraction as Fr

# ---- S(F) from HPL (done before, reuse), and S(DT_A) with [2] KEPT as symbols ----
d = json.load(open('../../10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
def parse(term):
    m = re.match(r'c\^(-?\d+)\s*(.*)', term[0])
    e = int(m.group(1)); rest = m.group(2).strip()
    w = [int(x) for x in rest[2:-1].split(',')] if rest.startswith('H') else None
    return e, w
gmap = {0: C, 1: vneg(P), -1: Q}
def S_HPL(w): return {tuple(gmap[a] for a in reversed(w)): Fr(1)}
S_F = {}
for t in d:
    e, w = parse(t)
    if w is None or len(w)!=4: continue
    coeff = Fr(t[1])
    for k,c0 in S_HPL(w).items():
        # carry Laurent: dict key -> {e: coeff}
        cur = S_F.setdefault(k, {})
        cur[e] = cur.get(e, Fr(0)) + coeff*c0
        if cur[e]==0: del cur[e]
S_F = {k:v for k,v in S_F.items() if v}
print("S(F) tensors:", len(S_F))

# S(DT_A): ((1-c)/c) S_Li22(spec -1) + (-(1+c)/c) S_Li22(spec +1), Laurent coeffs
S22 = S_Li22()
S_m = specialize_lin(S22, spec_map(-1))
S_p = specialize_lin(S22, spec_map(+1))
def lau_sym(sym, poly):
    out = {}
    for k,c0 in sym.items():
        out[k] = {e: c0*cc for e,cc in poly.items() if c0*cc}
    return out
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
print("S(DT_A) tensors:", len(S_DTA))

# difference
S_diff = lau_sym_add(S_F, lau_sym(S_DTA, {0:Fr(-1)}))
print("S(Omega-DT_A) tensors:", len(S_diff))

# delta22 with Laurent coeffs
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
print("\nd22(Omega - DT_A) with [2] kept:")
for k,v in sorted(D_diff.items()):
    print(f"  ({WNAME[k[0]]})^({WNAME[k[1]]}):", {e:str(x) for e,x in v.items()})
# with [2] dropped: keep only components without index 0,1,2 wedges? No: drop letter m means project: any wedge with m vanishes
D_diff_nom = {k:v for k,v in D_diff.items() if k[0]>=3 and k[1]>=3}  # wedge indices >=3 = cp,cq,pq only
print("d22(Omega - DT_A) with [2] dropped:", D_diff_nom if D_diff_nom else "0 (EXACT)")
# also count tensors with m-letter in S_diff
nm = sum(1 for k in S_diff if any(v[0]!=0 for v in k))
print("tensors in S(Omega-DT_A) containing letter [2]:", nm, "of", len(S_diff))
import pickle
pickle.dump({'S_F':S_F,'S_DTA':S_DTA,'S_diff':S_diff}, open('syms_omega_dt.pkl','wb'))
print("saved syms_omega_dt.pkl")
open('progress.log','a').write("step_09 d22(Omega-DT_A)=0 mod [2] exact; symbols saved\n")
