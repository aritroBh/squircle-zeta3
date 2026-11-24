# step_97: octagon delta invariance under Landen triple T1 and duplication, exact
from fractions import Fraction as Fr
def wedge(u,v):
    out={}
    for a,ca in u.items():
        for b,cb in v.items():
            if a==b: continue
            if a<b: key,s=(a,b),1
            else: key,s=(b,a),-1
            out[key]=out.get(key,Fr(0))+s*ca*cb
    return {k:c for k,c in out.items() if c}
def dterm(c, vz, vomz):
    out={}
    w=wedge(vz,vomz)
    for (a,b),cab in w.items():
        for l3,c3 in vz.items():
            out[((a,b),l3)] = out.get(((a,b),l3),Fr(0))+c*cab*c3
    return out
def vadd(*vs):
    r={}
    for v in vs:
        for k,c in v.items():
            r[k]=r.get(k,Fr(0))+c
            if r[k]==0: del r[k]
    return r
# letters r2 (sqrt2), u (1+sqrt2): vecs from field table
V = {"d": {"r2":Fr(1),"u":Fr(-1)}, "1-d": {"u":Fr(-1)},
     "k": {"u":Fr(-1)}, "1-k": {"r2":Fr(1),"u":Fr(-1)},
     "mr2h": {"r2":Fr(-1)}, "1-mr2h": {"r2":Fr(-1),"u":Fr(1)},
     "r2h": {"r2":Fr(-1)}, "1-r2h": {"r2":Fr(-1),"u":Fr(-1)},
     "x1": {"u":Fr(-1)}, "1-x1": {"r2":Fr(1)},
     "x2": {"r2":Fr(-1),"u":Fr(-1)}, "1-x2": {"r2":Fr(-1)}}
# Landen triples: T1 = {d, k, mr2h} (= {2-sqrt2, sqrt2-1, -sqrt2/2}); T2 = {r2h, x2, x1}
T1 = vadd(dterm(1,V["d"],V["1-d"]), dterm(1,V["k"],V["1-k"]), dterm(1,V["mr2h"],V["1-mr2h"]))
T2 = vadd(dterm(1,V["r2h"],V["1-r2h"]), dterm(1,V["x2"],V["1-x2"]), dterm(1,V["x1"],V["1-x1"]))
DUP = vadd(dterm(1,V["r2h"],V["1-r2h"]), dterm(1,V["mr2h"],V["1-mr2h"]))
print("delta(Landen T1 {d,k,mr2h}) =", T1 if T1 else 0)
print("delta(Landen T2 {r2h,x2,x1}) =", T2 if T2 else 0)
print("delta({r2h}+{mr2h}) =", DUP if DUP else 0, " [duplication partner, = 7/32 zeta(3) class]")
# so the octagon delta computed from the 6 coefficients is invariant under the whole relation lattice
