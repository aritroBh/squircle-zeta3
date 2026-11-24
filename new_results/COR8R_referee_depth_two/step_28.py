import json, re
import cor8r_sym as S
import sympy as sp
from fractions import Fraction
cc = S.cc
NAMES = ['-1','2','c','1-c','1+c']

# ---- build Omega symbol from the 56 weight-4 HPL terms ----
base = '<project>'
data = json.load(open(base+'/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
omega = {}
nterms = 0
for name, coeff in data:
    m = re.match(r'c\^(-?\d+)\s*H\[([^\]]*)\]', name)
    E = int(m.group(1))
    w = tuple(int(x) for x in m.group(2).split(','))
    if len(w)!=4: continue
    nterms += 1
    hs = S.hpl_sym(w)
    hs = S.S_scaleLP(hs, S.LP_scale(Fraction(coeff), S.LP_pow_c(E)))
    omega = S.S_add(omega, hs)
print("weight-4 terms used:", nterms, " symbol terms:", len(omega))
dOm = S.delta22(omega)
print("\ndelta_{2,2}(Omega) components ((pairI,pairJ) -> Laurent coeff):")
def show_pair(I):
    i,j = S.PAIRS[I]
    return f"[{NAMES[i]}]^[{NAMES[j]}]"
for k,v in sorted(dOm.items()):
    print("  ", (show_pair(k[0]), show_pair(k[1])), " -> ", v)

# ---- claimed formula ----
# (1-1/c)([c]^[1-c])^([1-c]^[1+c]) + (-1-1/c)([c]^[1+c])^([1-c]^[1+c])
Lc, Lm, Lp = (0,0,1,0,0),(0,0,0,1,0),(0,0,0,0,1)
P_cm = S.wedge2(Lc,Lm)  # [c]^[1-c]
P_cp = S.wedge2(Lc,Lp)
P_mp = S.wedge2(Lm,Lp)
def mk_key(P,Q):
    # returns (I,J) with I<J and sign
    (i1,), (i2,) = P.keys(), Q.keys()  # each is a single basis pair
    I,J = i1,i2
    return (I,J) if I<J else (J,I)
k1 = (P_cm.keys().__iter__().__next__(), P_mp.keys().__iter__().__next__())
k2 = (P_cp.keys().__iter__().__next__(), P_mp.keys().__iter__().__next__())
print("\nclaimed component keys:", k1, k2)
claim = {}
A = S.LP_add(S.ONE, S.LP_neg(S.LP_pow_c(-1)))     # 1 - 1/c
B = S.LP_neg(S.LP_add(S.ONE, S.LP_pow_c(-1)))     # -1 - 1/c
claim[k1 if k1[0]<k1[1] else (k1[1],k1[0])] = A
claim[k2 if k2[0]<k2[1] else (k2[1],k2[0])] = B
print("\nclaim == delta22(Omega)? ", claim == dOm)
