import mpmath as mp, pickle, re
from fractions import Fraction as Fr
mp.mp.dps = 50
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator

# check polylog at negative reals
for z in [Fr(-2), Fr(-80), Fr(-3)]:
    print("Li3(", z, ") =", mp.nstr(mp.polylog(3, M(z)), 20))

# recompute the symbol identity check directly: does bigM*c == t exactly?
import sympy as sp, itertools, sys
sys.path.insert(0,".")
from oct8 import *
sols = pickle.load(open("sols235.pkl","rb"))
vecfun = vec_of_rational
letters = [2,3,5]
Q = Quotient(letters)
d = pickle.load(open("symbol_identities_B.pkl","rb"))
G_names, D_names, sol_out = d["G_names"], d["D_names"], d["sol_out"]

def tensor2(u, v):
    S = {}
    for a, ca in u.items():
        for b, cb in v.items():
            S[(a,b)] = S.get((a,b), Fr(0)) + ca*cb
    return {k: x for k,x in S.items() if x != 0}

# rebuild G and D symbols
G_syms = [S_Li3(Fr(gn[4:-1]), vecfun) for gn in G_names]
def Dsym(dn):
    if dn.startswith("-Li2("):
        m = re.match(r"-Li2\((.*)\)\*log(\d+)", dn)
        w, p = Fr(m.group(1)), int(m.group(2))
        T = tensor2(vecfun(1-w), vecfun(w)); a = p
    else:
        ps = [int(x[3:]) for x in dn.split("*")]
        a, b, c2 = ps
        T = sadd(tensor2({b:Fr(1)},{c2:Fr(1)}), tensor2({c2:Fr(1)},{b:Fr(1)}))
    gensym = {}
    for (b,c_), cf in T.items():
        for t in [(a,b,c_),(b,a,c_),(b,c_,a)]:
            gensym[t] = gensym.get(t, Fr(0)) + cf
    return {k:v for k,v in gensym.items() if v != 0}
D_syms = [Dsym(dn) for dn in D_names]

def S_Li12_xy(x, y):
    L1 = vecfun(x*y-1); L2 = vecfun(x*(y-1)/(x-1)); L3 = vecfun(y)
    L4 = vecfun(y-1);   L5 = vecfun(x-1);           L6 = vecfun(x*y)
    L7 = vecfun(x/(x-1))
    return sadd(tensor3(L1,L2,L3), tensor3(L4,L5,L3), tensor3(L1,L6,L7), tensor3(L4,L3,L5))
targetsB = {"Lam1": (Fr(1,3), Fr(-1,3)), "Lam2": (Fr(3,1), Fr(-1,9)),
            "Lam3": (Fr(2,5), Fr(1,4)),  "Lam4": (Fr(5,2), Fr(1,10))}

for name in ["Lam1","Lam2"]:
    c = [Fr(x) for x in sol_out[name]]
    tot = {}
    for i, S in enumerate(G_syms):
        if c[i]: tot = sadd(tot, smul(S, c[i]))
    for j, S in enumerate(D_syms):
        if c[len(G_syms)+j]: tot = sadd(tot, smul(S, c[len(G_syms)+j]))
    targ = S_Li12_xy(*targetsB[name])
    diff = sadd(targ, smul(tot, -1))
    print(name, "symbol identity residual (should be 0):", diff if diff else "EXACT 0")
