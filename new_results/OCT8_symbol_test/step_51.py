import mpmath as mp, pickle, sympy as sp, itertools, re, sys
from fractions import Fraction as Fr
sys.path.insert(0,".")
from oct8 import *
mp.mp.dps = 300
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator
sols = pickle.load(open("sols235.pkl","rb"))
vecfun = vec_of_rational; letters = [2,3,5]
Q = Quotient(letters)
dB = pickle.load(open("symbol_identities_B.pkl","rb"))
G_names, D_names = dB["G_names"], dB["D_names"]
G_syms = [S_Li3(Fr(gn[4:-1]), vecfun) for gn in G_names]
def tensor2(u, v):
    S = {}
    for a, ca in u.items():
        for b, cb in v.items():
            S[(a,b)] = S.get((a,b), Fr(0)) + ca*cb
    return {k: x for k,x in S.items() if x != 0}
def Dsym(dn):
    m = re.match(r"-Li2\((.*)\)\*log(\d+)", dn)
    w, p = Fr(m.group(1)), int(m.group(2))
    T = tensor2(vecfun(1-w), vecfun(w)); a = p
    r = {}
    for (b,c_), cf in T.items():
        for t in [(a,b,c_),(b,a,c_),(b,c_,a)]:
            r[t] = r.get(t, Fr(0)) + cf
    return {k:v for k,v in r.items() if v != 0}
D_syms = [Dsym(dn) for dn in D_names]
gens = [Q.coord(S) for S in G_syms]; Dn = [Q.coord(S) for S in D_syms]
bigM = sp.Matrix.hstack(*[g.T for g in gens] + [d.T for d in Dn])
print("bigM:", bigM.rows, "x", bigM.cols, " rank", bigM.rank())
NS = bigM.nullspace()
print("nullity:", len(NS))
# value map
logv = {p: mp.log(p) for p in (2,3,5)}
Li3v = {gn: mp.polylog(3, M(Fr(gn[4:-1]))) for gn in G_names}
def Dval(dn):
    m = re.match(r"-Li2\((.*)\)\*log(\d+)", dn)
    w, p = Fr(m.group(1)), int(m.group(2))
    return -mp.polylog(2, M(w))*logv[p]
Dv = {dn: Dval(dn) for dn in D_names}
vals = []
for nvec in NS:
    tot = mp.mpf(0)
    for i, gn in enumerate(G_names):
        if nvec[i]: tot += M(sp.Rational(nvec[i]))*Li3v[gn]
    for j, dn in enumerate(D_names):
        if nvec[len(G_names)+j]: tot += M(sp.Rational(nvec[len(G_names)+j]))*Dv[dn]
    vals.append(tot)
# analyze the value space: PSLQ each value against {z3, z2log p} to see if all land there
z3, z2 = mp.zeta(3), mp.zeta(2)
inb = [z3, z2*logv[2], z2*logv[3], z2*logv[5]]
outside = []
for k, v in enumerate(vals):
    rel = mp.pslq([v]+inb, tol=mp.mpf('1e-270'), maxcoeff=10**40, maxsteps=100000)
    if rel is None:
        outside.append((k, v))
print("null vectors with value OUTSIDE {z3, z2log p}:", len(outside))
for k, v in outside[:10]:
    print("  idx", k, " value", mp.nstr(v, 15))
pickle.dump({"NS":[[str(x) for x in nv] for nv in NS]}, open("nullspace.pkl","wb"))
