import sys, os, pickle
sys.path.insert(0, os.getcwd())
import mpmath as mp
mp.mp.dps = 200
import sympy as sp
from dod3_R import *

vals200 = {l2: mp.log(2), l3: mp.log(3), lu: mp.log(2-mp.sqrt(3))}
def ev(e):
    e = e.subs(vals200)
    e = e.replace(Li2A, lambda z: mp.re(mp.polylog(2, mp.mpc(sp.N(z,300))))).replace(Li3A, lambda z: mp.re(mp.polylog(3, mp.mpc(sp.N(z,300)))))
    return mp.mpf(sp.N(e, 200))

(i,j,k) = (6,1,5)
xi = Xs[i]; a = pdiv(Xs[j],xi); b = pdiv(Xs[k],xi)
ka = pdiv(a, psub(a,ONE)); kb = pdiv(b, psub(b,ONE))
lo = psub(ONE,xi); da, db = pmul(ka,lo), pmul(kb,lo)
bnd = -logfac_pair(lo)**3/3
pieces = [bnd, Mcf(a,xi), Mcf(b,xi), QnonC(a,b,xi), Swrem(ka,kb), pneg(Swrem(da,db))]
tot = sp.expand(sum(pieces))
R2 = sp.expand(tot * sp_of(pdiv(ONE,xi)))
R1 = build_R(i,j,k)
print("sum pieces*xi^-1:", mp.nstr(ev(R2),16))
print("build_R         :", mp.nstr(ev(R1),16))
print("diff of symbolic expressions:", sp.simplify(R1-R2) != 0)
# check the actual numeric J path
SwpRN = None
print("R2*xi =", mp.nstr(ev(tot),16))
print("R1*xi =", mp.nstr(ev(R1)*pfloat(xi),16))
# hmm compare atoms of R1 vs pieces sum
A1 = R1.atoms(Li2A)|R1.atoms(Li3A); A2 = R2.atoms(Li2A)|R2.atoms(Li3A)
print("atoms equal:", A1==A2, len(A1), len(A2))
