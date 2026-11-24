import sys, os, pickle, sympy as sp
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, wd)
from esum_engine import esum, Hf, H2f, k, x, L1, L2, L3, AA, z2
O5 = '<project>/research/08_theory/OCT6_octagon_closed/oct5'
Jloc = pickle.load(open(O5+'/Jloc.pkl','rb'))
print("Jloc keys:", len(Jloc), sorted(Jloc.keys())[:8])
# rebuild: Jloc(n,a,b) = moment at l=n+1
z3 = sp.Symbol('z3')
def Hh(l): return sp.harmonic(l)
def H2h(l): return sp.harmonic(l,2)
def JLOC(n,a,b):
    l = n+1
    if (a,b)==(0,0): return 1/sp.Integer(l)
    if (a,b)==(1,0): return -1/sp.Integer(l)**2
    if (a,b)==(2,0): return 2/sp.Integer(l)**3
    if (a,b)==(3,0): return -6/sp.Integer(l)**4
    if (a,b)==(0,1): return -Hh(l)/l
    if (a,b)==(1,1): return Hh(l)/l**2 - sp.pi**2/(6*l) + H2h(l)/l
    if (a,b)==(2,1): return sp.expand(sp.nsimplify(0)*l) # placeholder, derive below
    if (a,b)==(0,2): return (Hh(l)**2 + H2h(l))/l
    raise ValueError((n,a,b))
# need (2,1),(1,2),(0,3),(3,0) etc: derive via series
# int t^(l-1) (log t)^a (log(1-t))^b dt = d^a/dalpha^a d^b/dbeta^b B(l,1) -- use exact expansion
# (2,1): int t^(l-1) log^2 t log(1-t) = -sum_j 2/(j*(l+j)^3)?? expand log(1-t) = -sum t^j/j:
# int = -sum_{j>=1} 1/j * int t^(l+j-1) log^2 t = -sum 2/(j (l+j)^3)  (exact sum: harmonic/zeta combination)
mism = 0
for key in sorted(Jloc.keys()):
    n,a,b = key
    if (a,b) in [(0,0),(1,0),(2,0),(3,0),(0,1),(1,1),(0,2)]:
        mine = JLOC(n,a,b)
        mine = mine.subs(sp.pi**2, 6*sp.Symbol('z2x')) 
        want = Jloc[key].subs(sp.Symbol('z2'), sp.Symbol('z2x')) if Jloc[key].has(sp.Symbol('z2')) else Jloc[key]
        d = sp.simplify(mine - want)
        if d != 0:
            mism += 1
            print("MISMATCH", key, "mine:", mine, "want:", want)
print("first-pass mismatches:", mism, "of", len(Jloc))
# what (a,b) types are in Jloc?
print(sorted(set((a,b) for (n,a,b) in Jloc)))
