import sys, os, pickle, time
base = "<project>/research/08_theory/"
sys.path.insert(0, base+"DOD5_dodecagon_termwise/")
import mpmath as mp
mp.mp.dps = 620
import sympy as sp
t0=time.time()
sd = pickle.load(open(base+"DOD1_dodecagon/S_atom_decomp.pkl","rb"))
coefs = sd['coefs']
print("60 atom coef keys families:", __import__('collections').Counter(k[0] for k in coefs))
r3 = sp.sqrt(3)
def pq_of(e):
    """exact e in Q(sqrt3) -> (p,q) Fractions with e = p + q*sqrt3, via rationalization"""
    e2 = sp.radsimp(sp.expand(e))
    p = sp.simplify((e2 + e2.subs(r3,-r3))/2)
    q = sp.simplify((e2 - e2.subs(r3,-r3))/(2*r3))
    return sp.Rational(p), sp.Rational(q)
out = {}
for k in coefs:
    p,q = pq_of(coefs[k])
    out[k] = (p,q)
print("rationalized, t=%.0fs"%(time.time()-t0))
# verify Lam/T3a match atom_coefs_pq
ac = pickle.load(open(base+"DOD2_dodecagon_finish/atom_coefs_pq.pkl","rb"))
bad = 0
for k in ac:
    p,q = ac[k]
    p2,q2 = out[k]
    if sp.Rational(p)!=p2 or sp.Rational(q)!=q2:
        bad += 1; print("MISMATCH", k, (p,q), (p2,q2))
print("Lam/T3a coefs match atom_coefs_pq:", bad==0)
# show the depth-1 atom coefficients
fams = {}
for k in out:
    fams.setdefault(k[0], {})[k[1]] = out[k]
for fam in sorted(fams):
    print("\n==", fam, len(fams[fam]))
    for kk in sorted(fams[fam]):
        p,q = fams[fam][kk]
        print("   %-10s (%s) + (%s) sqrt3"%(kk, p, q))
pickle.dump(out, open('dod5_atom_coefs_all.pkl','wb'))
print("saved dod5_atom_coefs_all.pkl, t=%.0fs"%(time.time()-t0))
