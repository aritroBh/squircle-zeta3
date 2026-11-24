import sys, os, pickle, itertools, time
sys.path.insert(0, os.getcwd())
from fractions import Fraction as Fr
from dod3_symbol import *

t0 = time.time()
def enum_sunits(B):
    """all z = +-P2^a P3^b EPS^c, |a|,|b|,|c|<=B, with 1-z also S-unit. Return dict z->(fac_z, fac_w)."""
    # precompute powers
    pw2 = {a: ppow(P2,a) for a in range(-B,B+1)}
    pw3 = {b: ppow(P3,b) for b in range(-B,B+1)}
    pwe = {c: ppow(EPS,c) for c in range(-B,B+1)}
    sols = {}
    for a in range(-B,B+1):
        for b in range(-B,B+1):
            ab = pmul(pw2[a], pw3[b])
            for c in range(-B,B+1):
                base = pmul(ab, pwe[c])
                for sgn in (1,-1):
                    z = base if sgn==1 else pneg(base)
                    if peq(z, ONE) or peq(z, P(0,0)): continue
                    w = psub(ONE, z)
                    fz = factor_sunit(z); fw = factor_sunit(w)
                    if fz and fw and fz[0]!='NOT_SUNIT' and fw[0]!='NOT_SUNIT':
                        sols[z] = (fz, fw)
    return sols

counts = {}
for B in (4,6,8):
    s = enum_sunits(B)
    counts[B] = len(s)
    print("B=%2d: %d solutions   (t=%.1fs)"%(B, len(s), time.time()-t0))
    pickle.dump(s, open('solsF_B%d.pkl'%B,'wb'))

sols = pickle.load(open('solsF_B8.pkl','rb'))
print("\nstable 6->8?", counts)
# organize into anharmonic orbits
def orbit(z):
    oz = psub(ONE,z)
    iz = pdiv(ONE,z)
    return frozenset([z, oz, iz, pdiv(ONE,oz), pdiv(z, psub(z,ONE)), pdiv(psub(z,ONE), z)])
seen = set(); orbits = []
for z in sols:
    if z in seen: continue
    o = orbit(z) & set(sols.keys())
    # full orbit should be inside sols
    ofull = orbit(z)
    orbits.append(ofull)
    seen |= ofull
print("number of anharmonic orbits:", len(orbits))
print("all orbits S3-closed:", all(len(o) in (3,6) for o in orbits))
sizes = sorted(len(o) for o in orbits)
print("orbit sizes:", sizes)
# show representatives with small float in (0,1) where possible
reps = []
for o in orbits:
    cand = sorted(o, key=lambda z: (pfloat(z)<=0 or pfloat(z)>=1, abs(pfloat(z)-0.5)))
    reps.append(cand[0])
def show(z):
    f = factor_sunit(z)
    return "(%+d) p2^%d p3^%d eps^%d"%f
print("\nrepresentatives:")
for r in sorted(reps, key=pfloat):
    print("   z = %-14s ~ %+.10f  fac %s"%(str(r), pfloat(r), show(r)))
