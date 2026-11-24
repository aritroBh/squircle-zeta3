import sys, itertools
sys.path.insert(0,'.')
from fractions import Fraction as F

def deg_monomials(d, n=5):
    out=[]
    def rec(i, rem, cur):
        if i==n-1:
            out.append(tuple(cur+[rem])); return
        for k in range(rem+1):
            rec(i+1, rem-k, cur+[k])
    rec(0,d,[])
    return out

def revex(e): return e[::-1]

mons7 = deg_monomials(7)
print("deg-7 monomials in 5 vars:", len(mons7), "(expect 330)")
# group into reversal orbits
seen=set(); orbits=[]
for e in mons7:
    if e in seen: continue
    r=revex(e)
    if r==e:
        orbits.append((e,)); seen.add(e)
    else:
        orbits.append((e,r)); seen.add(e); seen.add(r)
print("reversal orbits:", len(orbits), "(expect 170)")
print("singleton orbits:", sum(1 for o in orbits if len(o)==1))

# representative exponent per orbit (use min for determinism)
orbit_reps = [min(o) for o in orbits]
# map each monomial to its orbit index
orbit_of = {}
for i,o in enumerate(orbits):
    for e in o:
        orbit_of[e]=i

def eval_orbit_monomials(a):
    """value of each orbit-basis polynomial (sum of monomials in orbit) at point a (length5)."""
    # precompute powers
    vals=[F(0)]*len(orbits)
    for e in mons7:
        i=orbit_of[e]
        p=F(1)
        for x,ei in zip(a,e):
            if ei: p*= x**ei
        vals[i]+=p
    return vals

# quick sanity: at a=(1,1,1,1,1), each orbit value = size of orbit
v=eval_orbit_monomials((1,1,1,1,1))
print("check at all-ones: orbit values equal orbit sizes:", all(v[i]==len(orbits[i]) for i in range(len(orbits))))
