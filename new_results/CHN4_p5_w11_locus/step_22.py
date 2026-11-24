import sys, random, pickle
sys.path.insert(0,'.')
from fractions import Fraction as F
from flint import fmpz_mod_mat, fmpz_mod_ctx

def deg_monomials(d, n=5):
    out=[]
    def rec(i, rem, cur):
        if i==n-1: out.append(tuple(cur+[rem])); return
        for k in range(rem+1): rec(i+1, rem-k, cur+[k])
    rec(0,d,[]); return out
mons7 = deg_monomials(7)
seen=set(); orbits=[]
for e in mons7:
    if e in seen: continue
    r=e[::-1]; orbits.append((e,) if r==e else (e,r)); seen.add(e); seen.add(r)
orbit_of={}
for i,o in enumerate(orbits):
    for e in o: orbit_of[e]=i
NORB=len(orbits)
def orbit_vals_int(a):
    vals=[0]*NORB
    for e in mons7:
        p=1
        for x,ei in zip(a,e):
            if ei: p*= x**ei
        vals[orbit_of[e]]+=p
    return vals

p=2**62-57
ctx=fmpz_mod_ctx(p)
random.seed(123)
got=None
for attempt in range(20):
    pts=[tuple(random.randint(1,6) for _ in range(5)) for _ in range(NORB)]
    M=fmpz_mod_mat(ctx, NORB, NORB, [x%p for a in pts for x in orbit_vals_int(a)])
    if M.rank()==NORB:
        got=pts; print("invertible set found on attempt", attempt); break
assert got is not None
pickle.dump({'sel':got,'orbits':orbits,'mons7':mons7}, open('interp_points.pkl','wb'))
print("saved", len(got), "points")
print("sample:", got[:3])
