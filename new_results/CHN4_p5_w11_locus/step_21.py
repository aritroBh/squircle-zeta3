import sys, itertools, random, pickle
sys.path.insert(0,'.')
from fractions import Fraction as F

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
    r=e[::-1]
    orbits.append((e,) if r==e else (e,r)); seen.add(e); seen.add(r)
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

# pool of small positive-integer points
random.seed(7)
pool=set()
while len(pool)<4000:
    pool.add(tuple(random.randint(1,5) for _ in range(5)))
pool=list(pool)

# select 170 independent rows over GF(p)
from flint import fmpz_mod_mat, fmpz_mod_ctx
p=2**63-25
ctx=fmpz_mod_ctx(p)
sel=[]; rows_mod=[]
rank=0
for a in pool:
    v=orbit_vals_int(a)
    cand=rows_mod+[v]
    M=fmpz_mod_mat(ctx, len(cand), NORB, [x%p for r in cand for x in r])
    rk=M.rank()
    if rk>rank:
        sel.append(a); rows_mod.append(v); rank=rk
        if rank==NORB: break
print("selected points:", len(sel), "rank:", rank)
pickle.dump({'sel':sel,'orbits':orbits,'mons7':mons7}, open('interp_points.pkl','wb'))
print("saved interp_points.pkl")
