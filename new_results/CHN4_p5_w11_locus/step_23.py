import sys, random, pickle
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
    r=e[::-1]; orbits.append((e,) if r==e else (e,r)); seen.add(e); seen.add(r)
orbit_of={}
for i,o in enumerate(orbits):
    for e in o: orbit_of[e]=i
NORB=len(orbits)
def orbit_vals_int(a):
    vals=[0]*NORB
    for e in mons7:
        pw=1
        for x,ei in zip(a,e):
            if ei: pw*= x**ei
        vals[orbit_of[e]]+=pw
    return vals

def rank_mod(rows, p):
    M=[ [x%p for x in r] for r in rows ]
    nr=len(M); nc=len(M[0]); r=0; 
    for c in range(nc):
        piv=None
        for i in range(r,nr):
            if M[i][c]%p: piv=i; break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]
        inv=pow(M[r][c]%p,p-2,p)
        M[r]=[(v*inv)%p for v in M[r]]
        for i in range(nr):
            if i!=r and M[i][c]%p:
                f=M[i][c]%p
                M[i]=[(M[i][j]-f*M[r][j])%p for j in range(nc)]
        r+=1
        if r==nr: break
    return r

p=(1<<61)-1
random.seed(123)
got=None
for attempt in range(30):
    pts=[tuple(random.randint(1,6) for _ in range(5)) for _ in range(NORB)]
    rows=[orbit_vals_int(a) for a in pts]
    if rank_mod(rows,p)==NORB:
        got=pts; print("invertible set on attempt",attempt); break
assert got
pickle.dump({'sel':got,'orbits':orbits,'mons7':mons7}, open('interp_points.pkl','wb'))
print("saved", len(got),"points; sample:",got[:2])
