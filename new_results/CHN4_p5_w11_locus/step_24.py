import pickle
d=pickle.load(open('interp_points.pkl','rb'))
sel=d['sel']; orbits=d['orbits']; mons7=d['mons7']
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
    M=[[x%p for x in r] for r in rows]; nr=len(M); nc=len(M[0]); r=0
    for c in range(nc):
        piv=next((i for i in range(r,nr) if M[i][c]%p), None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]
        inv=pow(M[r][c]%p,p-2,p); M[r]=[(v*inv)%p for v in M[r]]
        for i in range(nr):
            if i!=r and M[i][c]%p:
                f=M[i][c]%p; M[i]=[(M[i][j]-f*M[r][j])%p for j in range(nc)]
        r+=1
        if r==nr: break
    return r
rows=[orbit_vals_int(a) for a in sel]
for p in [(1<<61)-1, (1<<63)-25, 10**18+9]:
    print(f"rank mod {p}:", rank_mod(rows,p))
# Investigate the 15-attempt issue: check rank of a few fresh random sets
import random
random.seed(999)
fails=0
for t in range(10):
    pts=[tuple(random.randint(1,9) for _ in range(5)) for _ in range(NORB)]
    rk=rank_mod([orbit_vals_int(a) for a in pts], (1<<61)-1)
    if rk<NORB: fails+=1
print("fresh sets failing full rank (of 10, larger range):", fails)
