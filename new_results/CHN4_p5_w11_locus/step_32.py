import random, pickle
def deg_monomials(d, n=5):
    out=[]
    def rec(i, rem, cur):
        if i==n-1: out.append(tuple(cur+[rem])); return
        for k in range(rem+1): rec(i+1, rem-k, cur+[k])
    rec(0,d,[]); return out
mons11 = deg_monomials(11)
print("deg-11 monomials:", len(mons11), "(expect 1365)")
seen=set(); orbits11=[]
for e in mons11:
    if e in seen: continue
    r=e[::-1]; orbits11.append((e,) if r==e else (e,r)); seen.add(e); seen.add(r)
print("reversal orbits deg-11:", len(orbits11), "(expect 693)")
orbit_of11={}
for i,o in enumerate(orbits11):
    for e in o: orbit_of11[e]=i
N11=len(orbits11)
def orbit_vals11(a):
    vals=[0]*N11
    for e in mons11:
        pw=1
        for x,ei in zip(a,e):
            if ei: pw*=x**ei
        vals[orbit_of11[e]]+=pw
    return vals
def rank_mod(rows,p):
    M=[[x%p for x in r] for r in rows]; nr=len(M); nc=len(M[0]); r=0
    for c in range(nc):
        piv=next((i for i in range(r,nr) if M[i][c]%p),None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; inv=pow(M[r][c]%p,p-2,p); M[r]=[(v*inv)%p for v in M[r]]
        for i in range(nr):
            if i!=r and M[i][c]%p:
                f=M[i][c]%p; M[i]=[(M[i][j]-f*M[r][j])%p for j in range(nc)]
        r+=1
        if r==nr: break
    return r
# select 693 points, no reversal-duplicates, rank 693
random.seed(31337)
p=(1<<61)-1
got=None
for attempt in range(15):
    seenpts=set(); pts=[]
    while len(pts)<N11:
        a=tuple(random.randint(1,10) for _ in range(5))
        if a in seenpts or a[::-1] in seenpts: continue
        seenpts.add(a); pts.append(a)
    rows=[orbit_vals11(a) for a in pts]
    if rank_mod(rows,p)==N11:
        got=pts; print("rank-693 set on attempt",attempt); break
if got:
    pickle.dump({'sel':got,'orbits11':orbits11,'mons11':mons11}, open('det11_points.pkl','wb'))
    print("saved det11_points.pkl with", len(got), "points")
else:
    print("no full-rank set found")
