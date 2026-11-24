import os, sys, time, random, pickle
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
from flint import nmod_mat
def deg_monomials(d, n=5):
    out=[]
    def rec(i, rem, cur):
        if i==n-1: out.append(tuple(cur+[rem])); return
        for k in range(rem+1): rec(i+1, rem-k, cur+[k])
    rec(0,d,[]); return out
mons12=deg_monomials(12); mons8=deg_monomials(8)
def make_orbits(mons):
    seen=set(); orbits=[]
    for e in mons:
        if e in seen: continue
        r=e[::-1]; orbits.append((e,) if r==e else (e,r)); seen.add(e); seen.add(r)
    return orbits
orbits12=make_orbits(mons12); orbits8=make_orbits(mons8)
orbit_of12={}; 
for i,o in enumerate(orbits12):
    for e in o: orbit_of12[e]=i
orbit_of8={}
for i,o in enumerate(orbits8):
    for e in o: orbit_of8[e]=i
N12=len(orbits12); N8=len(orbits8)
def orbit_vals(a, mons, orbit_of, N, p=None):
    vals=[0]*N
    for e in mons:
        pw=1
        for x,ei in zip(a,e):
            if ei: pw*=x**ei
        vals[orbit_of[e]]+=pw
    if p: vals=[v%p for v in vals]
    return vals
def rank_nmod(rows, p, nc):
    nr=len(rows)
    flat=[x for r in rows for x in r]
    M=nmod_mat(nr,nc,flat,p)
    return M.rank()
random.seed(20251029)
p1=(1<<61)-1; p2=(1<<63)-25
t0=time.time()
got=None
for attempt in range(10):
    seenpts=set(); pts=[]
    while len(pts)<N12:
        a=tuple(random.randint(1,10) for _ in range(5))
        if a in seenpts or a[::-1] in seenpts: continue
        seenpts.add(a); pts.append(a)
    rows=[orbit_vals(a,mons12,orbit_of12,N12,p1) for a in pts]
    if rank_nmod(rows,p1,N12)==N12:
        got=pts; print(f"rank-924 set on attempt {attempt} [{time.time()-t0:.0f}s]"); break
if got:
    # verify with second prime
    rows2=[orbit_vals(a,mons12,orbit_of12,N12,p2) for a in got]
    r2=rank_nmod(rows2,p2,N12)
    print(f"rank mod p2: {r2}")
    # choose interpolation subset for degree-8: greedy independent rows mod p1
    rows8=[orbit_vals(a,mons8,orbit_of8,N8,p1) for a in got]
    # greedy: add point if it increases rank
    sel_idx=[]; cur=[]
    r=0
    for i,row in enumerate(rows8):
        cur.append(row)
        rr=rank_nmod(cur,p1,N8)
        if rr>r:
            sel_idx.append(i); r=rr
            if r==N8: break
    print(f"deg-8 interpolation subset: {len(sel_idx)} points, rank {r}")
    pickle.dump({'sel':got,'orbits12':orbits12,'mons12':mons12,'orbits8':orbits8,'mons8':mons8,
                 'interp_idx':sel_idx}, open('det12_points.pkl','wb'))
    print("saved det12_points.pkl", f"[{time.time()-t0:.0f}s]")
