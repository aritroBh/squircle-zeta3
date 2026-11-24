import pickle, random
import sympy as sp
d=pickle.load(open('interp_points.pkl','rb')); orbits=d['orbits']; mons7=d['mons7']
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
random.seed(999)
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
# find a failing set and get exact left-nullspace
fails=[]
for t in range(3):
    pts=[tuple(random.randint(1,9) for _ in range(5)) for _ in range(NORB)]
    rows=[orbit_vals_int(a) for a in pts]
    if rank_mod(rows,(1<<61)-1)<NORB:
        fails.append((pts,rows))
pts,rows=fails[0]
Msp=sp.Matrix(rows)
ns=Msp.T.nullspace()  # left nullspace
print("left nullity:", len(ns))
for vec in ns:
    # which point-indices are involved
    nz=[i for i,x in enumerate(vec) if x!=0]
    print(" kernel vector nonzero at point-indices:", nz[:10], "count", len(nz))
# Check: are there duplicate or reversal-duplicate points?
for i in range(len(pts)):
    for j in range(i+1,len(pts)):
        if pts[i]==pts[j] or pts[i]==pts[j][::-1]:
            print("dup/revdup at",i,j,pts[i])
