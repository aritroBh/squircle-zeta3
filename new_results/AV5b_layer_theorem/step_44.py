import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
step = r'''
import mpmath as mp
import json
mp.mp.dps = 80

def load(path):
    data = json.load(open(path))
    pts = [(mp.mpf(d), mp.mpf(v)) for d, v in data]
    pts.sort(key=lambda x: -x[0])
    return pts

def neville(pts):
    n = len(pts)
    Q = [[mp.mpf(0)]*n for _ in range(n)]
    for i in range(n): Q[i][0] = pts[i][1]
    for j in range(1, n):
        for i in range(n-j):
            di, dj = pts[i][0], pts[i+j][0]
            Q[i][j] = (dj*Q[i][j-1] - di*Q[i+1][j-1])/(dj - di)
    errs = [abs(Q[0][k]-Q[0][k-1]) for k in range(1, n)]
    best = min(range(1,n), key=lambda k: errs[k-1])
    return Q[0][best], errs[best-1]

def extract(pts, order=4):
    An = []
    rem = [(d, v) for d, v in pts]
    for k in range(1, order+1):
        fk = [(d, v/d**k) for d, v in rem]
        val, err = neville(fk)
        An.append((val, err))
        rem = [(d, v - val*d**k) for d, v in rem]
    return An

z3 = mp.zeta(3); pi2 = mp.pi**2; ln2 = mp.log(2)
pts = load("ckpt_true_m3.json")
An = extract(pts, 4)
print("m=3 true:")
for k,(v,e) in enumerate(An,1): print("  A%d = %s (err %s)" % (k, mp.nstr(v,30), mp.nstr(e,3)))
print("  A1 target:", mp.nstr(2*mp.sqrt(3)*ln2, 30))
A3 = An[2][0]
basis = [z3, pi2*ln2, ln2**3, pi2, ln2**2, ln2, mp.mpf(1)]
M = [A3/mp.sqrt(3)] + basis
rel = mp.pslq(M, tol=mp.mpf('1e-25'), maxcoeff=10**10, maxsteps=10000)
print("  PSLQ A3/sqrt3:", rel)
if rel: print("   resid:", mp.nstr(mp.fsum(mp.mpf(r)*v for r,v in zip(rel,M)),6))
tgtP3 = -9*mp.sqrt(3)/8*z3 + mp.sqrt(3)*(3*pi2*ln2/8 + 2*ln2**3 - 3*pi2/4 - 12*ln2**2 + 12*ln2)
print("   A3 target:", mp.nstr(tgtP3/6, 30), " resid:", mp.nstr(A3-tgtP3/6, 6))
'''
with open("av5b_fit_true.py","w") as f: f.write(step)
print(open("av5b_fit_true.py").read()[:80])
exec(step)
