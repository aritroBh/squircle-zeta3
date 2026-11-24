
import mpmath as mp
import json
mp.mp.dps = 80

def load(mm, full=0):
    data = json.load(open(f"ckpt_arc_m{mm}_f{full}.json"))
    # keep the geometric series (10 values)
    pts = [(mp.mpf(d), mp.mpf(v)) for d, v in data]
    pts.sort(key=lambda x: -x[0])
    return pts[:10]

def neville(pts):
    # extrapolate f(d) to d=0; pts list of (d, f); returns (value, errorest)
    n = len(pts)
    Q = [[mp.mpf(0)]*n for _ in range(n)]
    for i in range(n): Q[i][0] = pts[i][1]
    for j in range(1, n):
        for i in range(n-j):
            di, dj = pts[i][0], pts[i+j][0]
            Q[i][j] = (dj*Q[i][j-1] - di*Q[i+1][j-1])/(dj - di)
    return Q[0][n-1], abs(Q[0][n-1]-Q[0][n-2]), Q

def extract(pts, order=4):
    # returns A1..A_order with error estimates
    out = []
    rem = [(d, v) for d, v in pts]
    An = []
    for k in range(1, order+1):
        fk = [(d, v/d**k) for d, v in rem]
        val, err, Q = neville(fk)
        An.append((val, err))
        rem = [(d, v - val*d**k) for d, v in rem]
    return An

z3 = mp.zeta(3); pi2 = mp.pi**2; ln2 = mp.log(2)
for mm in [2, 3]:
    pts = load(mm)
    An = extract(pts, 4)
    print(f"=== m={mm} local ===")
    for k,(v,e) in enumerate(An, start=1):
        print(f"  A{k} = {mp.nstr(v, 30)}  (err {mp.nstr(e,4)})")
    A3 = An[2][0]
    print("  6*A3 =", mp.nstr(6*A3, 30))
    # PSLQ for zeta(3) coefficient with algebraic unit
    if mm == 2:
        tgt = mp.sqrt(2)*(z3 - pi2/3)
        print("  target sqrt2(z3-pi2/3):", mp.nstr(tgt, 30), " resid:", mp.nstr(A3-tgt, 6))
    if mm == 3:
        unit = mp.sqrt(3)
        basis = [z3, pi2*ln2, ln2**3, pi2, ln2**2, ln2, mp.mpf(1)]
        M = [A3/unit] + basis
        rel = mp.pslq(M, tol=mp.mpf('1e-22'), maxcoeff=10**12, maxsteps=10000)
        print("  PSLQ A3/sqrt3:", rel)
        if rel: print("   resid:", mp.nstr(mp.fsum(mp.mpf(r)*v for r,v in zip(rel,M)),6))
