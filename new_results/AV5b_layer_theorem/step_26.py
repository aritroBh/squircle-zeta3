
import mpmath as mp
import json
mp.mp.dps = 80

def load(mm, full):
    data = json.load(open(f"ckpt_arc_m{mm}_f{full}.json"))
    return [(mp.mpf(d), mp.mpf(v)) for d, v in data]

def fit_coeffs(pts, n=5):
    # solve for A1..An in DeltaP = sum Ak d^k
    M = mp.matrix(len(pts), n)
    b = mp.matrix(len(pts), 1)
    for i,(d,v) in enumerate(pts):
        for k in range(1, n+1):
            M[i, k-1] = d**k
        b[i] = v
    if len(pts) == n:
        return mp.lu_solve(M, b)
    else:
        return mp.lu_solve(M.T*M, M.T*b)

z3 = mp.zeta(3)
for mm, full in [(2,0),(3,0),(4,0)]:
    pts = load(mm, full)
    A = fit_coeffs(pts, 5)
    A1,A2,A3,A4,A5 = [A[i] for i in range(5)]
    print(f"=== m={mm} full={full} ===")
    print("  A1..A5:", [mp.nstr(x, 18) for x in [A1,A2,A3,A4,A5]])
    print("  6*A3 (model Ppp/1) =", mp.nstr(6*A3, 25))
    print("  6*A3/zeta(3) =", mp.nstr(6*A3/z3, 22))
    C_m = 48*mm*mp.sin(mp.pi/(2*mm))**4*mp.cos(mp.pi/(2*mm))*(3*mp.sin(mp.pi/(2*mm))**2-1)
    print("  C_m =", mp.nstr(C_m, 22), "  6A3 - C_m*z3 =", mp.nstr(6*A3 - C_m*z3, 8))
