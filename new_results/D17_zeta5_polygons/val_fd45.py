
import mpmath as mp, sys, json, time
mp.mp.dps = 120
sys.path.insert(0, "research/08_theory/D17_zeta5_polygons"); sys.path.insert(0, "research/08_theory/D15b_polygon_law")
import poly_coeffs5 as p5
from poly_coeffs import P_direct, regular_zeros

def fd_derivs(m, h, K):
    zeros = regular_zeros(m); w=[1]*m; lams=[0]*m
    xs = [j*h for j in range(-K,K+1)]
    ys = [P_direct(zeros, w, lams, 1+x) for x in xs]
    M = mp.matrix([[x**j for j in range(2*K+1)] for x in xs])
    cf = mp.lu_solve(M, mp.matrix(ys))
    from math import factorial
    return [cf[j]*factorial(j) for j in range(6)]

print("dps:", mp.mp.dps)
t0=time.time()
for m in [2,3,4]:
    A=mp.pi/m
    # integral values P^{(k)}(1) = k! * 4m int_0^{A/2} D_k
    pc = p5.regular_P_coeffs5(m, 5)
    from math import factorial
    integ = [pc[k]*factorial(k) for k in range(6)]
    # FD
    fdv = fd_derivs(m, mp.mpf('0.006'), 8)
    print("m=%d (%.0fs):"%(m, time.time()-t0))
    for k in range(1,6):
        rel = abs(fdv[k]-integ[k])/max(mp.mpf(1),abs(integ[k]))
        print("  P^(%d)(1): integral=%s  FD=%s  relres=%s"%(k, mp.nstr(integ[k],12), mp.nstr(fdv[k],12), mp.nstr(rel,2)))
