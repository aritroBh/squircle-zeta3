code = r'''
from Uq import *
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, fsum, mpf, zeta, gamma, fac, matrix, lu_solve
setdps(150)
s = 4
def U1_formula(p, s, nmax=4):
    tot = mpf(0)
    for n in range(1, nmax+1):
        tot += ((-1)**(n+1))/n * zetaEZ(s, p*n)
    return tot
u1_pred = gamma(s-1)*zeta(s-1)*zeta(s)*(1-power(2,1-s))
print("u1 prediction:", nstr(u1_pred, 30))
vals = []
for p_ in [16, 24, 32, 48, 64, 96]:
    u = U1_formula(p_, s)
    vals.append((p_, u))
    print(f"p={p_}: U1={nstr(u, 35)}   p^3 U1={nstr(u*p_**(s-1), 25)}")
def lstsq(xs, ys, Q):
    n = len(xs)
    A = matrix(n, Q+1)
    for i in range(n):
        xp = mpf(1)
        for q in range(Q+1):
            A[i,q] = xp; xp *= xs[i]
    AtA = matrix(Q+1, Q+1); Atb = matrix(Q+1, 1)
    for a in range(Q+1):
        for b in range(Q+1):
            AtA[a,b] = sum(A[i,a]*A[i,b] for i in range(n))
        Atb[a,0] = sum(A[i,a]*ys[i] for i in range(n))
    c = lu_solve(AtA, Atb)
    return [c[q,0] for q in range(Q+1)]
# fit U1(p) = p^{-(s-1)} (u0 + u1/p + u2/p^2 + u3/p^3)
xs = [mpf(1)/p_ for p_,_ in vals]
ys = [u*power(p_, s-1) for p_, u in vals]
c = lstsq(xs, ys, 4)
print("fit u0:", nstr(c[0], 30), " diff:", nstr(c[0]-u1_pred, 6))
print("fit u1:", nstr(c[1], 25))
print("fit u2:", nstr(c[2], 22))
print("fit u3:", nstr(c[3], 18))
print("fit u4:", nstr(c[4], 14))
'''
exec(code)
with open("progress.log","a") as f: f.write("step_41: U_1(p) asymptotic fit at dps=150\n")
