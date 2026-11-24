code = r'''
from squircle import setdps
from mpmath import mp, nstr, fabs, log, power, mpf, zeta, gamma, exp, fsum, pi, matrix, lu_solve
setdps(80)
s = 4
def J(u):
    tot = pi**2*u/12 - log(2)/2 + 1/(24*u)
    n = 0
    while True:
        n += 1
        t = log(1+exp(-2*pi**2*u*n))
        tot -= t
        if t < mpf(10)**(-75): break
    return tot
def inner(k, p):
    kf = mpf(k)
    return fsum([log(1+power(1-mpf(d)/kf, p)) for d in range(1, k)])
# for each u in {0.5, 1, 1.5, 2, 3, 4}: measure at p = 256, 512, 768 (k = u*p)
def extractA(u):
    ys = []
    for p in [256, 512, 768]:
        k = int(u*p)
        d0 = inner(k, p) - J(u)
        ys.append((mpf(1)/p, d0))
    # fit d0 = A1 x + A2 x^2 + A3 x^3
    A = matrix(3,3); b = matrix(3,1)
    for i,(x,y) in enumerate(ys):
        A[i,0] = x; A[i,1] = x**2; A[i,2] = x**3
        b[i,0] = y
    c = lu_solve(A, b)
    return c[0,0], c[1,0], c[2,0]
def M2(u): return fsum([d*d/(exp(d/u)+1) for d in range(1, 80)])
print("u     A1 meas        A1=-M2/(2u^2)    A2 meas")
for u in [0.5, 1.0, 1.5, 2.0, 3.0, 4.0]:
    a1, a2, a3 = extractA(u)
    a1m = -M2(u)/(2*u*u)
    print(f"u={u}: {nstr(a1, 12)}  {nstr(a1m, 12)}  A2={nstr(a2, 12)}")
'''
exec(code)
