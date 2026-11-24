
# step_03.py : (a) refined check of corner-layer identity at p=100
# (b) TASK 2 core: extract deficit coefficients d1..d7 from A(p) series at 100 dps
#     and test  d2 = -zeta(2)/2,  d3 = +zeta(3)/6,  d4 = 0,  d5 = -(1/2880) int (Phi^4-u^4).
import mpmath as mp
mp.mp.dps = 105

def A_series(p, tol=mp.mpf('1e-101')):
    p = mp.mpf(p); s = mp.mpf(0); m = 0
    while True:
        m += 1
        t = mp.gamma(1 + m/p)**2 / (m**2 * mp.gamma(1 + 2*m/p))
        s += t
        if m > 20 and abs(t) < tol: break
    return s

def Phi(u): return -u + mp.log(1 - mp.exp(-u))

def Fser(eps, u, N=30):
    # F(eps,u) via log(-Phi/u) + L(eps*Phi) - L(-eps*u), L(z)=log((e^z-1)/z)
    # computed stably: F = log((1-exp(eps*Phi))/(1-exp(-eps*u)))
    a = eps*Phi(u); b = -eps*u
    # use sinh form: (e^a-1)/a = e^{a/2} sinh(a/2)/(a/2)
    def L(z):
        if abs(z) < mp.mpf('1e-8'):
            return z/2 + z**2/24 - z**4/2880 + z**6/181440
        return mp.log(mp.expm1(z)/z)
    return mp.log(-Phi(u)/u) + L(a) - L(b)

def deficit_int(p):
    eps = mp.mpf(1)/p
    f = lambda u: Fser(eps,u)
    pts = [mp.mpf('1e-10'), mp.mpf('1e-6'), mp.mpf('1e-3'), mp.mpf('0.1'), 1, 5, 20, mp.inf]
    return eps*mp.quad(f, pts)

p = 100
d1 = mp.zeta(2) - A_series(p); d2 = deficit_int(p)
print('refined check p=100: |series - cornerint| =', mp.nstr(abs(d1-d2),3))

# ---- extract coefficients. A(1/eps) = zeta2 - sum_{k>=1} d_k eps^k.
# Use eps_j = base/j and least-squares-free Richardson: evaluate deficit at many eps,
# fit d1..d8 (with d4 expected 0) by solving linear system on successive eps.
ps = [12,14,16,18,20,24,28,32,40,48,56,64,80,96,128]
Ds = [mp.zeta(2) - A_series(p) for p in ps]
es = [mp.mpf(1)/p for p in ps]
# solve for d1..d7 using consecutive eps values (Vandermonde, well-conditioned enough at 105 dps with small matrix)
import itertools
def extract(eps_list, def_list, k):
    M = mp.matrix(len(eps_list), k)
    for i,e in enumerate(eps_list):
        for j in range(k):
            M[i,j] = e**(j+1)
    b = mp.matrix(def_list)
    sol = mp.lu_solve(M, b)
    return [sol[j] for j in range(k)]

for k,sl in [(4, slice(0,4)), (6, slice(3,9)), (7, slice(8,15))]:
    sol = extract(es[sl], Ds[sl], k)
    print(f"k={k} fit:", [mp.nstr(s,18) for s in sol])

# candidate exact values
cand = {1: None, 2: -mp.zeta(2)/2, 3: mp.zeta(3)/6, 4: mp.mpf(0)}
print('zeta(2)/2 =', mp.nstr(mp.zeta(2)/2,20))
print('zeta(3)/6 =', mp.nstr(mp.zeta(3)/6,20))
# d5 candidate: -(1/2880) int_0^inf (Phi^4 - u^4) du  and MZV guess -zeta(5)/3 + zeta(2)zeta(3)/3
I = mp.quad(lambda u: Phi(u)**4 - u**4, [mp.mpf('1e-9'), mp.mpf('1e-4'), mp.mpf('1e-2'), 1, 5, 20, mp.inf])
print('-(1/2880) int(Phi^4-u^4) =', mp.nstr(-I/2880, 25))
print('(-zeta(5)+zeta(2)zeta(3))/3 =', mp.nstr((-mp.zeta(5)+mp.zeta(2)*mp.zeta(3))/3, 25))
