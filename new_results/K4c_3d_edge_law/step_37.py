
pslq_code = r'''"""pslq_lib.py -- PSLQ integer relation finder (Bailey-Broadhurst), mpmath."""
import mpmath as mp
from mpmath import mpf

def pslq(x, tol=None, gam=None, maxsteps=100000):
    n = len(x)
    if gam is None: gam = mp.sqrt(mpf(4)/3)
    if tol is None: tol = mpf(10)**(-mp.mp.dps*2/3)
    nx = mp.sqrt(sum(xi**2 for xi in x))
    y = [xi/nx for xi in x]
    A = [[mpf(1) if i==j else mpf(0) for j in range(n)] for i in range(n)]
    B = [[mpf(1) if i==j else mpf(0) for j in range(n)] for i in range(n)]
    s = [mpf(0)]*n
    for k in range(n-1, -1, -1):
        s[k] = mp.sqrt(sum(y[i]**2 for i in range(k, n)))
    H = [[mpf(0)]*(n-1) for _ in range(n)]
    for j in range(n-1):
        for i in range(n):
            if i < j: H[i][j] = mpf(0)
            elif i == j: H[i][j] = s[j+1]/s[j]
            else: H[i][j] = -y[i]*y[j]/(s[j]*s[j+1])
    def swap_rows(M, r):
        M[r], M[r+1] = M[r+1], M[r]
    for step in range(maxsteps):
        best = -1; r = 0
        for i in range(n-1):
            val = gam**i * abs(H[i][i])
            if val > best: best = val; r = i
        y[r], y[r+1] = y[r+1], y[r]
        swap_rows(H, r); swap_rows(A, r)
        for k in range(n):
            B[k][r], B[k][r+1] = B[k][r+1], B[k][r]
        if r < n-2:
            t0 = mp.sqrt(H[r][r]**2 + H[r][r+1]**2)
            if t0 == 0: return None, None
            t1 = H[r][r]/t0; t2 = H[r][r+1]/t0
            for i in range(r, n):
                t3, t4 = H[i][r], H[i][r+1]
                H[i][r] = t1*t3 + t2*t4
                H[i][r+1] = -t2*t3 + t1*t4
        for i in range(r+1, n):
            for j in range(min(i, n-2), -1, -1):
                if H[j][j] == 0: continue
                t = mp.nint(H[i][j]/H[j][j])
                if t != 0:
                    y[j] += t*y[i]
                    for k in range(j+1):
                        H[i][k] -= t*H[j][k]
                    for k in range(n):
                        A[i][k] -= t*A[j][k]
                        B[k][j] += t*B[k][i]
        # termination
        imin = min(range(n), key=lambda i: abs(y[i]))
        if abs(y[imin]) < tol:
            rel = [B[j][imin] for j in range(n)]
            g = rel[0]
            for t in rel[1:]:
                from math import gcd
                g = gcd(int(g), int(t))
            if g != 0:
                rel = [t//g for t in rel]
            return [int(t) for t in rel], abs(y[imin])
    return None, None
'''
with open('pslq_lib.py','w') as f: f.write(pslq_code)
import mpmath as mp
mp.mp.dps = 50
from pslq_lib import pslq
# test: ln4 - 2 ln2 = 0
x = [mp.log(4), mp.log(2), mpf(1)]
rel, res = pslq(x, tol=mpf('1e-40'))
print("test1 (ln4,ln2,1):", rel, "resid:", mp.nstr(res,3))
# test: zeta3 relation? x = [zeta(3), pi*... just check no false positive
rel, res = pslq([mp.zeta(3), mp.pi**3, mpf(1)], tol=mpf('1e-38'))
print("test2 (zeta3,pi^3,1):", rel, "resid:", mp.nstr(res,3) if res else None)
# test: q for cubo should give (2, -1, -4, 0): q = 1/2 + 2 ln2? q=1.8862944
import json
J = json.load(open("research/08_theory/K4b_cubo_finish/tetra_Sdk_final.json"))
S = [mpf(s_) for s_ in J["Sdk"]]
q = S[1]/S[0]
rel, res = pslq([q, mpf(1), mp.log(2), mp.log(3)], tol=mpf('1e-40'))
print("test3 cubo q on (1,ln2,ln3):", rel, "resid:", mp.nstr(res,3))
