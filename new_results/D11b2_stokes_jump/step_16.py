from mpmath import mp, mpf, mpc, pi, nstr, findroot
import numpy as np
mp.dps = 100
def pade_coeffs(a, L, M):
    A = mp.matrix(M, M); rhs = mp.matrix(M, 1)
    for i in range(M):
        for j in range(M):
            k = L + i - j
            A[i, j] = a[k] if k >= 0 else mpf(0)
        rhs[i] = -a[L + 1 + i]
    q = (A**-1)*rhs
    den = [mpf(1)] + [q[j] for j in range(M)]
    num = [sum(den[j]*a[n-j] for j in range(0, min(M, n)+1)) for n in range(L+1)]
    return num, den
def peval(c, t):
    v = c[0]
    for k in range(1,len(c)): v += c[k]*t**k
    return v
def refine_roots(den, tol=mpf('1e-40')):
    c = [complex(den[k]) for k in range(len(den)-1,-1,-1)]
    rts = np.roots(c)
    out = []
    f = lambda t: peval(den, t)
    for r in rts:
        try:
            rt = findroot(f, mpc(r.real, r.imag), solver='mnewton', tol=1e-30, maxsteps=50)
            out.append(rt)
        except Exception:
            out.append(('FAIL', r))
    return out
tp = 2/pi
N = 25
a1 = [mpf(0)] + [mpc(0,1)**n * mpf(n)**mpf('-1.5')*tp**n for n in range(1,2*N+1)]
num, den = pade_coeffs(a1, N, N)
rr = refine_roots(den)
ok = [r for r in rr if not isinstance(r, tuple)]
fail = [r for r in rr if isinstance(r, tuple)]
print("1-cut Li_{3/2} [%d/%d]: refined %d roots, %d failed" % (N,N,len(ok),len(fail)))
prs = sorted(ok, key=lambda r: abs(r))
for r in prs[:8]:
    z = mpc(0,2)*r/pi
    resid = abs(peval(den, r))
    print("   t = %s  |t|=%s  z=%s  |den|=%s" % (nstr(r,14), nstr(abs(r),8), nstr(z,8), nstr(resid,2)))
