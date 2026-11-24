from mpmath import mp, mpf, mpc, pi, sqrt, nstr, polyroots
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

# 1-cut model: Li_{3/2}(z), z = 2it/pi;  expect poles on z in [1,inf) i.e. t = -i*pi*z/2
tp = 2/pi
N = 25
a1 = [mpf(0)] + [mpc(0,1)**n * mpf(n)**mpf('-1.5')*tp**n for n in range(1,2*N+1)]
num, den = pade_coeffs(a1, N, N)
try:
    rts = polyroots([den[k] for k in range(len(den)-1,-1,-1)], maxsteps=200)
    prs = sorted(rts, key=lambda r: abs(r))
    print("1-cut model Li_{3/2}(2it/pi) [%d/%d] poles (mpmath polyroots, dps100):"%(N,N))
    for r in prs[:8]:
        z = 2j*r/pi
        print("   t = %s   |t|=%s   z=2it/pi=%s" % (nstr(r,14), nstr(abs(r),10), nstr(z,10)))
except Exception as e:
    print("polyroots failed:", e)
