import json, time
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, log, nstr, polyroots
mp.dps = 120
f450 = [mpf(v) for v in json.load(open('step2_num450_N180.json'))['f450']]
sq2 = sqrt(2)
f1 = [f450[n] - sq2*((-1)**n) for n in range(len(f450))]
b1 = [f1[n]/factorial(n) for n in range(len(f1))]
b  = [f450[n]/factorial(n) for n in range(len(f450))]

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

t0=time.time()
mp.dps = 200
print("building PB1[70/70] and [90/90] from 450-digit data...")
num70, den70 = pade_coeffs(b1[:141], 70, 70)
num90, den90 = pade_coeffs(b1[:181], 90, 90)
print("built, %.1fs"%(time.time()-t0))

def poles_zeros(num, den):
    pr = polyroots([den[k] for k in range(len(den)-1,-1,-1)], maxsteps=400)
    zr = polyroots([num[k] for k in range(len(num)-1,-1,-1)], maxsteps=400)
    return pr, zr

for tag,(num,den) in [('70/70',(num70,den70)),('90/90',(num90,den90))]:
    pr, zr = poles_zeros(num, den)
    zr = list(zr)
    prs = sorted(pr, key=lambda r: abs(r))
    print(f"\n==== PB1[{tag}] poles nearest origin (with dist to nearest zero) ====")
    for r in prs[:8]:
        dz = min(abs(r-z) for z in zr)
        print("   pole %s   |t|=%s   d_zero=%s" % (nstr(r,16), nstr(abs(r),12), nstr(dz,3)))
json.dump({'num':[nstr(v,180) for v in num90],'den':[nstr(v,180) for v in den90]}, open('step_09_pade90.json','w'))
print("saved step_09_pade90.json, %.1fs total"%(time.time()-t0))
