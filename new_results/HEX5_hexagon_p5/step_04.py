import sys, time
from mpmath import mp, mpf, cos, sin, atan2, pi, log, sqrt, exp, factorial, tanh, sinh, cosh, nstr

# ---- series machinery (copied from ../CX3k_review_draft8/polygon_engine.py) ----
def s_mul(a, b, N):
    out = [mpf(0)]*(N+1)
    for i in range(min(len(a), N+1)):
        ai = a[i]
        if ai == 0: continue
        for j in range(min(len(b), N+1-i)):
            out[i+j] += ai*b[j]
    return out

def s_log(a, N):
    out = [mpf(0)]*(N+1)
    out[0] = log(a[0])
    x = [mpf(0)] + [a[i]/a[0] for i in range(1, N+1)]
    xk = [mpf(0)]*(N+1); xk[0] = mpf(1)
    for k in range(1, N+1):
        xk = s_mul(xk, x, N)
        c = mpf(1)/k if k % 2 == 1 else -mpf(1)/k
        for i in range(1, N+1):
            out[i] += c*xk[i]
    return out

def s_exp(a, N):
    out = [mpf(0)]*(N+1); out[0] = mpf(1)
    for k in range(1, N+1):
        acc = mpf(0)
        for j in range(1, k+1):
            acc += j*a[j]*out[k-j]
        out[k] = acc/k
    return out

def s_pow_series_exp(a, w, N):
    lg = s_log(a, N)
    t = s_mul(w, lg, N)
    t0 = t[0]
    t[0] = mpf(0)
    e = s_exp(t, N)
    return [exp(t0)*c for c in e]

VS = None
N = 7
def integrand_series(th):
    A = [mpf(0)]*(N+1); Ap = [mpf(0)]*(N+1)
    ct, st = cos(th), sin(th)
    tiny = mpf(2)**(-mp.prec + 10)
    for (vx, vy) in VS:
        c = vx*ct + vy*st
        cp = -vx*st + vy*ct
        sgn = mpf(1) if c >= 0 else mpf(-1)
        g = abs(c)
        if g < tiny: g = tiny
        gp = cp*sgn
        el = log(g)
        elk = mpf(1); elkm1 = mpf(1); fk = mpf(1)
        for k in range(N+1):
            if k > 0:
                elkm1 = elk; elk *= el; fk *= k
            A[k] += g*elk/fk
            Ap[k] += (gp*elk + (k*elkm1*gp if k > 0 else mpf(0)))/fk
    w = [mpf(-1)] + [mpf((-1)**(k+1)) for k in range(1, N+1)]
    r = s_pow_series_exp(A, w, N)
    wm1 = list(w); wm1[0] = w[0] - 1
    Swm1 = s_pow_series_exp(A, wm1, N)
    rp = s_mul(s_mul(w, Swm1, N), Ap, N)
    tot = [r[i]*0 for i in range(N+1)]
    tot = s_mul(r, r, N)
    tmp = s_mul(rp, rp, N)
    tot = [tot[i]+tmp[i] for i in range(N+1)]
    lg = s_log(tot, N)
    half = [c/2 for c in lg]
    h0 = half[0]; half[0] = mpf(0)
    e = s_exp(half, N)
    return [exp(h0)*c for c in e]

def one_wedge(kmax, eps, logf):
    a = pi/6; b = pi/2
    mid = (a+b)/2; half = (b-a)/2
    h = mpf(2); prev = None; diff = None
    for level in range(20):
        t0 = time.time()
        h = h/2
        acc = [mpf(0)]*(kmax+1)
        k = 0; npts = 0
        while True:
            t = k*h
            arg = (pi/2)*sinh(t)
            u = tanh(arg)
            dx = half*(pi/2)*cosh(t)/(cosh(arg)**2)*h
            if k == 0:
                F = integrand_series(mid); mult = 1; npts += 1
            else:
                F = integrand_series(mid + half*u); mult = 2; npts += 2
            for j in range(kmax+1):
                acc[j] += mult*dx*F[j]
            contrib = max(abs(dx*f) for f in F)
            k += 1
            if contrib < eps and k > 8: break
            if k > 5000000: raise RuntimeError("no convergence")
        if prev is not None:
            diff = max(abs(acc[j]-prev[j]) for j in range(kmax+1))
        else:
            diff = None
        logf("  level %2d: pts=%6d  acc[5]=%s  diff=%s  (%.1fs)" % (
            level, npts, nstr(acc[5], 12), (nstr(diff,4) if diff else '-'), time.time()-t0))
        if diff is not None and diff < eps:
            return acc, diff
        prev = acc
    return prev, diff

def main():
    global VS
    DPS = int(sys.argv[1]); tag = sys.argv[2]
    mp.dps = DPS + 15
    VS = [(mpf(1)/2, sqrt(3)/2), (-mpf(1)/2, sqrt(3)/2), (mpf(-1), mpf(0))]
    eps = mpf(10)**(-(DPS - 5))
    t0 = time.time()
    logf = lambda s: (print(s, flush=True))
    logf("DPS=%d N=%d  computing one wedge (pi/6, pi/2), x6 symmetry" % (DPS, N))
    acc, diff = one_wedge(5, eps, logf)
    res = [6*factorial(k)*acc[k] for k in range(6)]
    out = ["DPS=%d elapsed=%.1fs final level diff(k<=5)=%s" % (DPS, time.time()-t0, nstr(diff,6))]
    for k in range(6):
        out.append("P^(%d)(1) = %s" % (k, nstr(res[k], DPS)))
    out.append("P5/sqrt3 = " + nstr(res[5]/sqrt(3), DPS))
    txt = "\n".join(out)
    logf(txt)
    open(tag, 'w').write(txt + "\n")

if __name__ == '__main__':
    main()
