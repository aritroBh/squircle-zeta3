
from mpmath import mp, mpf, mpc, pi, sqrt, exp, log, nstr, sin, cos, gamma
mp.dps = 100

# ---------- exact series coefficients of the model G(eps) = (1-eps)^b e^{A/(1-eps)} ----------
def model_coeffs(A, b, K):
    # series arithmetic: G = e^A (1-eps)^b exp(A v), v = eps/(1-eps)
    # work with mpf/mpc polynomials of degree K
    one = mpc(1)
    # v = eps + eps^2 + ... 
    v = [mpc(0)]*(K+2)
    for j in range(1,K+2): v[j] = one
    # exp(A v) by series: e = 1 + Av + (Av)^2/2! + ...
    Av = [A*c for c in v]
    e = [mpc(0)]*(K+2); e[0] = 1
    term = [mpc(0)]*(K+2); term[0] = 1
    for m in range(1, K+2):
        # term *= Av/m  (Cauchy product)
        new = [mpc(0)]*(K+2)
        for i in range(K+1):
            if term[i] == 0: continue
            for j in range(K+2-i):
                if Av[j] == 0: continue
                new[i+j] += term[i]*Av[j]
        term = [c/m for c in new]
        for i in range(K+2): e[i] += term[i]
        if all(term[i] == 0 for i in range(K+2)): break
    # multiply by e^A
    e = [exp(A)*c for c in e]
    # multiply by (1-eps)^b = sum_j (-1)^j C(b,j) eps^j
    cb = [mpc(1)]
    for j in range(1,K+1):
        cb.append(cb[-1]*(-(b-(j-1))/j))
    g = [mpc(0)]*(K+1)
    for i in range(K+1):
        for j in range(K+1-i):
            g[i+j] += e[i]*cb[j]
    return g  # g[k] = [eps^k] G

A = mpc(0, pi/2)
K = 400
for b in [mpf('0.5'), mpf('1.5')]:
    g = model_coeffs(A, b, K)
    print(f"--- b = {b} ---")
    for k in [50, 100, 200, 400]:
        N = mpf(k+1)
        # saddle formula: (1/(2 i sqrt pi)) A^{(2b+1)/4} e^{A/2} N^{-(2b+3)/4} e^{2 sqrt(A N)} (decaying root)
        sqrtA = -sqrt(A)  # decaying branch: 2 sqrt(A N) = -sqrt(2 pi N) e^{i pi/4}
        two_sqrt_AN = 2*sqrtA*sqrt(N)
        pred = 1/(2*mpc(0,1)*sqrt(pi)) * A**((2*b+1)/4) * exp(A/2) * N**(-(2*b+3)/4) * exp(two_sqrt_AN)
        ratio = g[k]/pred
        print(f"  k={k:3d}: |g_k|={nstr(abs(g[k]),4)}  |pred|={nstr(abs(pred),4)}  ratio={nstr(ratio,8)}")
