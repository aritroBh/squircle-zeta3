
"""Direct U_q(p,s) = sum_{k>j} k^{-s} ln^q(1+(j/k)^p), boundary-layer direct sum + EM tail."""
from mpmath import mp, mpf, zeta, bernoulli, fac, binomial, fsum, log, power, gamma, loggamma, exp

def _lnq_derivs_1(p, q, R):
    """(d/dx)^m ln^q(1+x^p) at x=1, m=1..R. Via u = ((1+t)^p-1)/2, L = ln(1+u), psi = L^q."""
    bk = [mpf(0)]*(R+1); b = mpf(1)
    for k in range(R+1):
        if k > 0: b *= (p-(k-1))/k
        bk[k] = b
    u = [bki/2 for bki in bk]; u[0] = mpf(0)
    # L = ln(1+u): L' (1+u) = u': coefficients: (n+1)L_{n+1} + sum_{j=1}^n (n+1-j) u_j L_{n+1-j} = (n+1) u_{n+1}
    L = [mpf(0)]*(R+1)
    for n in range(R):
        acc = (n+1)*u[n+1]
        for j in range(1, n+1):
            acc -= (n+1-j)*u[j]*L[n+1-j]
        L[n+1] = acc/(n+1)
    # psi = L^q: series power: psi' = q L' psi / L: L psi' = q L' psi
    P = [mpf(0)]*(R+1); P[0] = L[0]**q if L[0] != 0 else mpf(0)
    # L[0] = ln((1+1)/2)=0? u[0]=0 -> L[0]=ln1=0. Handle via L = ln2 + L~ hmm
    # Actually psi(x)=ln^q(1+x^p): at t=0: ln^q(2) != 0. Redo: psi = (ln2 + L)^q with L[0]=0:
    # psi = sum_{i=0}^q binom(q,i) (ln2)^{q-i} L^i : compute powers of L series
    ln2 = log(2)
    def series_mul(A, B, R):
        C = [mpf(0)]*(R+1)
        for i in range(R+1):
            if A[i] == 0: continue
            for j in range(R+1-i):
                if B[j] != 0: C[i+j] += A[i]*B[j]
        return C
    psi = [mpf(0)]*(R+1); psi[0] = ln2**q
    Lpow = [mpf(0)]*(R+1); Lpow[0] = 1
    for i in range(1, q+1):
        Lpow = series_mul(Lpow, L, R)
        for k in range(R+1):
            psi[k] += binomial(q, i)*ln2**(q-i)*Lpow[k]
    return [psi[m]*fac(m) for m in range(R+1)]

def Uq_direct(p, s, q, K=None, R=10, info=False):
    p = mpf(p); s = mpf(s); q = int(q)
    if K is None:
        K = max(400, int(6*p))
    I0q = fsum([_stirling_ratio(n, q) / (n*p+1) for n in range(q, 300)])
    # wait: int_0^1 ln^q(1+x^p) dx = sum_{n>=q} (-1)^{n-q} |s(n,q)| q!/n! /(np+1)
    psi0 = mpf(0); psi1 = log(2)**q
    d1 = _lnq_derivs_1(p, q, 2*R+1)
    def HT(K, sig):
        return zeta(sig) - fsum([power(m, -sig) for m in range(1, K+1)])
    P = mpf(0)
    for k in range(2, K+1):
        inner = mpf(0)
        kf = mpf(k)
        # sum over d = k-j (boundary layer): ln^q(1+(1-d/k)^p)
        for d in range(1, k):
            inner += log(1+power(1-mpf(d)/kf, p))**q
        P += power(kf, -s)*inner
    tail = I0q*HT(K, s-1) - (psi0+psi1)/2*HT(K, s)
    for r in range(1, R+1):
        tail += bernoulli(2*r)/fac(2*r) * d1[2*r-1] * HT(K, s-1+2*r)
    # cusp terms: q! (-1)^{n-q} |s(n,q)|/n! zeta(-np) HT(K, s+np), n=1..4
    for n in range(q, q+6):
        if n < q: continue
        C = ((-1)**(n-q)) * stirling1u(n, q)*fac(q)/fac(n) * zeta(-n*p)
        if C != 0:
            tail += C*HT(K, s+n*p)
    T = P + tail
    if info: return T, dict(P=P, tail=tail, I0q=I0q)
    return T

def stirling1u(n, q):
    if n == q: return 1
    S = [[0]*(q+1) for _ in range(n+1)]
    S[0][0] = 1
    for a in range(1, n+1):
        for b in range(1, min(a,q)+1):
            S[a][b] = S[a-1][b-1] + (a-1)*S[a-1][b]
    return S[n][q]

def _stirling_ratio(n, q):
    # |s(n,q)| * q!/n! with sign (-1)^{n-q}
    return ((-1)**(n-q)) * stirling1u(n, q) * fac(q)/fac(n)
