
import mpmath as mp
mp.mp.dps = 120
from math import factorial

N = 90
Zv = {n: mp.zeta(n) for n in range(2, N+3)}

def series_mul(a, b, N):
    c = [mp.mpf(0)]*(N+1)
    for i in range(N+1):
        if a[i]==0: continue
        for j in range(N+1-i):
            if b[j]!=0: c[i+j]+=a[i]*b[j]
    return c

def series_exp(a, N):
    # exp of series a with a[0]=0
    out = [mp.mpf(1)]+[mp.mpf(0)]*N
    for m in range(1, N+1):
        out[m] = sum(j*a[j]*out[m-j] for j in range(1, m+1) if a[j]!=0)/m
    return out

def rhombus_fn(q, N):
    """f_n = [eta^n] I(eta) for the rhombus with parameter q=(A^2-B^2)/(A^2+B^2)."""
    # g(z) = sqrt(cosh z + q sinh z): h(z) = sum a_n z^n, a_n = (n even ? 1 : q)/n!
    a = [(1 if n%2==0 else q)/mp.factorial(n) for n in range(N+1)]
    g = [mp.mpf(0)]*(N+1); g[0]=mp.mpf(1)
    for k in range(1, N+1):
        g[k] = (a[k] - sum(g[i]*g[k-i] for i in range(1,k)))/2
    # logA(eta,s) = sum_{n>=2} (-1)^n zeta(n)/n [ (eta/2+s)^n + (eta/2-s)^n - eta^n ]
    # mu_j(eta) = d_s^j logA(eta,0), even j only:
    #   mu_0 = sum_{n>=2} (-1)^n zeta(n)/n (2^{1-n} - 1) eta^n
    #   mu_j = 2 sum_{n>=j} (-1)^n zeta(n) (n-1)!/(n-j)! (eta/2)^{n-j}   (j>=2 even)
    mu = {}
    mu[0] = [mp.mpf(0)]*(N+1)
    for n in range(2, N+1):
        mu[0][n] = ((-1)**n)*Zv[n]/n*(mp.mpf(2)**(1-n)-1)
    for j in range(2, N+1, 2):
        arr = [mp.mpf(0)]*(N+1)
        for n in range(max(2,j), N+1):
            c = 2*((-1)**n)*Zv[n]*factorial(n-1)/factorial(n-j)/mp.mpf(2)**(n-j)
            arr[n-j] = c
        mu[j] = arr
    # A(eta,s) = exp(mu_0 + sum_{j even>=2} mu_j s^j/j!)
    A0 = series_exp(mu[0], N)
    # d_s^k A(eta,0) for even k: Bell polynomials in mu_2,mu_4,... times A0
    # use series in s: logA_s(s) = sum_{j even} mu_j s^j/j! ; A(s) = A0 * exp(rest)
    dsA = {}
    rest0 = [mp.mpf(0)]*(N+1)
    exprest = {0: [mp.mpf(1)]+[mp.mpf(0)]*N}
    # build A(eta,s) as series in s iteratively
    # exp(rest) where rest = sum_{j>=2 even} mu_j s^j/j!
    # do it via complete Bell: A_k = d_s^k A |0 = A0 * Y_k(mu_2/2!*2!..)
    # simpler: A(s) series in s with coefficients in eta-series; compute exp of rest
    KMAX = N
    exprest_s = [[mp.mpf(1)]+[mp.mpf(0)]*N]  # [s^0]
    # exp(rest): coefficients e_k(eta): e_0=1, e_k = (1/k) sum_{j even<=k} j * (mu_j/j!) * e_{k-j}
    e = [[mp.mpf(1)]+[mp.mpf(0)]*N]
    for k in range(1, KMAX+1):
        ek = [mp.mpf(0)]*(N+1)
        for j in range(2, k+1, 2):
            muj = mu[j]
            ekj = e[k-j]
            for m in range(N+1):
                if muj[m]!=0: ek[m] += j*muj[m]/factorial(j)*ekj[m]*0  # placeholder
        e.append(ek)
    # redo properly: e_k = (1/k) sum_{j even, 2<=j<=k} mu_j/(j-1)! * e_{k-j}  (convolution in eta)
    e = [[mp.mpf(1)]+[mp.mpf(0)]*N]
    for k in range(1, KMAX+1):
        ek = [mp.mpf(0)]*(N+1)
        for j in range(2, k+1, 2):
            conv = series_mul(mu[j], e[k-j], N)
            for m in range(N+1):
                ek[m] += conv[m]/factorial(j-1)
        e.append([v/k for v in ek])
    # d_s^k A(eta,0) = A0 * e_k * 1  (since exp(rest) = sum e_k s^k, and rest has no s^0,s^1)
    # careful: A(eta,s) = A0(eta) * (sum_k e_k(eta) s^k); so d_s^k A = k! A0 e_k
    dsA = {k: [factorial(k)*v for v in series_mul(A0, e[k], N)] for k in range(0, KMAX+1)}
    # I(eta) = sum_{k even} g_k eta^k dsA_k  (odd k vanish since A even in s)
    I = [mp.mpf(0)]*(N+1)
    for k in range(0, N+1):
        if k%2==1: continue
        gk = g[k]
        if gk==0: continue
        ser = dsA[k]
        for m in range(N+1-k):
            I[k+m] += gk*ser[m]
    return I, g

for A, B in [(1, 2), (2, 1)]:
    q = mp.mpf(A*A-B*B)/mp.mpf(A*A+B*B)
    r = mp.sqrt(A*A+B*B)
    L = mp.log(mp.mpf(B)/A)
    rho = mp.sqrt(L**2 + mp.pi**2/4)
    I, g = rhombus_fn(q, N)
    f = I
    # cross-check f_0..f_3 vs direct numerical differentiation of I(eta) by quadrature
    def Idirect(eta):
        return mp.quad(lambda t: (t*(1-t))**(eta/2)*mp.sqrt(mp.cosh(eta*mp.log(t/(1-t)))+q*mp.sinh(eta*mp.log(t/(1-t)))), [0, mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('0.75'), 1])
    h = mp.mpf('1e-20')
    d1 = (Idirect(h)-Idirect(-h))/(2*h)
    print(f"A={A},B={B}: q={mp.nstr(q,6)} rho={mp.nstr(rho,15)}")
    print(f"  f_0={mp.nstr(f[0],10)} f_1={mp.nstr(f[1],10)} (numeric d/d eta: {mp.nstr(d1,10)}) f_2={mp.nstr(f[2],15)} f_3={mp.nstr(f[3],15)}")
    # law: f_n ~ -(sqrt2/(r sqrt(pi))) rho^{1/2-n} n^{-3/2} [B (-1)^n cos((n-1/2)th0) - A sin((n-1/2)th0)] n!
    th0 = mp.atan2(mp.pi/2, L)  # arg of L + i pi/2
    def law(n):
        bracket = B*((-1)**n)*mp.cos((n-mp.mpf('0.5'))*th0) - A*mp.sin((n-mp.mpf('0.5'))*th0)
        return -(mp.sqrt(2)/(r*mp.sqrt(mp.pi)))*rho**(mp.mpf('0.5')-n)*n**(-mp.mpf('1.5'))*bracket*mp.factorial(n)
    print("  n, f_n, law, ratio:")
    for n in [10, 20, 40, 60, 80]:
        lv = law(n)
        print(f"  n={n}: f={mp.nstr(f[n],10)} law={mp.nstr(lv,10)} ratio={mp.nstr(f[n]/lv, 18)}")
    # root test
    bn = [f[n]/mp.factorial(n) for n in range(N+1)]
    for n in [40, 60, 80]:
        print(f"  |b_{n}|^(1/n) = {mp.nstr(abs(bn[n])**(mp.mpf(1)/n), 15)}   1/rho = {mp.nstr(1/rho, 15)}")
