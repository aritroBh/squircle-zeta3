# 08_Q1.py -- first eps-derivative of the vertex function Q(eps) = -H_3(eps)/Gamma(1+2eps).
#   Q'(0) = -H_3'(0) + 2 gamma Q0,   H_3'(0) = (3/sqrt pi) int_0^oo k(b)^2 k1(b) db,
#   k1(b) = d/deps h_eps(b) at eps=0
#         = -log^2 b/2 - (1-gamma/2) log b + gamma/2 + 3gamma^2/8 + pi^2/16
#           - sum_{n>=1} (-1)^n Gamma(n/2) b^n /(2 n!) [ (n+1) psi(n/2)/2 + 1 + (n+1) log b ]
#   (Mellin-Barnes residues of K1(s) = -(1/2) Gamma(-s/2) Gamma(s) psi(s) (1-s)), and for large b
#   k1(b) ~ sum_{m>=1} (-1)^m (2m-1)! (2m-1) psi(2m) b^(-2m)/m!.
# Check: F'(0) = (2/sqrt pi) int k k1 db + gamma F(0) must equal -2 M1 (squircle, Theorem 1 of main.tex).
import mpmath as mp, sys
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 60
mp.mp.dps = DPS + 15
import importlib
k0 = importlib.import_module("01_kfun_Q0").k

def k1_series(b):
    b = mp.mpf(b)
    with mp.workdps(mp.mp.dps + int(float(b*b)/4/2.3) + 15):
        lb = mp.log(b); g = mp.euler
        s = -lb**2/2 - (1 - g/2)*lb + g/2 + 3*g**2/8 + mp.pi**2/16
        n = 1; tot = 0; small = 0
        tol = mp.mpf(10)**(-(DPS + 25))
        while True:
            t = (-1)**n*mp.gamma(mp.mpf(n)/2)*b**n/(2*mp.factorial(n))*((n+1)*mp.digamma(mp.mpf(n)/2)/2 + 1 + (n+1)*lb)
            tot += t
            if n > 5 and abs(t) < tol:
                small += 1
                if small > 3: break
            else: small = 0
            n += 1
        v = s - tot
    return +v

def k1_asym(b):
    b = mp.mpf(b); s = 0; m = 1; prev = None
    while True:
        t = (-1)**m*mp.factorial(2*m-1)*(2*m-1)*mp.digamma(2*m)/(mp.factorial(m)*b**(2*m))
        if prev is not None and abs(t) > abs(prev): break
        s += t
        if abs(t) < mp.mpf(10)**(-(mp.mp.dps+5)): break
        prev = t; m += 1
    return s

def k1(b):
    return k1_asym(b) if b > 3.1*mp.sqrt(mp.mp.dps) else k1_series(b)

if __name__ == "__main__":
    for b in ['0.3', '2', '9']:
        bb = mp.mpf(b); t = 1/bb**2
        q = mp.quad(lambda x: (mp.log(x)*(-mp.expm1(-t*x*x)) - 2*t*x*x*mp.log(x)*mp.exp(-t*x*x))*mp.exp(-x)/x, [0, bb, mp.inf])
        print("k1(%s) = %s   quad diff %s" % (b, mp.nstr(k1(bb), 25), mp.nstr(k1(bb)-q, 3)))
    bsw = 3.1*mp.sqrt(mp.mp.dps)
    pts = [0, mp.mpf(1)/8, mp.mpf(1)/2, 1, 2, 4, 8, 16, bsw, 2*bsw, mp.inf]
    kk = {}
    def K(b):
        key = mp.nstr(b, mp.mp.dps)
        if key not in kk: kk[key] = (k0(b), k1(b))
        return kk[key]
    I21 = mp.quad(lambda b: K(b)[0]*K(b)[1], pts)
    I31 = mp.quad(lambda b: K(b)[0]**2*K(b)[1], pts)
    F0 = 2*mp.sqrt(2)*mp.log(1+mp.sqrt(2)) - 2*mp.log(2)
    Fp0 = 2*I21/mp.sqrt(mp.pi) + mp.euler*F0
    r2 = mp.sqrt(2); L = mp.log(1+r2); lam = mp.log(2)
    M1 = ((19*r2/48 - mp.mpf(1)/24)*mp.pi**2 + lam**2/2 + 2*r2*lam*L - (1 + 5*r2/2)*L**2 - 5*r2*mp.polylog(2, r2-1))
    print("F'(0) =", mp.nstr(Fp0, DPS), "\n-2M1  =", mp.nstr(-2*M1, DPS), "\nresidual", mp.nstr(Fp0 + 2*M1, 3))
    Q0 = mp.mpf(open("Q0_value.txt").read().strip())
    H3p = 3*I31/mp.sqrt(mp.pi)
    Q1 = -H3p + 2*mp.euler*Q0
    print("Q1 = Q'(0) =", mp.nstr(Q1, DPS))
    open("Q1_value_%d.txt" % DPS, "w").write(mp.nstr(Q1, DPS) + "\n")
