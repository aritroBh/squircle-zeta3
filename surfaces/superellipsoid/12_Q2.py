# 12_Q2.py -- second eps-derivative: Q2 = [eps^2] Q(eps), and s4.
#   k2(b) = d^2/deps^2 h_eps(b) at eps=0, Mellin transform K2(s) = -(1/2) Gamma(-s/2) Gamma''(s) (1-s)^2,
#   small-b series: k2 = Res_0 + sum_{n>=1} (-1)^n/n! g_n''(0),  g_n(z) = -(1/2) Gamma((n-z)/2) (1+n-z)^2 b^(n-z)
#   (Res_0 from residue_k2.py), large b: k2 ~ -sum_{m>=1} (-1)^m/m! Gamma''(2m) (2m-1)^2 b^(-2m).
# Checks: [eps^2] F(eps) = -2 M2 (squircle Theorem 2).
# Output: Q2 and s4 = 48 M3 - 8 pi^2 M1 + 96 zeta(3) M0 - (3/5) pi^4 - 8 Q2   (54 zeta(4) = 3 pi^4/5).
import mpmath as mp, sys, importlib
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 45
mp.mp.dps = DPS + 15
k0 = importlib.import_module("01_kfun_Q0").k
m08 = importlib.import_module("08_Q1"); k1 = m08.k1
mp.mp.dps = DPS + 15

def k2_series(b):
    b = mp.mpf(b)
    with mp.workdps(mp.mp.dps + int(float(b*b)/4/2.3) + 15):
        l = mp.log(b); g = mp.euler; pi2 = mp.pi**2; z3 = mp.zeta(3)
        res0 = (-l**3/3 - 2*l**2 + g*l**2/2 - 2*l - pi2*l/24 - g**2*l/4 + 2*g*l
                - 7*g*pi2/48 - pi2/12 - 7*z3/12 - g**2/2 - 7*g**3/24 + g)
        tot = 0; n = 1; small = 0; tol = mp.mpf(10)**(-(DPS+25))
        while True:
            h = mp.mpf(n)/2
            Gm = mp.gamma(h); p1 = mp.digamma(h); p2 = mp.psi(1, h)
            G0, G1, G2 = Gm, -Gm*p1/2, Gm*(p1**2 + p2)/4
            P0, P1, P2 = (1+n)**2, -2*(1+n), 2
            E0 = b**n; E1 = -l*E0; E2 = l*l*E0
            d2 = G2*P0*E0 + G0*P2*E0 + G0*P0*E2 + 2*G1*P1*E0 + 2*G1*P0*E1 + 2*G0*P1*E1
            t = (-1)**n/mp.factorial(n)*(-d2/2)
            tot += t
            if n > 6 and abs(t) < tol:
                small += 1
                if small > 3: break
            else: small = 0
            n += 1
        v = res0 + tot
    return +v

def k2_asym(b):
    b = mp.mpf(b); s = 0; m = 1; prev = None
    while True:
        G2m = mp.factorial(2*m-1)*(mp.digamma(2*m)**2 + mp.psi(1, 2*m))
        t = -(-1)**m/mp.factorial(m)*G2m*(2*m-1)**2/b**(2*m)
        if prev is not None and abs(t) > abs(prev): break
        s += t
        if abs(t) < mp.mpf(10)**(-(mp.mp.dps+5)): break
        prev = t; m += 1
    return s

def k2(b):
    return k2_asym(b) if b > 3.1*mp.sqrt(mp.mp.dps) else k2_series(b)

if __name__ == "__main__":
    for b in ['0.3', '2', '9']:
        bb = mp.mpf(b); t = 1/bb**2
        # d^2/deps^2 [(1-exp(-t x^(2-2e))) x^(e-1)] at e=0
        def integrand(x):
            lx = mp.log(x); E = mp.exp(-t*x*x); X2 = x*x
            # f(e) = (1-exp(-t x^{2-2e})) x^{e-1};  f'' = x^{-1}[ lx^2 (1-E) + 2 lx * (-2 t X2 lx E) + (4 t X2 lx^2 E - 4 t^2 X2^2 lx^2 E) ]
            return (lx**2*(1-E) - 4*t*X2*lx**2*E + 4*t*X2*lx**2*E - 4*t*t*X2*X2*lx**2*E)*mp.exp(-x)/x
        q = mp.quad(integrand, [0, bb, mp.inf])
        print("k2(%s) = %s   quad diff %s" % (b, mp.nstr(k2(bb), 25), mp.nstr(k2(bb)-q, 3)))
    bsw = 3.1*mp.sqrt(mp.mp.dps)
    pts = [0, mp.mpf(1)/8, mp.mpf(1)/2, 1, 2, 4, 8, 16, bsw, 2*bsw, mp.inf]
    cache = {}
    def K(b):
        key = mp.nstr(b, mp.mp.dps)
        if key not in cache: cache[key] = (k0(b), k1(b), k2(b))
        return cache[key]
    sp_ = mp.sqrt(mp.pi)
    H2pp = mp.quad(lambda b: 2*K(b)[1]**2 + 2*K(b)[0]*K(b)[2], pts)/sp_
    H3pp = mp.quad(lambda b: 6*K(b)[0]*K(b)[1]**2 + 3*K(b)[0]**2*K(b)[2], pts)/sp_
    H2p = mp.quad(lambda b: 2*K(b)[0]*K(b)[1], pts)/sp_
    H3p = mp.quad(lambda b: 3*K(b)[0]**2*K(b)[1], pts)/sp_
    H2 = mp.quad(lambda b: K(b)[0]**2, pts)/sp_
    H3 = mp.quad(lambda b: K(b)[0]**3, pts)/sp_
    g = mp.euler; pi2 = mp.pi**2
    F2c = H2pp/2 + g*H2p + (g**2/2 - pi2/12)*H2
    import squircle_consts
    M0, M1, M2 = squircle_consts.consts()
    print("[eps^2]F =", mp.nstr(F2c, DPS), "\n   -2M2  =", mp.nstr(-2*M2, DPS), "\n residual", mp.nstr(F2c + 2*M2, 3))
    Q0 = -H3; Q1 = -(H3p + 2*g*H3)
    Q2 = -(H3pp/2 + 2*g*H3p + (2*g**2 - pi2/3)*H3)
    print("Q0 =", mp.nstr(Q0, DPS)); print("Q1 =", mp.nstr(Q1, DPS)); print("Q2 =", mp.nstr(Q2, DPS))
    M3 = mp.mpf('0.02277052149178657745673050974232780025018369172921726137909')
    s4 = 48*M3 - 8*pi2*M1 + 96*mp.zeta(3)*M0 - 3*mp.pi**4/5 - 8*Q2
    print("s4 =", mp.nstr(s4, 50))
    open("Q2_value_%d.txt" % DPS, "w").write("Q0 %s\nQ1 %s\nQ2 %s\ns4 %s\n" % tuple(mp.nstr(x, DPS) for x in (Q0, Q1, Q2, s4)))
