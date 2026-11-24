# 09_Q0_direct.py -- Q0 by two routes that do not use J, P1, P2a, P2c:
#  (a) Q0 = -(1/sqrt(pi)) int_0^oo k0(y)^3 dy,  k0(y) = int_0^oo (1-e^{-x^2/y^2}) e^{-x} dx/x
#      with k0(y) = -log y - gamma/2 + (pi/2) erfi(y/2) - (y^2/4) 2F2(1,1;3/2,2;y^2/4)  (checked vs quadrature)
#      and the asymptotic series for large y.
#  (b) the simplex form Q0 = int_Delta [sum_{i<j} sqrt(s_i^2+s_j^2) - 1 - rho]/(s1 s2 s3) dsigma by 2-D quadrature
#      (6-fold symmetry, polar coordinates about a vertex).
import mpmath as mp, sys, json
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 75
mode = sys.argv[2] if len(sys.argv) > 2 else "k"
mp.mp.dps = DPS + 20

def k_closed(y):
    return -mp.log(y) - mp.euler/2 + mp.pi/2*mp.erfi(y/2) - y**2/4*mp.hyp2f2(1, 1, 1.5, 2, y**2/4)
def k_asym(y):
    s = mp.mpf(0); n = 1; prev = None
    while True:
        t = (-1)**(n + 1)*mp.factorial(2*n - 1)/(mp.factorial(n)*y**(2*n))
        if prev is not None and abs(t) > abs(prev): raise ValueError("asymptotic series not accurate")
        s += t
        if abs(t) < mp.mpf(10)**(-(mp.mp.dps + 5)): return s
        prev = t; n += 1
def k0(y):
    y = mp.mpf(y)
    if y > 3.2*mp.sqrt(mp.mp.dps): return k_asym(y)
    with mp.workdps(mp.mp.dps + int(y*y/4/2.3) + 10):
        return +k_closed(y)

if mode == "k":
    for yy in ['0.1', '1', '3.7', '12']:
        kq = mp.quad(lambda x: (1 - mp.exp(-x*x))*mp.exp(-mp.mpf(yy)*x)/x, [0, 1, 10, mp.inf])
        print("k0(%s): closed vs quadrature %s" % (yy, mp.nstr(abs(k0(yy) - kq), 3)))
    bsw = 3.2*mp.sqrt(mp.mp.dps)
    pts = [0, mp.mpf(1)/16, mp.mpf(1)/4, 1, 2, 4, 8, 16, bsw, 2*bsw, 4*bsw, mp.inf]
    I1 = mp.quad(k0, pts)/mp.sqrt(mp.pi)
    I2 = mp.quad(lambda y: k0(y)**2, pts)/mp.sqrt(mp.pi)
    I3 = mp.quad(lambda y: k0(y)**3, pts)/mp.sqrt(mp.pi)
    print("(1/sqrt pi) int k0   - 1 =", mp.nstr(I1 - 1, 3))
    print("(1/sqrt pi) int k0^2 - (2 sqrt2 L - 2 log2) =", mp.nstr(I2 - (2*mp.sqrt(2)*mp.log(1 + mp.sqrt(2)) - 2*mp.log(2)), 3))
    print("Q0 (k0^3) =", mp.nstr(-I3, DPS))
    json.dump({"Q0_k3": mp.nstr(-I3, DPS + 5)}, open("Q0_k3.json", "w"))
else:
    # simplex: region u >= v >= w (one sixth), vertex polar coordinates about (1,0,0):
    # v = r t, w = r (1-t), u = 1 - r, t in [1/2,1], r in [0, 1/(1+t)]; dsigma = r dr dt
    def R0(u, v, w):
        return (mp.sqrt(u*u + v*v) + mp.sqrt(v*v + w*w) + mp.sqrt(w*w + u*u) - mp.sqrt(u*u + v*v + w*w) - 1)/(u*v*w)
    def inner(t):
        def f(r):
            # cancellation near r=0: numerator ~ r^3; raise precision there
            extra = int(max(0, -3*mp.log10(r))) + 5 if r > 0 else 0
            with mp.workdps(mp.mp.dps + extra):
                v = r*R0(1 - r, r*t, r*(1 - t))
            return +v
        return mp.quad(f, [0, mp.mpf(1)/(4*(1 + t)), 1/(1 + t)])
    val = 6*mp.quad(inner, [mp.mpf(1)/2, mp.mpf(3)/4, 1])
    print("Q0 (simplex, 2-D quadrature) =", mp.nstr(val, DPS))
    json.dump({"Q0_simplex": mp.nstr(val, DPS + 5)}, open("Q0_simplex_%d.json" % DPS, "w"))
