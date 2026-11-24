# 01_kfun_Q0.py
# k(b) = int_0^oo (1-exp(-y^2)) exp(-b y) dy / y
#      = -log b - gamma/2 + (pi/2) erfi(b/2) - (b^2/4) 2F2(1,1;3/2,2;b^2/4)
# Checks: (1) closed form vs quadrature, (2) F(0)=(1/sqrt pi) int k^2 = 2 sqrt2 L - 2 log2,
#         (3) Q0 = -(1/sqrt pi) int k^3  vs direct 2D simplex quadrature of R_0.
import mpmath as mp, sys
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 50
mp.mp.dps = DPS + 20

def k_closed(b):
    b = mp.mpf(b)
    return -mp.log(b) - mp.euler/2 + mp.pi/2*mp.erfi(b/2) - b**2/4*mp.hyp2f2(1, 1, 1.5, 2, b**2/4)

def k_asym(b, N=None):
    # k(b) ~ sum_{n>=1} (-1)^(n+1) (2n-1)!/(n! b^(2n))  (asymptotic, b large)
    b = mp.mpf(b); s = mp.mpf(0); n = 1
    term_prev = None
    while True:
        t = (-1)**(n+1)*mp.factorial(2*n-1)/(mp.factorial(n)*b**(2*n))
        if term_prev is not None and abs(t) > abs(term_prev): break
        s += t
        if abs(t) < mp.mpf(10)**(-(mp.mp.dps+5)): break
        term_prev = t; n += 1
    return s

def k(b):
    b = mp.mpf(b)
    if b > 3.1*mp.sqrt(mp.mp.dps):  # asymptotic series is accurate to ~exp(-b^2/4)
        return k_asym(b)
    with mp.workdps(mp.mp.dps + int(b*b/4/2.3) + 10):
        return +k_closed(b)

if __name__ == "__main__":
    # (1) closed form vs quadrature
    for b in ['0.1', '1', '3.7', '12']:
        kq = mp.quad(lambda y: (1-mp.exp(-y*y))*mp.exp(-mp.mpf(b)*y)/y, [0, 1, 10, mp.inf])
        print("b=%s  k=%s  |k-quad|=%s" % (b, mp.nstr(k(b), 25), mp.nstr(abs(k(b)-kq), 3)))
    bsw = 3.1*mp.sqrt(mp.mp.dps)
    print("switch point", mp.nstr(bsw, 5), " closed vs asym at switch:",
          mp.nstr(abs(k_asym(bsw*1.0001) - k_closed(bsw*1.0001)), 3))
    pts = [0, mp.mpf(1)/4, 1, 2, 4, 8, bsw, 2*bsw, mp.inf]
    F0 = mp.quad(lambda b: k(b)**2, pts)/mp.sqrt(mp.pi)
    L = mp.log(1+mp.sqrt(2)); lam = mp.log(2)
    F0ex = 2*mp.sqrt(2)*L - 2*lam
    print("F(0)  =", mp.nstr(F0, DPS))
    print("2sqrt2 L - 2log2 =", mp.nstr(F0ex, DPS), " residual", mp.nstr(F0-F0ex, 3))
    Q0 = -mp.quad(lambda b: k(b)**3, pts)/mp.sqrt(mp.pi)
    print("Q0 =", mp.nstr(Q0, DPS))
    with open("Q0_value.txt", "w") as f:
        f.write(mp.nstr(Q0, DPS) + "\n")
