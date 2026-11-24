
from mpmath import mp, mpf, mpc, pi, sqrt, exp, log, nstr, gamma, quad
mp.dps = 80

def load_mk(fn):
    out = {}
    for line in open(fn):
        p = line.split("\t")
        if len(p) < 2: continue
        out[int(p[0])] = mpf(p[1].strip().split(" ")[0].strip("[]"))
    return out
Mk = load_mk("../D13b_Mk_constant/mk_deg9.tsv")

def H(s):
    return sqrt(pi)/4 * gamma(s/2)/gamma(s/2+mpf('1.5'))

def Mk_main(k, dps=60):
    # (2/pi) Im int_0^inf e^{i pi s/2} (s/(1+s))^N H(s) ds/s, contour s = e^{i pi/4} x
    old = mp.dps; mp.dps = dps
    N = mpf(k+1)
    rot = exp(mpc(0,pi/4))
    def integrand(x):
        s = rot*x
        return rot * exp(mpc(0,1)*pi*s/2) * (s/(1+s))**N * H(s)/s
    try:
        v = quad(integrand, [0, mpf('0.5'), 2, 10, 50, mp.inf])
    finally:
        mp.dps = old
    return (2/pi)*v.imag

print("k, M_k (table), M_k^main (1D integral), difference, e^{-sqrt(2 pi N)} scale:")
for k in [50, 100, 150, 200, 250, 300, 400]:
    mm = Mk_main(k)
    N = mpf(k+1)
    dif = Mk[k]-mm
    scale = exp(-sqrt(2*pi*N))
    print("  k=%3d: Mk=%s  main=%s  diff=%s  (diff/scale=%s)" % (
        k, nstr(Mk[k],8), nstr(mm,8), nstr(dif,4), nstr(dif/scale,4)))
