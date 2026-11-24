"""Shared numerics for RSTS. M_k table copied from D13f_sweep/mk_deg9.tsv (arb balls, ~80 digits)."""
import mpmath as mp, re, os
pi = mp.pi
HERE = os.path.dirname(os.path.abspath(__file__))
def H(s):  # H(s) = (1/2) B(s/2, 3/2)
    return mp.beta(s/2, mp.mpf(3)/2)/2
def load_M():
    out = {}
    for line in open(os.path.join(HERE, "mk_deg9.tsv")):
        p = line.rstrip("\n").split("\t")
        if p and p[0].isdigit():
            out[int(p[0])] = mp.mpf(re.match(r"\[([-0-9.eE+]+)", p[1]).group(1))
    return out
def ray(f, scale):
    om = mp.expjpi(mp.mpf(1)/4)
    g = lambda r: f(r*om)*om
    pts = [0, scale/8, scale/2, scale, 2*scale, 4*scale, 10*scale, 40*scale, mp.inf]
    return mp.quad(g, pts)
def Mmain(k):
    """M_k^main = (2/pi) Im int_0^inf H(s) e^{i pi s/2} s^{-1} (s/(1+s))^N ds  (CX6 eq. after (7)); ray arg=pi/4."""
    N = k+1
    f = lambda s: H(s)*mp.expj(pi*s/2)*(s/(1+s))**N/s
    return 2/pi*ray(f, mp.sqrt(2*N/pi)).imag
def S2(k):
    """S_k = (2/pi) Im int_0^inf H(s) e^{i pi s/2} A2_k(s) ds, A2_k = -2 (s+1) s^{k-1}/(s+2)^{k+1}."""
    N = k+1
    f = lambda s: H(s)*mp.expj(pi*s/2)*(-2)*(s+1)*s**(k-1)/(s+2)**(k+1)
    return 2/pi*ray(f, 2*mp.sqrt(N/pi)).imag
def lead2(k):
    """Predicted leading term (sqrt(pi)/N) e^{-th2} cos(th2 - pi/4), th2 = sqrt(2 pi N)."""
    N = k+1; th = mp.sqrt(2*pi*N)
    return mp.sqrt(pi)/N*mp.exp(-th)*mp.cos(th-pi/4)
def ray_fine(f, rs, width):
    """Integrate f along arg s = pi/4 with dense breakpoints of spacing `width` around |s| = rs."""
    om = mp.expjpi(mp.mpf(1)/4)
    g = lambda r: f(r*om)*om
    pts = [mp.mpf(0)] + [rs*x for x in (mp.mpf(1)/16, mp.mpf(1)/8, mp.mpf(1)/4, mp.mpf(1)/2)]
    j = -int(rs/(2*width))
    while True:
        p = rs + j*width
        if p > rs/2: pts.append(p)
        if j > 0 and p > 3*rs: break
        j += 1
    pts += [4*rs, 8*rs, 32*rs, mp.inf]
    return mp.quad(g, pts)
def Z2(k):
    N = k+1; lam = mp.sqrt(N/pi); Lam = mp.sqrt(pi*N)
    f = lambda s: H(s)*mp.expj(pi*s/2)*(-2)*(s+1)*s**(k-1)/(s+2)**(k+1)
    return ray_fine(f, 2*lam, 2*lam/mp.sqrt(Lam)/2)
