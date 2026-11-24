
import mpmath as mp
mp.mp.dps = 80

def S(phi, m, p):
    return mp.fsum(abs(mp.cos(phi - i*mp.pi/m))**p for i in range(1, m+1))
def rho(phi, m, p):
    return S(phi, m, p)**(-mp.mpf(1)/p)
def drho(phi, m, p):
    return mp.diff(lambda ph: rho(ph, m, p), phi)
def speed(phi, m, p):
    r = rho(phi, m, p); dr = drho(phi, m, p)
    return mp.sqrt(r*r + dr*dr)

def P_of_p(m, p, verbose=False):
    # kink directions (vertex dirs): pi/2 + pi*i/m ; edge-mid dirs: pi*i/m
    ks = [mp.pi/2 + mp.pi*i/m for i in range(2*m)]
    ms = [mp.pi*i/m for i in range(2*m)]
    pts = sorted(set([x % (2*mp.pi) for x in ks+ms]))
    pts = [x for x in pts if mp.mpf('1e-9') < x < 2*mp.pi - mp.mpf('1e-9')]
    b = [mp.mpf(0)] + pts + [2*mp.pi]
    tot = mp.mpf(0)
    for a, bb in zip(b[:-1], b[1:]):
        mid = (a+bb)/2
        seg = mp.quad(lambda ph: speed(ph, m, p), [a, mid, bb])
        if verbose: print(f"  [{mp.nstr(a,6)}, {mp.nstr(bb,6)}] -> {mp.nstr(seg,25)}")
        tot += seg
    return tot

print("m=3 P(1):", mp.nstr(P_of_p(3, mp.mpf(1), True), 30))
print("exact 2 sqrt3:", mp.nstr(2*mp.sqrt(3), 30))
print("m=2 P(1):", mp.nstr(P_of_p(2, mp.mpf(1)), 30), " exact 4 sqrt2:", mp.nstr(4*mp.sqrt(2), 30))
