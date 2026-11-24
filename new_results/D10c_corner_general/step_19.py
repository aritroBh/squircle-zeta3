"""step_19.py -- expanded formula, fixed rho=1 row regularization (row lower limit = 0)."""
import mpmath as mp
from corner_series import corner_coeffs
mp.mp.dps = 90
log2 = mp.log(2); z2 = mp.pi**2/6; ipi = mp.pi*1j

def M1_expanded(a):
    c = mp.cos(a/2); s = mp.sin(a/2); ca = mp.cos(a)
    r1, r2, t1 = 2*c - 1, -1 - 2*c, 1 + 2*s
    poles = [mp.mpf(1), mp.mpf(-1), ca, r1, r2]
    ares = [mp.mpf(0), mp.mpf(-2), mp.mpf(1), -c, c]
    bres = [mp.mpf(-1), mp.mpf(1), mp.mpf(-1), 2*c, -2*c]
    W = {}
    for rho, ar, br in zip(poles, ares, bres):
        W[(rho, mp.mpf(1))] = ar; W[(rho, mp.mpf(-1))] = ar
        W[(rho, ca)] = -(ar + br); W[(rho, r1)] = br; W[(rho, r2)] = br
    def Cfun(rho, sig):
        return mp.log(rho - sig) if rho > sig else mp.log(sig - rho) - ipi
    tot = mp.mpc(0); nterms = 0; zseen = mp.mpc(0)
    for (rho, sig), w in W.items():
        if w == 0: continue
        if rho == sig:
            tot += w*(mp.log(t1 - rho)**2 - (mp.log(1 - rho)**2 if rho != 1 else 0))/2
            continue
        endpoints = [(t1, 1)] if rho == 1 else [(t1, 1), (mp.mpf(1), -1)]
        for tt, sgn in endpoints:
            x = (tt - rho)/(sig - rho)
            if abs(x) < mp.mpf('1e-80'): continue
            if abs(x - 1) < mp.mpf('1e-80'):
                tot += sgn*w*(-z2); zseen += -sgn*w; continue
            tot += sgn*w*(-mp.polylog(2, x)); nterms += 1
        if rho == 1:
            tot += w*Cfun(rho, sig)*mp.log(t1 - rho)
        else:
            tot += w*Cfun(rho, sig)*(mp.log(t1 - rho) - mp.log(1 - rho))
    logw = {rho: -(ar + br) for rho, ar, br in zip(poles, ares, bres)}
    for rho, lw in logw.items():
        if lw == 0: continue
        if rho == 1: tot += lw*log2*mp.log(t1 - rho)
        else: tot += lw*log2*(mp.log(t1 - rho) - mp.log(1 - rho))
    return tot, nterms, zseen

print("expanded formula vs corner_series:")
for nm, a in [("2pi/3", 2*mp.pi/3), ("pi/2", mp.pi/2), ("pi/3", mp.pi/3), ("3pi/4", 3*mp.pi/4),
              ("pi/6", mp.pi/6), ("5pi/6", 5*mp.pi/6), ("pi/12", mp.pi/12), ("1.3", mp.mpf('1.3'))]:
    v, nterm, zseen = M1_expanded(a)
    ref = corner_coeffs(a, 1)[1]
    res = mp.re(v) - ref
    dig = float(-mp.log10(abs(res))) if res != 0 else 999.0
    print(f"  {nm:>6}: resid {mp.nstr(res,2)} ({dig:.1f} dig), Im {mp.nstr(abs(mp.im(v)),2)}, Li2 terms {nterm}, zeta2 weight {mp.nstr(zseen,3)}")
