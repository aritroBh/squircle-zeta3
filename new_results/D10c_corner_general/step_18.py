"""step_18.py -- expanded closed form M1(alpha) valid for all alpha in (0,pi):
weight matrix W(rho,sig), Li2 argument tables (functions of c,s), verification.
Convention: C(rho,sig) = log|rho-sig| - i*pi if rho<sig else log(rho-sig)  [D10b].
"""
import mpmath as mp
from corner_series import corner_coeffs
mp.mp.dps = 90
log2 = mp.log(2); z2 = mp.pi**2/6; ipi = mp.pi*1j

def build(a):
    c = mp.cos(a/2); s = mp.sin(a/2); ca = mp.cos(a)
    r1, r2, t1 = 2*c - 1, -1 - 2*c, 1 + 2*s
    poles = [mp.mpf(1), mp.mpf(-1), ca, r1, r2]
    ares = [mp.mpf(0), mp.mpf(-2), mp.mpf(1), -c, c]
    bres = [mp.mpf(-1), mp.mpf(1), mp.mpf(-1), 2*c, -2*c]
    # weight matrix W[rho][sig], sig in [1,-1,ca,r1,r2]
    W = {}
    for rho, ar, br in zip(poles, ares, bres):
        W[(rho, mp.mpf(1))] = ar; W[(rho, mp.mpf(-1))] = ar
        W[(rho, ca)] = -(ar + br); W[(rho, r1)] = br; W[(rho, r2)] = br
    logw = {rho: -(ar + br) for rho, ar, br in zip(poles, ares, bres)}
    return c, s, ca, r1, r2, t1, poles, W, logw

def Cfun(rho, sig):
    return mp.log(rho - sig) if rho > sig else mp.log(sig - rho) - ipi

def M1_expanded(a):
    c, s, ca, r1, r2, t1, poles, W, logw = build(a)
    tot = mp.mpc(0); nterms = 0; zeta2_seen = mp.mpc(0)
    for (rho, sig), w in W.items():
        if w == 0: continue
        if rho == sig:
            tot += w*(mp.log(t1 - rho)**2 - mp.log(1 - rho)**2)/2
            continue
        for tt, sgn in [(t1, 1), (mp.mpf(1), -1)]:
            x = (tt - rho)/(sig - rho)
            if abs(x) < mp.mpf('1e-80'):
                continue
            if abs(x - 1) < mp.mpf('1e-80'):
                tot += sgn*w*(-z2); zeta2_seen += -sgn*w
                continue
            tot += sgn*w*(-mp.polylog(2, x))
            nterms += 1
        tot += w*Cfun(rho, sig)*(mp.log(t1 - rho) - mp.log(1 - rho))
    # regularized rho=1 row: lower limit contributes 0; but the -log2 term at rho=1
    # needs the row limit too: handle -log2 terms per row with regularization at rho=1
    for rho, lw in logw.items():
        if lw == 0: continue
        if rho == 1:
            tot += lw*log2*mp.log(t1 - rho)   # regularized: log(1-rho) dropped (row limit 0)
        else:
            tot += lw*log2*(mp.log(t1 - rho) - mp.log(1 - rho))
    return tot, nterms, zeta2_seen

print("expanded formula vs corner_series (60+ digits target):")
for nm, a in [("2pi/3", 2*mp.pi/3), ("pi/2", mp.pi/2), ("pi/3", mp.pi/3), ("3pi/4", 3*mp.pi/4),
              ("pi/6", mp.pi/6), ("5pi/6", 5*mp.pi/6), ("pi/12", mp.pi/12), ("1.3", mp.mpf('1.3'))]:
    v, nterm, zseen = M1_expanded(a)
    ref = corner_coeffs(a, 1)[1]
    res = v - ref
    dig = float(-mp.log10(abs(res))) if res != 0 else 999.0
    print(f"  {nm:>6}: resid {mp.nstr(res,2)} ({dig:.1f} dig), Im {mp.nstr(abs(mp.im(v)),2)}, Li2 terms {nterm}, zeta2-weight {mp.nstr(zseen,3)}")

# Li2 argument table (symbolic, functions of c,s): upper endpoint t1 = 1+2s
import sympy as sp
cs, ss = sp.symbols('c s')
cas = 2*cs**2 - 1; r1s, r2s, t1s = 2*cs - 1, -1 - 2*cs, 1 + 2*ss
poless = [sp.Integer(1), sp.Integer(-1), cas, r1s, r2s]
aress = [0, -2, 1, -cs, cs]; bress = [-1, 1, -1, 2*cs, -2*cs]
print("\nsymbolic Li2 arguments at t1 (rows rho, cols sig) with weights W:")
sigl = [sp.Integer(1), sp.Integer(-1), cas, r1s, r2s]
for rho, ar, br in zip(poless, aress, bress):
    Wrow = [ar, ar, -(ar + br), br, br]
    for sig, w in zip(sigl, Wrow):
        if w == 0 or rho == sig: continue
        x1 = sp.simplify((t1s - rho)/(sig - rho))
        print(f"  W={str(w):>4}  Li2({sp.simplify(x1)})   [rho={rho}, sig={sig}]")
