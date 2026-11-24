"""report_formula_check.py -- re-types the formulas exactly as written in result.md and checks them against the
stored series values (318 digits) at all 9 angles.  Guards against transcription errors in the report.
 (F1) M1 classical (Thm B1)      (F2) M2 classical uniform form (Thm B2')    (F3) M2(2pi/3) minimal form (basis A)
 (F4) M2 = z2/2 - Li2(c^2)/2 - z2 (1-c^2)/(2c) artanh c + S(c) + S(-c)  (HPL form, Thm B2), 60 digits.
"""
import mpmath as mp
import hpl_eval as he

mp.mp.dps = 120
Li2 = lambda z: mp.polylog(2, z)
Li3 = lambda z: mp.polylog(3, z)
z2, z3, pi2, lam = mp.zeta(2), mp.zeta(3), mp.pi**2, mp.log(2)


def M0(c):
    return -(1 + c)*mp.log(1 + c) - (1 - c)*mp.log(1 - c)


def M1(c):  # (F1)
    return z2 - lam**2 - lam*M0(c) - sum((1 - e*c)*(Li2(e*c) + Li2((1 - e*c)/2)) for e in (1, -1))


def Phi(c):
    return mp.mpf(1)/2 - (1 - c*c)/(2*c)*mp.atanh(c)


def M2(c):  # (F2)
    lp, lm, l0 = mp.log(1 + c), mp.log(1 - c), mp.log(c)
    q = (1 - c*c)/(2*c)
    v = z2*Phi(c) - Li2(c) - Li2(-c)
    v += ((1 + c)*Li3(c) + (1 - c)*Li3(-c) + (1 - c)/c*Li3(1 - c) - (1 + c)/c*Li3(1/(1 + c))
          - 2*c*Li3((1 + c)/2) + (1 - c)*(Li3((1 - c)/(1 + c)) + Li3(-(1 - c)/(1 + c))) - 2*c*Li3(2*c/(1 + c)))
    v += Li2((1 + c)/2)*((1 + c)*lp - (1 - c)*lm) + q*(Li2(c) - Li2(-c))*(lm + lp)
    v += (-(1 - c)*lm**2*lp + (1 - c)*lm*lp**2 + (4*c*c - c + 1)/(6*c)*lp**3
          + q*l0*lm**2 - (1 - c)*l0*lm*lp - (1 + c*c)/(2*c)*l0*lp**2)
    v += lam*((1 - c)*lm**2 - (1 - c)*lm*lp - c*lp**2) + lam**2*(1 - c)/2*(lp - lm)
    v += -pi2*((1 - c)/(6*c)*lm + (1 + c)/(6*c)*lp)
    v += (mp.mpf(7)/4 + 2*c)*z3 - c/6*pi2*lam + c/3*lam**3
    return v


def M2_2pi3():  # (F3)
    l3 = mp.log(3); d = Li2(mp.mpf(1)/3)
    W2 = d + lam**2 - lam*l3 + l3**2/2
    W3 = (mp.mpf(15)/2*Li3(mp.mpf(1)/3) + mp.mpf(17)/2*Li3(-mp.mpf(1)/2) + mp.mpf(3)/16*z3
          - mp.mpf(5)/12*pi2*lam + pi2*l3/4 - lam**3/12 + mp.mpf(5)/4*lam**2*l3 - mp.mpf(7)/4*lam*l3**2
          + mp.mpf(5)/8*l3**3 - mp.mpf(7)/2*d*lam + mp.mpf(15)/4*d*l3)
    return W2 + W3


def S(c):   # HPL half of T3 (F4); S(-c) uses H_w(-c) = (-1)^{#nonzero letters} H_{-w}(c)
    H = he.H
    return ((1 - c)*(H((0, 0, 1), c) + H((0, 1, -1), c) - H((1, -1, 1), c))
            + (1 - c*c)/(2*c)*H((1, 0, 1), c) - (1 - c)**2/(2*c)*H((1, 0, -1), c))


def Sneg(c):
    H = he.H
    return ((1 + c)*(-H((0, 0, -1), c) + H((0, -1, 1), c) + H((-1, 1, -1), c))
            - (1 - c*c)/(2*c)*H((-1, 0, -1), c) + (1 + c)**2/(2*c)*H((-1, 0, 1), c))


ANG = {'pi_2': mp.pi/2, '2pi_3': 2*mp.pi/3, 'pi_3': mp.pi/3, '3pi_4': 3*mp.pi/4, 'pi_6': mp.pi/6,
       '5pi_6': 5*mp.pi/6, 'pi_4': mp.pi/4, 'gen1p3': mp.mpf(13)/10, 'gen2p0': mp.mpf(2)}
for tag, a in ANG.items():
    c = mp.cos(a/2)
    r1, r2 = he.stored(tag, 1), he.stored(tag, 2)
    print(f'{tag:>7}: |F1 - M1| = {mp.nstr(abs(M1(c) - r1), 3)}  |F2 - M2| = {mp.nstr(abs(M2(c) - r2), 3)}', flush=True)
print('2pi/3 : |F3 - M2| =', mp.nstr(abs(M2_2pi3() - he.stored('2pi_3', 2)), 3))
mp.mp.dps = 60
for tag in ['pi_2', '2pi_3', 'pi_6']:
    c = mp.cos(ANG[tag]/2)
    he._cache.clear()
    F4 = z2/2 - Li2(c*c)/2 - z2*(1 - c*c)/(2*c)*mp.atanh(c) + S(c) + Sneg(c)
    print(f'{tag:>7}: |F4 (HPL form) - M2| = {mp.nstr(abs(F4 - he.stored(tag, 2)), 3)}', flush=True)
