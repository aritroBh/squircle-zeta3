# 15_Q1_highprec.py -- Q1 = 2 Q0 - 3 Zq + 6 Za + 6 Zb to high precision, using the closed form of
# Psi1 (psi1_closed.py) and the proved closed form of Q0 (13_proof_Q0.py).
import mpmath as mp, sys, json
from psi1_closed import Psi1_closed as Psi1
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 100
mp.mp.dps = DPS + 30
def Q0_closed():
    r2 = mp.sqrt(2); lam, l3, L = mp.log(2), mp.log(3), mp.log(1+r2); pi = mp.pi; b = r2-1
    al = mp.atan(1/r2)
    return (-pi**2/8 + r2*pi**2/16 + 3*r2/2*pi*al - 3*lam**2 - 9*r2*l3*L + 6*L**2 + 6*r2*L**2
            - 6*r2*mp.polylog(2, b) + 9*r2/4*mp.polylog(2, -b**4))
def extra(x):
    return int(max(0, -mp.log10(abs(x)))) + 10 if x != 0 else 0
P1v = Psi1(1)
def fq(s):
    s = mp.mpf(s); m = s*(1-s)
    with mp.workdps(mp.mp.dps + 2*extra(m)):
        q = mp.sqrt(1-2*m); v = (Psi1(q) - P1v)/m
    return +v
def fa(t):
    t = mp.mpf(t)
    with mp.workdps(mp.mp.dps + 2*extra(1-t)):
        v = (Psi1(t) - P1v)/(1-t)
    return +v
def fb(t):
    t = mp.mpf(t)
    with mp.workdps(mp.mp.dps + 2*extra(t)):
        v = (Psi1(t) + 1)/t
    return +v
Zq = 2*mp.quad(fq, [0, mp.mpf(1)/16, mp.mpf(1)/4, mp.mpf(1)/2])
Za = mp.quad(fa, [0, mp.mpf(1)/4, mp.mpf(1)/2, mp.mpf(3)/4, 1])
Zb = mp.quad(fb, [0, mp.mpf(1)/256, mp.mpf(1)/16, mp.mpf(1)/4, 1])
Q0 = Q0_closed()
Q1 = 2*Q0 - 3*Zq + 6*Za + 6*Zb
print("Q0 =", mp.nstr(Q0, DPS)); print("Zq =", mp.nstr(Zq, DPS)); print("Za =", mp.nstr(Za, DPS)); print("Zb =", mp.nstr(Zb, DPS))
print("Q1 =", mp.nstr(Q1, DPS))
for fn in ("Q1_value_105.txt", "Q1_value_50.txt"):
    try:
        ref = mp.mpf(open(fn).read().strip()); print("  vs k-route %s: %s" % (fn, mp.nstr(Q1-ref, 3)))
    except Exception: pass
json.dump({k: mp.nstr(v, DPS) for k, v in dict(Q0=Q0, Zq=Zq, Za=Za, Zb=Zb, Q1=Q1, Psi1_1=P1v).items()},
          open("Q1_hp_%d.json" % DPS, "w"), indent=1)
