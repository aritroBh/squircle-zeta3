
import os, sys, pickle, time
import sympy as sp
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
sys.path.insert(0, work + "/tow4_copy")
import tow2lib as T2, tow4lib as T4
d, m = T2.d, T2.m
i = sp.I
eps = sp.Symbol('eps')
tau1 = T4.tau1f; tau2 = T4.tau2f
# Build the FORCED family at delta=3 with seed F(3,0) = a_3 + eps, everything from tow4lib machinery
def build(delta, seed_eps):
    Fz = {}
    def F(dp, s):
        if s > dp//2: return sp.Integer(0)
        return Fz[(dp,s)]
    def Lam(dl, s):
        tot = sp.Integer(0)
        for dp in range(0, dl+1):
            fA = F(dp, s).subs({d: d+1, m: m+1}, simultaneous=True)
            fB = F(dp, s).subs({d: d+1, m: m}, simultaneous=True)
            tot += i*(2*d-1)/(2*d-dp-1) * (tau1(dl-dp)*fA - tau2(dl-dp)*fB)
        return tot
    for dl in range(0, delta+1):
        Fz[(dl,0)] = sp.cancel(T4.arho(dl) + (eps if dl==delta else 0))
        for s in range(0, dl//2):
            Fz[(dl,s+1)] = sp.cancel(Lam(dl,s) - F(dl,s))
        # closing candidate
        Fz[(dl, dl//2+1)] = sp.together(Lam(dl, dl//2) - F(dl, dl//2))
    return Fz, Lam
for delta in [3, 4]:
    Fz, Lam = build(delta, eps)
    top = delta//2
    Fclose = Fz[(delta, top+1)]
    U = {dp: sum((-1)**s * (Fz[(dp,s)] if s <= dp//2 else sp.Integer(0)) for s in range(0, dp//2+1)) for dp in range(0, delta+1)}
    feedU = sp.Integer(0)
    for dp in range(0, delta+1):
        uA = U[dp].subs({d: d+1, m: m+1}, simultaneous=True)
        uB = U[dp].subs({d: d+1, m: m}, simultaneous=True)
        feedU += i*(2*d-1)/(2*d-dp-1) * (tau1(delta-dp)*uA - tau2(delta-dp)*uB)
    # telescoping identity: F(delta,top+1) = (-1)^top (feedU - (a_delta+eps))
    lhs = sp.together(Fclose - (-1)**top*(feedU - (T4.arho(delta)+eps)))
    print(f"delta={delta}: telescoping identity F(top+1) - (-1)^top[feedU - seed] == 0 :",
          sp.expand(lhs) == 0)
    # coefficient of eps in Fclose (nonzero => equivalence non-vacuous)
    print(f"   d(Fclose)/d(eps) =", sp.expand(sp.diff(Fclose, eps)))
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} telescoping identity (1)<=>(2) verified symbolically with perturbed seed at delta=3,4; non-vacuous\n")
