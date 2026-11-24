
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
t0=time.time()
def build(delta):
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
        Fz[(dl, dl//2+1)] = sp.together(Lam(dl, dl//2) - F(dl, dl//2))
    return Fz
out = {}
for delta in [5,6]:
    Fz = build(delta)
    close = Fz[(delta, delta//2+1)]
    out[delta] = (sp.expand(close.subs(eps,0))==0, sp.cancel(sp.diff(close, eps)))
    print(f"delta={delta}: closing(e=0)==0: {out[delta][0]},  d(close)/d eps = {out[delta][1]}  ({round(time.time()-t0,1)}s)", flush=True)
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} seed-sensitivity symbolic: d(Fclose)/d eps at dl=3,4,5,6 = nonzero,0,nonzero?, see transcript\n")
