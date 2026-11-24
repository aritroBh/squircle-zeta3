
import os, sys, pickle, time
import sympy as sp
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
sys.path.insert(0, work + "/tow4_copy")
import tow2lib as T2, tow4lib as T4
d, m = T2.d, T2.m
i = sp.I
tau1 = T4.tau1f; tau2 = T4.tau2f
def fam(delta, dv, mv, eps=None):
    a = {dl: T4.arho(dl) for dl in range(0, delta+1)}
    aval = {}
    for dl in range(0, delta+1):
        for t in range(0, delta+2):
            for u in [0,1]:
                aval[(dl,t,u)] = a[dl].subs({d: dv+t, m: mv+u})
    t1 = {n: tau1(n) for n in range(0, delta+1)}
    t2 = {n: tau2(n) for n in range(0, delta+1)}
    F = {}
    for dl in range(0, delta+1):
        for t in range(0, delta+2):
            for u in [0,1]:
                F[(dl,0,t,u)] = aval[(dl,t,u)] + (eps if (eps is not None and dl==delta) else sp.Integer(0))
        for s in range(0, dl//2+1):
            for t in range(0, delta+1-s):
                for u in [0,1]:
                    tot = sp.Integer(0)
                    for dp in range(0, dl+1):
                        fA = F[(dp,s,t+1,1)] if s <= dp//2 else sp.Integer(0)
                        fB = F[(dp,s,t+1,u)] if s <= dp//2 else sp.Integer(0)
                        tot += sp.I*(2*(dv+t)-1)/(2*(dv+t)-dp-1) * (t1[dl-dp]*fA - t2[dl-dp]*fB)
                    F[(dl,s+1,t,u)] = sp.expand(tot - F[(dl,s,t,u)])
    return F
dv, mv = 13, 5
F = fam(6, dv, mv)
ok = all(sp.expand(F[(dl,s,0,0)] - T4.forced(dl,s,dv,mv))==0 for dl in range(0,7) for s in range(0,dl//2+1))
print("F matches T4.forced at (13,5), dl<=6:", ok)
print("closings dl<=6 at (13,5):", [sp.simplify(F[(dl,dl//2+1,0,0)]) for dl in range(0,7)])
e = sp.Symbol('e')
for delta in [5,6,7,8]:
    Fp = fam(delta, dv, mv, eps=e)
    close = Fp[(delta, delta//2+1, 0, 0)]
    print(f"delta={delta}: closing(e=0) = {sp.simplify(close.subs(e,0))},  d/d eps = {sp.simplify(sp.diff(close, e))}")
