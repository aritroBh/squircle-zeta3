
import os, sys, pickle, time
import sympy as sp
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
sys.path.insert(0, work + "/tow4_copy")
import tow2lib as T2, tow4lib as T4
d, m = T2.d, T2.m
i = sp.I
tau1 = T4.tau1f; tau2 = T4.tau2f
# numeric (d,m) probe with perturbed seed at each depth: does closing change?
import random
random.seed(1)
def build_numeric(delta, seed_eps, dv, mv):
    # F family evaluated at (dv,mv) with seed perturbation at (delta,0); parents need evaluation at (dv+1,...) too
    # => build rational functions of eps only, at fixed (d,m), with shifts: we need F at (d+1,m+1),(d+1,m) too.
    # So build THREE families: at (d,m),(d+1,m+1),(d+1,m). Recursive in depth; seed enters all via its own (d,m).
    pass
# simpler: symbolic in eps, numeric in (d,m), families at the 3 needed points per depth step.
# But shifts nest: Lam at (d,m) needs F at (d+1,*); Lam at (d+1,*) needs F at (d+2,*) ... up to d+delta.
def fam(delta, eps, dv, mv):
    # returns dict F[(dp,s)] as numbers (complex rationals) at base point (dv,mv), computing shifted tables recursively
    from fractions import Fraction as Fr
    # precompute arho values at all needed points
    pts = {}
    for dl in range(0, delta+1):
        for s in range(0, dl//2+2):
            pts[(dl,s)] = {}
    # F tables at shift t = 0..delta: point (dv+t, mv+min(t,1)) -- need (d+t, m+1) and (d+t, m) both.
    # We'll compute F at every point (dv+t, mv+u) for t=0..delta, u in {0,1}.
    a = {dl: T4.arho(dl) for dl in range(0, delta+1)}
    aval = {}
    for dl in range(0, delta+1):
        for t in range(0, delta+1):
            for u in [0,1]:
                aval[(dl,t,u)] = a[dl].subs({d: dv+t, m: mv+u})
    t1 = {n: tau1(n) for n in range(0, delta+1)}
    t2 = {n: tau2(n) for n in range(0, delta+1)}
    F = {}
    for dl in range(0, delta+1):
        for t in range(0, delta-dl+1):
            for u in [0,1]:
                F[(dl,0,t,u)] = aval[(dl,t,u)] + (eps if dl==delta else sp.Integer(0))
        for s in range(0, dl//2+1):
            for t in range(0, delta-dl+1):
                for u in [0,1]:
                    tot = sp.Integer(0)
                    for dp in range(0, dl+1):
                        fA = F.get((dp,s,t+1,1), sp.Integer(0)) if s <= dp//2 else sp.Integer(0)
                        fB = F.get((dp,s,t+1,u), sp.Integer(0)) if s <= dp//2 else sp.Integer(0)
                        tot += sp.I*(2*(dv+t)-1)/(2*(dv+t)-dp-1) * (t1[dl-dp]*fA - t2[dl-dp]*fB)
                    F[(dl,s+1,t,u)] = sp.expand(tot - F[(dl,s,t,u)])
    return F
for delta in [5,6,7,8]:
    dv, mv = 13, 5
    Fp = fam(delta, sp.Symbol('e'), dv, mv)
    close = Fp[(delta, delta//2+1, 0, 0)]
    deriv = sp.diff(close, sp.Symbol('e'))
    print(f"delta={delta}: closing at (13,5) = {sp.simplify(close.subs(sp.Symbol('e'),0))}, d/d eps = {sp.simplify(deriv)}")
