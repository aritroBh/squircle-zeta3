
import os, sys, pickle, time
import sympy as sp
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
sys.path.insert(0, work + "/tow4_copy")
import tow2lib as T2, tow4lib as T4
d, m = T2.d, T2.m
i = sp.I
tau1 = T4.tau1f; tau2 = T4.tau2f
def fam(delta, dv, mv):
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
                F[(dl,0,t,u)] = aval[(dl,t,u)]
        for s in range(0, dl//2+1):
            for t in range(0, delta-dl+1):
                for u in [0,1]:
                    tot = sp.Integer(0)
                    for dp in range(0, dl+1):
                        fA = F[(dp,s,t+1,1)] if s <= dp//2 else sp.Integer(0)
                        fB = F[(dp,s,t+1,u)] if s <= dp//2 else sp.Integer(0)
                        tot += sp.I*(2*(dv+t)-1)/(2*(dv+t)-dp-1) * (t1[dl-dp]*fA - t2[dl-dp]*fB)
                    F[(dl,s+1,t,u)] = sp.expand(tot - F[(dl,s,t,u)])
    return F
# test at delta=4: compare F values with T4.forced and closing with 0
dv, mv = 13, 5
F = fam(4, dv, mv)
for dl in range(0,5):
    for s in range(0, dl//2+1):
        ref = T4.forced(dl, s, dv, mv)
        got = F[(dl,s,0,0)]
        print((dl,s), "match:", sp.expand(got-ref)==0, end="  ")
    print()
print("closings:", [(dl, sp.simplify(F[(dl, dl//2+1, 0, 0)])) for dl in range(0,5)])
