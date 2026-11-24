import sys, time
import sympy as sp
sys.path.insert(0, '<project>/research/08_theory/TOW4_tower_depth8')
import tow4lib as T4
d, m = T4.d, T4.m
i = sp.I
tau1 = T4.tau1f; tau2 = T4.tau2f

def U(dp):
    return sum((-1)**s * T4.forced(dp, s) for s in range(0, dp//2+1))

def feedU(delta):
    tot = sp.Integer(0)
    for dp in range(0, delta+1):
        uA = U(dp).subs({d: d+1, m: m+1}, simultaneous=True)
        uB = U(dp).subs({d: d+1, m: m}, simultaneous=True)
        tot += i*(2*d-1)/(2*d-dp-1) * (tau1(delta-dp)*uA - tau2(delta-dp)*uB)
    return tot

for delta in [7, 8]:
    t0=time.time()
    expr = sp.cancel(sp.together(feedU(delta) - T4.arho(delta)))
    print(f"master x=-1 coefficient delta={delta}: feed(U(-1)) - a_{delta} == 0 : {expr==0}  ({round(time.time()-t0,1)}s)", flush=True)
