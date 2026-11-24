import sys, os, time
work = '<project>/research/08_theory/TOW2R_review/tow2_copy'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow2lib as T
d, m = T.d, T.m
i = T.i

# V(3) non-adjacent: rho(3)/e00 + e(3,1)/e00 == sum_{dp=0..3} [feeds of parents' rho(dp)]
def feed(delta, s_top_of_parent):
    """RHS of the (delta, s=0)-identity: feeds of parents' e(delta', 0) = rho(delta')."""
    tot = 0
    for dp in range(0, delta+1):
        cA = T.tau1(delta-dp) * (m+1) / (2*d-dp-1)
        cB = T.tau2(delta-dp) * (d-m+1) / (2*d-dp-1)
        tot += cA * T.arho(dp).subs({d: d+1, m: m+1}) * T.rA0()
        tot += cB * T.arho(dp).subs({d: d+1, m: m}) * T.rB0()
    return tot

t0 = time.time()
rho3 = T.arho(3)
e31 = i*(d-3)/(2*d-3)
lhs3 = rho3 + e31
rhs3 = feed(3, 0)
rem3 = sp.cancel(sp.together(lhs3 - rhs3))
print('V(3) non-adjacent remainder:', rem3, ' (', round(time.time()-t0,1),'s )')
