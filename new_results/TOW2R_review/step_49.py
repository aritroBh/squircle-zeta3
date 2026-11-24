import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
import sympy as sp
import genintegrand as GI
x = sp.Symbol('x')
P = {0: sp.Integer(1), 1: x*(x-1)/2, 2: x*(x-1)*(x-2)*(3*x-5)/24,
     3: x*(x-1)*(x-2)**2*(x-3)**2/48,
     4: x*(x-1)*(x-2)*(x-3)*(x-4)*(15*x**3-150*x**2+485*x-502)/5760}
bad = []
for r in range(1,5):
    for n in range(0, 30):
        sv = GI.stirling2(n, n-r)   # returns 0 for negative k
        pv = P[r].subs(x, n)
        if sp.Integer(sv) != pv:
            bad.append((r, n, sv, pv))
print('P_r(n) vs genintegrand.stirling2(n, n-r), n=0..29, r<=4:', 'ALL EQUAL' if not bad else bad[:10])
