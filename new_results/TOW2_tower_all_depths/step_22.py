import sys, os
sys.dont_write_bytecode = True
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from towlib import *
from genintegrand import S_poly
from fractions import Fraction as Fr

tstar = QSI(Fr(1,2),0,Fr(1,2),0)

def S_val(d, m, l):
    """(d/dt)^l S_m(t) at t*, exact QSI (should be in Q(i))."""
    P = S_poly(d, m)
    out = ZERO
    for n, c in enumerate(P):
        if n < l: continue
        # falling factorial n!/(n-l)! * c * t*^(n-l)
        ff = 1
        for a in range(l): ff *= (n-a)
        out = out + QSI(Fr(c*ff)) * (tstar**(n-l))
    return out

# compute ratios T_l/T_0 for l=0..4, d=2..10, all m; print a few to see structure
for l in range(0,5):
    print(f'=== S^({l})(t*)/S(t*) sample values (should be Q(i))')
    for (d,m) in [(5,5),(5,3),(6,6),(6,4),(7,2),(8,8),(9,0)]:
        T0 = S_val(d,m,0); Tl = S_val(d,m,l)
        r = Tl/T0
        assert r.p1==0 and r.p3==0, (d,m,l,r)
        print(f'  (d,m)=({d},{m}): {r.p0} + {r.p2}*i')
