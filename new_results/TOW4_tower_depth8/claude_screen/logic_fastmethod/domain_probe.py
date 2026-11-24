# read-only probe of TOWN1 data (no writes): is a_delta(d,m) == rho(delta)/rho(0) beyond delta <= 2d-4 ?
import pickle, sys, os
sys.path.insert(0, 'lib')
sys.path.insert(0, '.')
import pw
from pw import G, Fu, a_val
from fractions import Fraction as Fr
PKL = '<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl'
H = pickle.load(open(PKL, 'rb')); meta = H.pop('_meta')
# Q = (a,b,c,d): (a + b r2) + i (c + d r2)  -> need ratio; use own tiny Q(sqrt2,i) arithmetic
class Q:
    def __init__(s,a=0,b=0,c=0,d=0): s.a,s.b,s.c,s.d=Fr(a),Fr(b),Fr(c),Fr(d)
    def __mul__(s,o):
        a,b,c,d=s.a,s.b,s.c,s.d; e,f,g,h=o.a,o.b,o.c,o.d
        return Q(a*e+2*b*f-c*g-2*d*h, a*f+b*e-c*h-d*g, a*g+2*b*h+c*e+2*d*f, a*h+b*g+c*f+d*e)
    def __sub__(s,o): return Q(s.a-o.a,s.b-o.b,s.c-o.c,s.d-o.d)
    def iszero(s): return s.a==0 and s.b==0 and s.c==0 and s.d==0
    def conj2(s): return Q(s.a,-s.b,s.c,-s.d)
    def inv(s):
        t=s.conj2(); n=s*t; assert n.b==0 and n.d==0
        den=n.a*n.a+n.c*n.c; g=Q(n.a/den,0,-n.c/den,0); return t*g
print('meta', meta if not isinstance(meta, dict) else {k: meta[k] for k in list(meta)[:6]})
res = {}
for (p, q) in [(3,2),(4,3),(5,4),(4,2),(3,3),(6,0),(5,1),(4,1),(3,1),(3,0),(4,0),(5,0),(5,2),(4,4)]:
    if (p,q) not in H: print((p,q),'absent'); continue
    d = p+q
    v = H[(p,q)]
    row=[]
    for dl in range(0, 9):
        if dl not in v or d not in v[dl]:
            row.append((dl,'nodata')); continue
        rho0 = Q(*v[0][d]); rho = Q(*v[dl][d])
        ratio = rho * rho0.inv()
        a = a_val(dl, G(d), G(p))
        # compare ratio (in Q(sqrt2,i)) with a (in Q(i)): need sqrt2-part zero and equal
        ok = (ratio.b==0 and ratio.d==0 and ratio.a==a.a and ratio.c==a.b)
        row.append((dl, 'ok' if ok else 'DIFF'))
    print((p,q),'d=',d,'2d-4=',2*d-4, row)
