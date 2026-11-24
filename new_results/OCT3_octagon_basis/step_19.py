import sympy as sp
from itertools import permutations
r2 = sp.sqrt(2)
pts = [0, 1, sp.oo, -(1+r2), 2+r2]
def cr(a,b,c,d):
    # (a-c)(b-d)/((a-d)(b-c)) with infinity handling
    def sub(x,y):
        if x is sp.oo and y is sp.oo: return sp.nan
        if x is sp.oo: return sp.oo
        if y is sp.oo: return -sp.oo
        return x-y
    num = sub(a,c)*sub(b,d); den = sub(a,d)*sub(b,c)
    if num is sp.nan or den is sp.nan: return None
    if den == 0: return sp.oo
    if num is sp.oo and den is sp.oo: return sp.Rational(1)  # both infinite: ratio of leading coeffs = 1
    if num is sp.oo or den is sp.oo: return sp.oo
    if den is -sp.oo: return sp.Rational(0)
    return sp.simplify(num/den)
vals=set()
for a,b,c,d in permutations(pts,4):
    v=cr(a,b,c,d)
    if v is None or v is sp.oo: continue
    v=sp.simplify(v)
    if v in (0,1): continue
    vals.add(v)
print("distinct cross-ratio values:", len(vals))
def in_unit(z):
    return sp.simplify(z)
mapped=set()
for v in vals:
    w=v
    if sp.N(w)>1 or sp.N(w)<-1:
        w=sp.simplify(1/w)
    mapped.add(sp.simplify(w))
print("after |z|>1 -> 1/z:", len(mapped))
for m in sorted(mapped, key=lambda t: sp.N(t)):
    print("  ", m, "=", sp.N(m,12))
# our 17
bs=r2-1
zs17 = [ -1, -r2/2, -bs, (1-r2)/2, -bs**2, (4-3*r2)/8, bs**4, (2-r2)/4,
         bs**2, 1-r2/2, bs, sp.Rational(1,2), 2-r2, r2/2, 2*bs, (2+r2)/4, 12*r2-16 ]
S17=set(sp.simplify(z) for z in zs17)
print("\ncross-ratios not in the 17:", [m for m in mapped if m not in S17])
print("17 not hit by cross-ratios:", [z for z in S17 if z not in mapped])
