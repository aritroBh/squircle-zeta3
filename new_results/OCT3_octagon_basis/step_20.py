import sympy as sp
from itertools import permutations
r2 = sp.sqrt(2)
H = [(0,1),(1,1),(1,0),(sp.simplify(-(1+r2)),1),(sp.simplify(2+r2),1)]
labels = ["0","1","oo","-(1+r2)","2+r2"]
def diff(a,c): return sp.simplify(a[0]*c[1]-c[0]*a[1])
vals=set()
for a,b,c,d in permutations(range(5),4):
    num = diff(H[a],H[c])*diff(H[b],H[d])
    den = diff(H[a],H[d])*diff(H[b],H[c])
    if den==0: continue
    v = sp.simplify(num/den)
    if v in (0,1): continue
    vals.add(v)
print("distinct cross-ratio values:", len(vals))
def mapped_rep(v):
    w=sp.simplify(v)
    nw=sp.N(w)
    if nw>1 or nw<-1: w=sp.simplify(1/w)
    return w
mapped=set(mapped_rep(v) for v in vals)
print("distinct after inversion into [-1,1):", len(mapped))
bs=r2-1
zs17 = [ -1, -r2/2, -bs, (1-r2)/2, -bs**2, (4-3*r2)/8, bs**4, (2-r2)/4,
         bs**2, 1-r2/2, bs, sp.Rational(1,2), 2-r2, r2/2, 2*bs, (2+r2)/4, 12*r2-16 ]
S17={sp.simplify(z): str(z) for z in zs17}
notin=[]
for m in mapped:
    hit = any(sp.simplify(m-z)==0 for z in S17)
    if not hit: notin.append(m)
print("cross-ratio reps NOT in the 17:", notin)
for m in sorted(mapped, key=lambda t: sp.N(t)):
    tag = "in17" if any(sp.simplify(m-z)==0 for z in S17) else "MISSING"
    print("  %-22s = %-14s %s" % (m, sp.N(m,10), tag))
