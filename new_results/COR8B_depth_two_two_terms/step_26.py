
engine = r'''
# cor8b_sym.py -- independent exact symbol engine for COR8B
# Letters: integer 4-vectors over basis (m=[2], c=[c], p=[1-c], q=[1+c]).
# Symbols: dict { (v1,v2,v3,...) : Fraction coeff }.
# delta_2,2(a b c d) = (a^b)^(c^d) in Lambda^2(Lambda^2 V).
from fractions import Fraction as Fr

M = (1,0,0,0); C = (0,1,0,0); P = (0,0,1,0); Q = (0,0,0,1)
ZERO = (0,0,0,0)

def vadd(a,b): return tuple(x+y for x,y in zip(a,b))
def vneg(a): return tuple(-x for x in a)
def vsub(a,b): return tuple(x-y for x,y in zip(a,b))
def vmul(k,a): return tuple(k*x for x in a)

def sym_add(s1, s2):
    out = dict(s1)
    for k,v in s2.items():
        out[k] = out.get(k, Fr(0)) + v
        if out[k] == 0: del out[k]
    return out
def sym_scale(s, a):
    a = Fr(a)
    return {k: v*a for k,v in s.items() if v*a != 0}

# ---- two-variable symbols over formal basis {X,Y,U,V,W} (Z=xy -> X+Y) ----
# basis indices: 0=X[x],1=Y[y],2=U[1-x],3=V[1-y],4=W[1-xy]
XV=(1,0,0,0,0); YV=(0,1,0,0,0); UV=(0,0,1,0,0); VV=(0,0,0,1,0); WV=(0,0,0,0,1)
def v5add(a,b): return tuple(x+y for x,y in zip(a,b))
def v5sub(a,b): return tuple(x-y for x,y in zip(a,b))
ZV = v5add(XV,YV)   # [xy]=[x]+[y]
YVmV = v5sub(YV,VV) # [y/(y-1)] = [y]-[1-y]

def T(*vs): return {tuple(vs): Fr(1)}      # elementary tensor
def tprod(s, v):                            # append letter v at the end
    return {k+(v,): c for k,c in s.items()}

# weight-2/3 symbols (convB, from verified dlog recursions)
S_Li1x  = { (vneg(UV),) : Fr(-1) }  # placeholder not used
def S_Li2(which):
    # S(Li2(z)) = -[1-z]^[z]; which in {'x','xy'}
    if which=='x': return { (UV, XV): Fr(-1) }
    if which=='xy': return { (WV, ZV): Fr(-1) }

def S_Li11():
    # dLi11 = -Li1(xy)[dlog y + dlog(1-x) - dlog(1-y)] - Li1(x) dlog(1-y)
    # S(-Li1(z)) = [1-z]
    s = {}
    for L,c in [(YV,1),(UV,1),(v5sub(ZERO5,VV),1)]:
        pass
    s = { (WV, YV): Fr(1), (WV, UV): Fr(1), (WV, VV): Fr(-1), (UV, VV): Fr(1) }
    return s
ZERO5 = (0,0,0,0,0)

def S_Li12():
    # dLi12 = -Li2(xy) dlog(1-x) + Li11 dlog y
    s = { (WV, ZV, UV): Fr(1) }       # -S(Li2(xy)) = +[1-xy][xy]; tensor [1-x]
    sL = S_Li11()
    for k,c in sL.items():
        s[k+(YV,)] = s.get(k+(YV,), Fr(0)) + c
    return {k:v for k,v in s.items() if v}

def S_Li21():
    # dLi21 = Li11 dlog x - Li2(xy) dlog(y/(y-1)) - Li2(x) dlog(1-y)
    s = {}
    for k,c in S_Li11().items():
        s[k+(XV,)] = s.get(k+(XV,), Fr(0)) + c
    s[(WV, ZV, YV)] = s.get((WV,ZV,YV), Fr(0)) + 1      # -S(Li2(xy)) ^ [y]
    s[(WV, ZV, VV)] = s.get((WV,ZV,VV), Fr(0)) - 1      # -S(Li2(xy)) ^ (-[1-y])
    s[(UV, XV, VV)] = s.get((UV,XV,VV), Fr(0)) + 1      # -S(Li2(x)) ^ [1-y]
    return {k:v for k,v in s.items() if v}

def S_Li22():
    # dLi22 = Li12 dlog x + Li21 dlog y
    s = {}
    for k,c in S_Li12().items():
        s[k+(XV,)] = s.get(k+(XV,), Fr(0)) + c
    for k,c in S_Li21().items():
        s[k+(YV,)] = s.get(k+(YV,), Fr(0)) + c
    return {k:v for k,v in s.items() if v}

def S_Li31():
    # dLi31 = Li21 dlog x - Li3(xy) dlog(y/(y-1)) - Li3(x) dlog(1-y)
    # S(Li3(z)) = -[1-z][z][z]
    s = {}
    for k,c in S_Li21().items():
        s[k+(XV,)] = s.get(k+(XV,), Fr(0)) + c
    s[(WV, ZV, ZV, YV)] = s.get((WV,ZV,ZV,YV), Fr(0)) + 1
    s[(WV, ZV, ZV, VV)] = s.get((WV,ZV,ZV,VV), Fr(0)) - 1
    s[(UV, XV, XV, VV)] = s.get((UV,XV,XV,VV), Fr(0)) + 1
    return {k:v for k,v in s.items() if v}

def S_Li13():
    # dLi13 = -Li3(xy) dlog(1-x) + Li12 dlog y
    s = { (WV, ZV, ZV, UV): Fr(1) }
    for k,c in S_Li12().items():
        s[k+(YV,)] = s.get(k+(YV,), Fr(0)) + c
    return {k:v for k,v in s.items() if v}

def S_Li4z(which):
    # S(Li4(z)) = -[1-z][z]^3
    if which=='x': return { (UV, XV, XV, XV): Fr(-1) }
    if which=='xy': return { (WV, ZV, ZV, ZV): Fr(-1) }
    if which=='y': return { (VV, YV, YV, YV): Fr(-1) }

# ---- specialization map: 5-letter -> (m,c,p,q) for y = 1-c or y = 1+c ----
def spec_map(sign, ximg=None):
    # x=1/2. sign=-1: y=1-c; sign=+1: y=1+c
    # returns dict mapping 5-letter basis vectors to 4-vectors
    if ximg is None:
        ximg = vneg(M)      # [1/2] = -[2]
        uimg = vneg(M)      # [1/2]
    if sign == -1:
        yimg = P; vimg = C; wimg = vsub(Q, M)   # [1-xy] = [(1+c)/2]
    else:
        yimg = Q; vimg = C; wimg = vsub(P, M)   # [1-xy] = [(1-c)/2]
    return {XV:ximg, YV:yimg, UV:uimg, VV:vimg, WV:wimg}

def specialize(sym5, smap):
    out = {}
    for k, c in sym5.items():
        nk = tuple(smap[v] for v in k)
        out[nk] = out.get(nk, Fr(0)) + c
    return {k:v for k,v in out.items() if v}

# ---- generic pair (x0,y0): letter images as 4-vectors (given directly) ----
def specialize_pair(sym5, imgs):
    # imgs: dict with keys 'X','U','Y','V','W' -> 4-vectors
    m = {XV:imgs['X'], UV:imgs['U'], YV:imgs['Y'], VV:imgs['V'], WV:imgs['W']}
    return specialize(sym5, m)

# ---- wedge algebra over (m,c,p,q): V dim 4 ----
WEDGE_BASIS = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]  # index pairs
def wedge2(a, b):
    # a,b 4-vectors -> dict over wedge basis index
    out = {}
    for idx,(i,j) in enumerate(WEDGE_BASIS):
        c = a[i]*b[j] - a[j]*b[i]
        if c: out[idx] = out.get(idx,0)+c
    return {k:v for k,v in out.items() if v}

def wedge2_add(w1,w2,sgn=1):
    out = dict(w1)
    for k,v in w2.items():
        out[k] = out.get(k,0)+sgn*v
        if out[k]==0: del out[k]
    return out

def wedge4(A, B):
    # A,B: dicts over wedge-basis indices -> antisymmetric product over 6-dim space
    # represented as dict {(i,j): coeff} with i<j over wedge indices
    out = {}
    for i,ci in A.items():
        for j,cj in B.items():
            if i==j: continue
            if i<j: out[(i,j)] = out.get((i,j),0) + ci*cj
            else:   out[(j,i)] = out.get((j,i),0) - ci*cj
    return {k:v for k,v in out.items() if v}

def delta22(sym):
    # sym: dict { (v1,v2,v3,v4): coeff } -> dict {(i,j): Fraction}
    out = {}
    for (a,b,c,d), cf in sym.items():
        w = wedge4(wedge2(a,b), wedge2(c,d))
        for k,v in w.items():
            out[k] = out.get(k, Fr(0)) + cf*v
            if out[k]==0: del out[k]
    return out

def d22_scale(d, a):
    a = Fr(a)
    return {k: v*a for k,v in d.items() if v*a}
def d22_add(d1,d2):
    out = dict(d1)
    for k,v in d2.items():
        out[k] = out.get(k, Fr(0)) + v
        if out[k]==0: del out[k]
    return out

WNAME = {0:'mc',1:'mp',2:'mq',3:'cp',4:'cq',5:'pq'}
def show_d22(d):
    parts = []
    for (i,j),v in sorted(d.items()):
        parts.append(f"({v})*({WNAME[i]}^{WNAME[j]})")
    return " + ".join(parts) if parts else "0"
'''
open('cor8b_sym.py','w').write(engine)
print("engine written")
