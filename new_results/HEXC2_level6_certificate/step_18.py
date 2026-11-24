import json
from fractions import Fraction as Fr
from itertools import permutations, product

# ---- exact Q(zeta), zeta = e^{i pi/3}, zeta^2 = zeta - 1 ----
# element: (a, b) = a + b*zeta, a,b Fraction
def zadd(x,y): return (x[0]+y[0], x[1]+y[1])
def zneg(x): return (-x[0], -x[1])
def zsub(x,y): return (x[0]-y[0], x[1]-y[1])
def zmul(x,y):
    a,b=x; c,d=y
    return (a*c-b*d, a*d+b*c+b*d)
def zinv(x):
    a,b=x; n = a*a+a*b+b*b
    return ((a+b)/n, -b/n)
def zdiv(x,y): return zmul(x, zinv(y))
def zpow(x,n):
    r=(Fr(1),Fr(0))
    for _ in range(n): r=zmul(r,x)
    return r
Z0=(Fr(0),Fr(0)); Z1=(Fr(1),Fr(0)); ZETA=(Fr(0),Fr(1))
MU6 = [zpow(ZETA,j) for j in range(6)]   # zeta^j, j=0..5
def zstr(x):
    for j,m in enumerate(MU6):
        if m==x: return f'z{j}'
    return str(x)

# Mobius map (A,B,C,D): z |-> (A z + B)/(C z + D), entries in Q(zeta), INF=None
INF=None
def mobius_from_triple(z1,w1,z2,w2,z3,w3):
    # solve via determinant condition; use formula with cross ratios:
    # f(z) = w such that cr(w,w1,w2,w3) = cr(z,z1,z2,z3)
    def cr(x,a,b,c):
        # (x-a)(b-c)/((x-c)(b-a))
        if x is INF: x_a, x_c = Z1, Z1
        else: x_a, x_c = zsub(x,a), zsub(x,c)
        b_c = zsub(b,c); b_a = zsub(b,a)
        return zdiv(zmul(x_a,b_c), zmul(x_c,b_a))
    # cr_z(z) = cr(z,z1,z2,z3); want cr_w(f(z)) = cr_z(z) => f = inv(cr_w) o cr_z
    # build both as Mobius maps and compose
    def cr_map(a,b,c):
        # x |-> (x-a)(b-c)/((x-c)(b-a)) : A=(b-c), B=-a(b-c), C=(b-a), D=-c(b-a)
        bc = zsub(b,c); ba = zsub(b,a)
        A = bc; B = zneg(zmul(a,bc)) if a is not INF else zneg(bc)
        C = ba; D = zneg(zmul(c,ba)) if c is not INF else zneg(ba)
        if a is INF: A, B = Z0, Z1  # handle INF: cr(INF,a,b,c) = (b-c)/(b-a)
        if a is INF:
            A = Z0; B = zdiv(zsub(b,c), zsub(b,a)) if False else bc  # careful
        return [A,B,C,D]
    # simpler: construct map sending z1->w1, z2->w2, z3->w3 directly by solving
    # f(z) = (A z + B)/(C z + D); use the standard explicit formula:
    def det2(p,q,r,s): return zsub(zmul(p,s), zmul(q,r))
    def val(pt, xy):
        # coordinate helper: represent INF
        return pt
    z1c = z1; z2c=z2; z3c=z3
    # Use: f(z) = w iff (w-w1)(w2-w3)(z-z2)(z1-z3)... use cross-ratio equality
    # cr(f(z), w1,w2,w3) = cr(z, z1,z2,z3)
    # Compose map M1: z->cr(z,z1,z2,z3) with inverse of M2: w->cr(w,w1,w2,w3)
    def cr_mobius(a,b,c):
        # map x |-> cr(x,a,b,c) = (x-a)(b-c)/((x-c)(b-a))
        if a is INF:
            # cr = (b-c)/(x-c)
            return [Z0, zsub(b,c), Z1, zneg(c)]
        if c is INF:
            # cr = (x-a)/(b-a)
            return [Z1, zneg(a), Z0, zsub(b,a)]
        bc = zsub(b,c); ba = zsub(b,a)
        return [bc, zneg(zmul(a,bc)), ba, zneg(zmul(c,ba))]
    def minv(M):
        A,B,C,D = M
        return [D, zneg(B), zneg(C), A]
    def mcomp(M,N):
        # M o N
        A,B,C,D = M; a,b,c,d = N
        return [zadd(zmul(A,a),zmul(B,c)), zadd(zmul(A,b),zmul(B,d)),
                zadd(zmul(C,a),zmul(D,c)), zadd(zmul(C,b),zmul(D,d))]
    M = mcomp(minv(cr_mobius(w1,w2,w3)), cr_mobius(z1,z2,z3))
    return M

def mapply(M, z):
    A,B,C,D = M
    if z is INF:
        if C==Z0: return INF
        return zdiv(A,C)
    num = zadd(zmul(A,z),B); den = zadd(zmul(C,z),D)
    if den==Z0: return INF
    return zdiv(num,den)

PTS = [Z0, (Fr(1,2),Fr(0)), Z1, (Fr(-1),Fr(0)), INF]
names = ['0','1/2','1','-1','inf']
maps = []
for w0,w1,wI in permutations(range(6),3):
    M = mobius_from_triple(Z0, MU6[w0], Z1, MU6[w1], INF, MU6[wI])
    imgs = [mapply(M,p) for p in PTS]
    if all(any(im==MU6[j] for j in range(6)) if im is not INF else False for im in imgs):
        exps = tuple(next(j for j in range(6) if MU6[j]==im) for im in imgs)
        maps.append((M, exps))
print('number of maps sending {0,1/2,1,-1,inf} into mu6:', len(maps))
for M, exps in maps:
    print('  exps(0,1/2,1,-1,inf) =', exps, '  map coeffs:', [zstr(c) for c in M])
json.dump([{'exps':list(e)} for M,e in maps], open('mobius_maps.json','w'), indent=1)
with open('progress.log','a') as f: f.write(f'step_02: {len(maps)} Mobius maps {{0,1/2,1,-1,inf}} -> mu6, exps saved\n')
