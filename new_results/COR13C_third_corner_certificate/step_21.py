import json
from fractions import Fraction as Fr
from itertools import permutations

# ---- exact Q(zeta), zeta = e^{i pi/3}, zeta^2 = zeta - 1 (HEXC2 step_18 style) ----
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
def zscal(x,s): return (x[0]*s, x[1]*s)
Z0=(Fr(0),Fr(0)); Z1=(Fr(1),Fr(0)); ZETA=(Fr(0),Fr(1))
MU6 = [zpow(ZETA,j) for j in range(6)]
INF=None

# points of the c0=1/3 corner configuration
TH = (Fr(1,3),Fr(0)); M1 = (Fr(-1),Fr(0))
PTS = [Z0, TH, Z1, M1, INF]
names = ['0','1/3','1','-1','inf']

def mobius_from_triple(z1,w1,z2,w2,z3,w3):
    def cr_mobius(a,b,c):
        # map x |-> cr(x,a,b,c) = (x-a)(b-c)/((x-c)(b-a))
        if a is INF:
            return [Z0, zsub(b,c), Z1, zneg(c)]
        if c is INF:
            return [Z1, zneg(a), Z0, zsub(b,a)]
        bc = zsub(b,c); ba = zsub(b,a)
        return [bc, zneg(zmul(a,bc)), ba, zneg(zmul(c,ba))]
    def minv(M):
        A,B,C,D = M
        return [D, zneg(B), zneg(C), A]
    def mcomp(M,N):
        A,B,C,D = M; a,b,c,d = N
        return [zadd(zmul(A,a),zmul(B,c)), zadd(zmul(A,b),zmul(B,d)),
                zadd(zmul(C,a),zmul(D,c)), zadd(zmul(C,b),zmul(D,d))]
    return mcomp(minv(cr_mobius(w1,w2,w3)), cr_mobius(z1,z2,z3))

def mapply(M, z):
    A,B,C,D = M
    if z is INF:
        if C==Z0: return INF
        return zdiv(A,C)
    num = zadd(zmul(A,z),B); den = zadd(zmul(C,z),D)
    if den==Z0: return INF
    return zdiv(num,den)

maps = []
for w0,w1,wI in permutations(range(6),3):
    M = mobius_from_triple(Z0, MU6[w0], Z1, MU6[w1], INF, MU6[wI])
    imgs = [mapply(M,p) for p in PTS]
    if all(im is not INF and any(im==MU6[j] for j in range(6)) for im in imgs):
        exps = tuple(next(j for j in range(6) if MU6[j]==im) for im in imgs)
        if len(set(exps))==5:
            maps.append((M, exps))
print('number of maps sending {0,1/3,1,-1,inf} into mu6 (injective):', len(maps))
omitted = {}
for M, exps in maps:
    om = next(j for j in range(6) if j not in exps)
    omitted[om] = omitted.get(om,0)+1
    print('  exps(0,1/3,1,-1,inf) =', exps, ' omits', om)
print('maps per omitted root:', sorted(omitted.items()))
# relaxed: images in mu6 u {0, inf}
cnt_relaxed = 0
for w0,w1,wI in permutations(range(6),3):
    M = mobius_from_triple(Z0, MU6[w0], Z1, MU6[w1], INF, MU6[wI])
    imgs = [mapply(M,p) for p in PTS]
    ok = all(im is INF or im==Z0 or any(im==MU6[j] for j in range(6)) for im in imgs)
    if ok: cnt_relaxed += 1
print('relaxed (mu6 u {0,inf}) count:', cnt_relaxed)
json.dump([{'exps':list(e)} for M,e in maps], open('mobius_maps_13.json','w'), indent=1)
with open('progress.log','a') as f:
    f.write(f'step_21: {len(maps)} Mobius maps {{0,1/3,1,-1,inf}} -> mu6 exact over Q(zeta6); relaxed={cnt_relaxed}\n')
