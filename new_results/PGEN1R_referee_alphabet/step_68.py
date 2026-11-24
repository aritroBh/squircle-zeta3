from fractions import Fraction as Fr
def wedge(u, v):
    out = {}
    for a, ca in u.items():
        for b, cb in v.items():
            if a == b: continue
            key, sg = ((a,b),1) if a < b else ((b,a),-1)
            out[key] = out.get(key, Fr(0)) + sg*ca*cb
    return {k:c for k,c in out.items() if c}
def delta_li3(vz, vomz):
    w = wedge(vz, vomz)
    out = {}
    for ab, cab in w.items():
        for l3, c3 in vz.items():
            out[(ab, l3)] = out.get((ab, l3), Fr(0)) + cab*c3
    return {k:c for k,c in out.items() if c}
def dadd(*ds):
    out = {}
    for d in ds:
        for k,c in d.items():
            out[k] = out.get(k, Fr(0))+c
            if out[k]==0: del out[k]
    return out
def dsc(d, c):
    return {k: v*c for k,v in d.items() if v*c}
def neg(v): return {k:-c for k,c in v.items() if c}
def sub(u,v):
    r = dict(u)
    for k,c in v.items():
        r[k] = r.get(k,Fr(0))-c
        if r[k]==0: del r[k]
    return r

# --- formal Landen triple: z, 1-z, z/(z-1) with letters [z]=Z, [1-z]=W ---
Z, W = {'Z':1}, {'W':1}
dz   = delta_li3(Z, W)                       # z, 1-z
d1z  = delta_li3(W, Z)                       # 1-z, z
# z/(z-1): letter Z - W ; 1 - z/(z-1) = (z-1-z)/(z-1) = -1/(z-1): letter -(W) -> W (sign dropped), i.e. [1-z/(z-1)] = -[W]: vector -W
dzz1 = delta_li3(sub(Z,W), neg(W))
print('Landen triple delta (formal):', dadd(dz, d1z, dzz1))
# --- formal duplication: z^2, z, -z : delta{z^2} - 4 delta{z} - 4 delta{-z} ---
# letters: [z]=Z, [1-z]=W, [1+z]=P ; z^2: 2Z ; 1-z^2 = (1-z)(1+z): W+P ; -z: Z (sign dropped), 1+z: P
Z2 = {'Z':2}
d_z2 = delta_li3(Z2, dadd(W, {'P':1}))
d_z  = delta_li3(Z, W)
d_mz = delta_li3(Z, {'P':1})
dup = dadd(d_z2, dsc(d_z, Fr(-4)), dsc(d_mz, Fr(-4)))
print('duplication delta (formal):', dup)
# --- the specific octagon relations ---
vz8  = {'d': {'r2':1,'u':-1}, 'k': {'u':-1}, 'mr2h': {'r2':-1}, 'r2h': {'r2':-1}, 'x1': {'u':-1}, 'x2': {'r2':-1,'u':-1}}
vo8  = {'d': {'u':-1}, 'k': {'r2':1,'u':-1}, 'mr2h': {'u':1,'r2':-1}, 'r2h': {'r2':-1,'u':-1}, 'x1': {'r2':1}, 'x2': {'r2':-1}}
T1 = dadd(delta_li3(vz8['k'],vo8['k']), delta_li3(vz8['d'],vo8['d']), delta_li3(vz8['mr2h'],vo8['mr2h']))
T2 = dadd(delta_li3(vz8['x2'],vo8['x2']), delta_li3(vz8['r2h'],vo8['r2h']), delta_li3(vz8['x1'],vo8['x1']))
DUP8 = dadd(delta_li3(vz8['r2h'],vo8['r2h']), delta_li3(vz8['mr2h'],vo8['mr2h']))
print('T1{k,d,mr2h} delta:', T1, ' T2{x2,r2h,x1} delta:', T2, ' dupl{r2h,mr2h} delta:', DUP8)
# wedge^3 of the octagon delta: trivially 0 (2 letters) -- will test on 12-gon
open('progress.log','a').write('step: formal Landen/duplication delta=0; octagon T1,T2,dup delta=0 exactly\n')
