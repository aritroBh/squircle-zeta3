# 05_sunit_search.py -- find x in the S-unit group of Q(sqrt2,sqrt3), S={2,3}, with 1-x also an S-unit,
# 0<x<1, and record the antisymmetric symbol x ^ (1-x) in the log basis
# (lr2=log sqrt2, lr3=log sqrt3, L=log(1+sqrt2), K=log(2+sqrt3), M=log(sqrt2+sqrt3)).
import itertools, mpmath as mp
mp.mp.dps = 40
def conj(e2, e3):
    r2 = e2*mp.sqrt(2); r3 = e3*mp.sqrt(3)
    return dict(r2=r2, r3=r3, a=1+r2, e=2+r3, m=r2+r3)
gens = ['r2', 'r3', 'a', 'e', 'm']
logs = [mp.log(mp.sqrt(2)), mp.log(mp.sqrt(3)), mp.log(1+mp.sqrt(2)), mp.log(2+mp.sqrt(3)), mp.log(mp.sqrt(2)+mp.sqrt(3))]
R = range(-3, 4)
found = {}
for ex in itertools.product(R, repeat=5):
    for sgn in (1, -1):
        vals = []
        for e2, e3 in ((1,1),(-1,1),(1,-1),(-1,-1)):
            c = conj(e2, e3); v = mp.mpf(sgn)
            for g, k in zip(gens, ex): v *= c[g]**k
            vals.append(v)
        x = vals[0]
        if not (0 < x < 1): continue
        N = 1
        for v in vals: N *= (1 - v)
        if abs(N) < 1e-30: continue
        # check N = +-2^i 3^j
        lg = mp.log(abs(N))
        ok = False
        for i in range(-40, 41):
            j = (lg - i*mp.log(2))/mp.log(3)
            if abs(j - mp.nint(j)) < mp.mpf(10)**-25: ok = True; break
        if not ok: continue
        # express log|1-x| in the log basis
        y = 1 - x
        if abs(mp.log(y)) < mp.mpf(10)**-25: continue
        rel2 = mp.pslq([mp.log(y)] + logs, maxcoeff=40, maxsteps=20000)
        if rel2 is None or abs(rel2[0]) != 1: continue
        cx = [k for k in ex]; cy = [-c*rel2[0] for c in rel2[1:]]
        # antisymmetric symbol x^y : coefficient of g_i ^ g_j (i<j)
        sym = tuple(cx[i]*cy[j] - cx[j]*cy[i] for i in range(5) for j in range(i+1, 5))
        key = round(float(x), 12)
        if key not in found:
            found[key] = (sgn, ex, cy, sym)
pairs = ["%s^%s" % (gens[i], gens[j]) for i in range(5) for j in range(i+1, 5)]
print("pairs:", pairs)
for x in sorted(found):
    sgn, ex, cy, sym = found[x]
    print("x=%.12f  sign=%d exps(r2,r3,a,e,m)=%s  1-x exps=%s  sym=%s" % (x, sgn, ex, cy, sym))
print(len(found), "solutions")
# rank of the symbol span
import sympy
M = sympy.Matrix([list(v[3]) for v in found.values()])
print("rank of symbol span:", M.rank())
