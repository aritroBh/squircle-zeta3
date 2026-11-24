import mpmath as mp, sympy as sp, time, random
mp.mp.dps = 450
n = 100
C = sp.Integer(10)**430
vals = [mp.mpf(1), mp.log(2), mp.log(3), mp.log(5)]
random.seed(2)
for i in range(n-4):
    vals.append(mp.log(2)**random.randint(0,3) * mp.log(3)**random.randint(0,2) * mp.pi**(2*random.randint(0,1)))
def to_int(v):
    return sp.Integer(int(mp.floor(v*C + mp.mpf('0.5'))))
rows = []
for i in range(n):
    r = [sp.Integer(0)]*n + [to_int(vals[i])]
    r[i] = sp.Integer(1)
    rows.append(r)
M = sp.Matrix(rows)
t0=time.time()
L = M.lll(delta=sp.Rational(3,4))
t1=time.time()
print(f"sympy LLL dim {n}x{n+1}, 430-digit entries: {t1-t0:.1f}s")
# inspect: rows with tiny last coordinate = candidate relations
Lrows = L.tolist()
cands = [r for r in Lrows if abs(int(r[-1])) < 10**390]
print("rows with |last|<1e390:", len(cands))
for r in cands[:3]:
    print("  max first-n coord:", max(abs(int(x)) for x in r[:-1]), "|last|:", abs(int(r[-1])))
# verify one: residual
r = cands[0]
resid = sum(sp.Integer(int(x))*sp.nsimplify(0) for x in r[:-1]) if False else None
import mpmath as mp2
res = mp2.fsum([int(r[i])*vals[i] for i in range(n)]) if any(r[:-1]) else 0
print("candidate residual:", mp2.nstr(abs(res),3))
print("norm scale of reduced basis (min over rows of max|first-n|):", min(max(abs(int(x)) for x in r[:-1]) for r in Lrows))
