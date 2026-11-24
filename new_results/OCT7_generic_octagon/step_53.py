import mpmath as mp, sympy as sp, time, random
mp.mp.dps = 450
n = 60
C = sp.Integer(10)**430
vals = [mp.mpf(1), mp.log(2), mp.log(3), mp.log(5)]
random.seed(1)
for i in range(n-4):
    vals.append(mp.log(2)**random.randint(0,3) * mp.log(3)**random.randint(0,2) * mp.pi**(2*random.randint(0,1)))
def to_int(v):
    return sp.Integer(mp.nstr(v*C, 460))
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
v0 = list(L[0])
print("first reduced vector: first-n norm ~", max(abs(int(x)) for x in v0[:-1]), "last coord zero?", v0[-1]==0)
