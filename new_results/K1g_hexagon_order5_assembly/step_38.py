import sympy as sp, mpmath as mp
mp.mp.dps=80
lA0 = sp.Symbol('lA0')
v = mp.log(mp.mpf(1)/2)
print("type v:", type(v))
f1 = sp.Float(v)
print("sp.Float(mpf):", f1, type(f1))
e = lA0**2 + 2*lA0
print("subs test:", e.subs({lA0: f1}))
# The issue might be precision. Try with explicit precision
f2 = sp.Float(v, 60)
print("sp.Float(mpf,60):", f2)
print("subs test2:", e.subs({lA0: f2}))
