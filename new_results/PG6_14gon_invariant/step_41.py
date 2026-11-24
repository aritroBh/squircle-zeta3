import sympy as sp, time
# Test: algebraic number alpha = 2cos(2pi/7), minpoly x^3+x^2-2x-1
al = sp.AlgebraicNumber(sp.Poly(sp.Symbol('x')**3+sp.Symbol('x')**2-2*sp.Symbol('x')-1), [1,0,0])  # = alpha
# Actually simpler: use a symbol with minpoly via QQ.algebraic_field? Let me try raw AlgebraicNumber arithmetic
t0=time.time()
x1 = al - 1          # 1-alpha ... wait
# Let's just test arithmetic
a = al**2 + al - 1   # = u
b = (1-al)**3
c = sp.expand(a*b + al**7)
print("type:", type(c), " =", c)
print("time %.3f"%(time.time()-t0))
# Test with a fresh symbol approach: polynomial reduction
t0=time.time()
als = sp.symbols('al')
MIN = als**3+als**2-2*als-1
def red(e):
    e = sp.expand(e)
    return sp.rem(sp.Poly(sp.expand(e), als), sp.Poly(MIN, als), als).as_expr()
e = red((1-als)**7 + (als**2+als-1)**3)
print("manual reduce:", e, " time %.3f"%(time.time()-t0))
