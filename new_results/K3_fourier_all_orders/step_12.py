import sympy as sp
import mpmath as mp
import pickle

x = sp.symbols('x')
N = 14

def trunc(poly, n):
    poly = sp.expand(poly)
    return sum(poly.coeff(x, i)*x**i for i in range(n+1))

# sech series (rational), via 1/cosh
cosh_ser = sum(x**(2*i)/sp.factorial(2*i) for i in range(N))
sech1 = sp.series(1/cosh_ser, x, 0, N).removeO().expand()
# sech^{1/2} via binomial in (sech-1), rational coefficients
u = sech1 - 1
sech12 = sp.Rational(0)
for j in range(N+1):
    sech12 += sp.binomial(sp.Rational(1,2), j) * u**j
    sech12 = trunc(sech12, N)
sech32 = trunc(sech1*sech12, N)
expx2 = sum((x/2)**i/sp.factorial(i) for i in range(N+1))
H32 = trunc(expx2*sech32, N)   # = 2^{3/2} H_true(x), rational coeffs
print("2^{3/2} H_true =", H32)

gd = sp.series(sp.asin(sp.tanh(x)), x, 0, N).removeO().expand()

# numeric check
mp.mp.dps = 50
fH = lambda xx: mp.e**(xx/2)*mp.sech(xx)**mp.mpf('1.5')
tay = mp.taylor(fH, 0, 10)
maxerr = max(abs(mp.mpf(str(H32.coeff(x,m))) - tay[m]) for m in range(11))
print("max |H32 - mpmath| =", mp.nstr(maxerr, 3))

MAXR, MAXM = 4, 9
e = {}
for r in range(MAXR+1):
    ser = trunc(gd**(2*r)*H32, MAXM)
    for m in range(MAXM+1):
        e[(r,m)] = sp.Rational(sp.expand(ser).coeff(x, m))
print("\ne[r][m] (= 2^{3/2}[x^m](gd)^{2r}H_true), rationals:")
for r in range(MAXR+1):
    print(r, [e[(r,m)] for m in range(MAXM+1)])

u2 = sp.symbols('u2')
def g(k2val, m):
    return sp.simplify(sum((-1)**r * (4*k2val)**r / sp.factorial(2*r) * e[(r,m)] for r in range(m//2+1)) * sp.Rational(2)**sp.Rational(-3,2))
print("\ng0 =", g(0,0), "(want sqrt2/4); g1 =", g(0,1), "(want sqrt2/8); g2 =", g(u2,2), "(want -5sqrt2/32 - sqrt2 u2/2)")
print("g3 =", g(u2,3))
with open("e_rm.pkl","wb") as f:
    pickle.dump({f"{r},{m}": str(e[(r,m)]) for r in range(MAXR+1) for m in range(MAXM+1)}, f)
print("saved e_rm.pkl")
