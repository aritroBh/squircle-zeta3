import sympy as sp, pickle
x = sp.symbols('x')
N = 20
# Gudermannian: gd(x) = asin(tanh x)
gd = sp.series(sp.asin(sp.tanh(x)), x, 0, N+2).removeO().expand()
# H(x) = 2^{-3/2} e^{x/2} sech^{3/2} x ; work with 2^{3/2} H = e^{x/2} sech^{3/2}
cosh_ser = sum(x**(2*i)/sp.factorial(2*i) for i in range(N))
sech1 = sp.series(1/cosh_ser, x, 0, N+1).removeO().expand()
u = sech1 - 1
sech12 = sp.Rational(0)
def trunc(p, n=N):
    p = sp.expand(p)
    return sum(p.coeff(x, i)*x**i for i in range(n+1))
for j in range(N+1):
    sech12 = trunc(sech12 + sp.binomial(sp.Rational(1,2), j)*u**j)
sech32 = trunc(sech1*sech12)
expx2 = sum((x/2)**i/sp.factorial(i) for i in range(N+1))
H32 = trunc(expx2*sech32)   # = 2^{3/2} H
MAXR, MAXM = 4, 12
e = {}
for r in range(MAXR+1):
    ser = trunc(gd**(2*r)*H32, MAXM)
    for m in range(MAXM+1):
        e[(r,m)] = sp.Rational(sp.expand(ser).coeff(x, m))
# compare with K3's e_rm.pkl
theirs = {}
with open('../K3_fourier_all_orders/e_rm.pkl','rb') as f:
    for key, val in pickle.load(f).items():
        r, m = map(int, key.split(',')); theirs[(r,m)] = sp.Rational(val)
mism = [(k, e[k], theirs[k]) for k in theirs if e.get(k) != theirs[k]]
print("e_rm mismatches vs K3:", mism if mism else "NONE -- exact match")
print("spot: e[0,2] =", e[(0,2)], " e[1,4] =", e[(1,4)], " e[2,6] =", e[(2,6)], " e[3,8] =", e[(3,8)])
print("claimed e_{r,2r}=1:", all(e[(r,2*r)]==1 for r in range(MAXR+1)))
print("claimed e_{r,2r+1}=1/2:", all(e[(r,2*r+1)]==sp.Rational(1,2) for r in range(MAXR+1) if 2*r+1<=MAXM))
print("claimed e_{r,2r+2}=-5/8-r/3:", all(e[(r,2*r+2)]==-sp.Rational(5,8)-sp.Rational(r,3) for r in range(MAXR+1) if 2*r+2<=MAXM))
pickle.dump({f"{r},{m}": str(e[(r,m)]) for r in range(MAXR+1) for m in range(MAXM+1)}, open('e_rm_mine.pkl','wb'))
