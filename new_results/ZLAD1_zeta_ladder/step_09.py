import time, sympy as sp
u, y, a, b, n = sp.symbols('u y a b n')
KMAX = 6
t0 = time.time()
A = -(n - 1)*a + n*b
l = b - a
logE = u*A + (1 - 2*n)*sp.log(1 + y*(sp.exp(u*l) - 1))
logE_ser = sp.series(logE, u, 0, KMAX + 1).removeO()
E_ser = sp.series(sp.exp(logE_ser), u, 0, KMAX + 1).removeO()
print('series time', time.time()-t0, flush=True)
t0 = time.time()
Ek = {}
for k in range(KMAX + 1):
    Ek[k] = sp.Poly(sp.expand(E_ser.coeff(u, k)), y, a, b)
    print(k, 'terms:', len(Ek[k].terms()), 't=', round(time.time()-t0,1), flush=True)
import pickle
pickle.dump({k: Ek[k].as_expr() for k in Ek}, open('Ek_polys_K6.pkl','wb'))
print('saved')
