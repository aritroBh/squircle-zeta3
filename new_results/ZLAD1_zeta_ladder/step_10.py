import time, sympy as sp
u, y, a, b = sp.symbols('u y a b')
KMAX = 6

def E_series_for_n(nv):
    A = -(nv - 1)*a + nv*b
    l = b - a
    logE = u*A + (1 - 2*nv)*sp.log(1 + y*(sp.exp(u*l) - 1))
    logE_ser = sp.series(logE, u, 0, KMAX + 1).removeO()
    E_ser = sp.expand(sp.series(sp.exp(logE_ser), u, 0, KMAX + 1).removeO())
    return [sp.Poly(sp.expand(E_ser.coeff(u, k)), y, a, b) for k in range(KMAX+1)]

t0=time.time()
for nv in [1,2,3]:
    Ek = E_series_for_n(nv)
    print(nv, [len(Ek[k].terms()) for k in range(KMAX+1)], round(time.time()-t0,2))
t0=time.time()
for nv in range(1, 151):
    Ek = E_series_for_n(nv)
print('150 n values:', round(time.time()-t0,1), 's')
