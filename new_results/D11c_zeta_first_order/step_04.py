# step_04.py -- D11c(1) rhombus: numerical check of a_3(q)=(1-q^2)(zeta(3)-pi^2/3), smooth t-integral + interpolation/Richardson
import mpmath as mp
mp.mp.dps = 50

def perim_rhombus_t(p, a, b):
    eta = 1/p - 1
    def f(t):
        return mp.sqrt(a**2*(1-t)**(2*eta) + b**2*t**(2*eta))
    return 4*(1+eta)*mp.quad(f, [0, mp.mpf('0.5'), 1])

def a3_interp(a, b, r, h, deg=6):
    npts = deg+1
    xs = [h*(i - npts//2) for i in range(npts)]
    ys = [perim_rhombus_t(1+x, a, b)/r for x in xs]
    A = mp.matrix(npts, npts)
    for i, x in enumerate(xs):
        for j in range(npts): A[i,j] = x**j
    cs = mp.lu_solve(A, mp.matrix(ys))
    return cs[3]

for (a,b) in [(1,1),(2,1),(3,2),(5,1)]:
    r = mp.sqrt(a**2+b**2); q = (a**2-b**2)/(a**2+b**2)
    v1 = a3_interp(a,b,r, mp.mpf('0.008'))
    v2 = a3_interp(a,b,r, mp.mpf('0.004'))
    a3 = v2 + (v2-v1)/15
    pred = (1-q**2)*(mp.zeta(3)-mp.pi**2/3)
    print(f"a={a},b={b}: a3 fit={mp.nstr(a3,22)}")
    print(f"          pred={mp.nstr(pred,22)}   err={mp.nstr(abs(a3-pred),4)}")
