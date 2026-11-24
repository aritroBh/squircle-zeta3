
import mpmath as mp, sympy as sp, av9core as C
h, e = C.selfdual_coeffs(8)
print("P_n(1,1) == 2P_n(1,0) for all n (crux of zeta(n) cancellation):")
for n in range(1,9):
    en = sp.expand(e[n])
    deg_n = sum(cf*C.lv**a*C.lw**b for (a,b),cf in sp.Poly(en, C.lv, C.lw).terms() if a+b==n)
    p11 = deg_n.subs({C.lv:1, C.lw:1}); p10 = deg_n.subs({C.lv:1, C.lw:0})
    print(f"n={n}: equal={sp.simplify(p11-2*p10)==0}   deg_n = {sp.nsimplify(deg_n)}")

mp.mp.dps = 80
def Vnum(etaval):
    q = 1+mp.mpf(etaval)
    f = lambda v: (2/q)*(v*(1-v))**(1/q-1)*(v**(q-1)+(1-v)**(q-1))**(1/q)
    return mp.quad(f, [mp.mpf('1e-12'), mp.mpf('1e-6'), mp.mpf('1e-3'), mp.mpf('0.5')])
print("\nDuality with careful quadrature:")
for ev in ['0.25','0.1','0.05']:
    a = Vnum(ev); b = Vnum(mp.mpf(1)/mp.mpf(ev))
    print(f"eta={ev}: V={mp.nstr(a,22)}  V(1/eta)={mp.nstr(b,22)}  diff={mp.nstr(abs(a-b),3)}")
