
# step_13b.py : (ii) top-weight check; (iii) duality numeric
import mpmath as mp, sympy as sp, av9core as C
h, e = C.selfdual_coeffs(8)
print("-- top-weight check (self-dual) --")
for n in range(1,9):
    en = sp.expand(e[n])
    deg_n = sum(cf*C.lv**a*C.lw**b for (a,b),cf in sp.Poly(en, C.lv, C.lw).terms() if a+b==n)
    pred = (C.lv**n + C.lw**n)/(2*sp.factorial(n))
    print(f"n={n}: deg-n part of e_n == (lv^n+lw^n)/(2 n!) : {sp.simplify(deg_n-pred)==0}")
mp.mp.dps = 60
def Vnum(etaval):
    q = 1+etaval
    f = lambda v: (2/q)*(v*(1-v))**(1/q-1)*(v**(q-1)+(1-v)**(q-1))**(1/q)
    return mp.quad(f, [0, mp.mpf('0.5'), 1])
print("-- V(eta) vs V(1/eta) --")
for ev in [mp.mpf('0.25'), mp.mpf('0.1'), mp.mpf('0.05'), mp.mpf('3')]:
    a, b = Vnum(ev), Vnum(1/ev)
    print(f"eta={ev}: diff={mp.nstr(abs(a-b),3)}")
