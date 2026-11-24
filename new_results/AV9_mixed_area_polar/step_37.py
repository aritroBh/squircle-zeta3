
fix = r'''
import mpmath as mp, sympy as sp, av9core as C
h, e = C.selfdual_coeffs(4)
for n in [1,3]:
    en = sp.expand(e[n])
    deg_n = sum(cf*C.lv**a*C.lw**b for (a,b),cf in sp.Poly(en, C.lv, C.lw).terms() if a+b==n)
    pred = (C.lv**n + C.lw**n)/(2*sp.factorial(n))
    print(f"n={n}: deg_n =", sp.expand(deg_n))
    print(f"     pred  =", sp.expand(pred))
    print(f"     ratio:", sp.simplify(deg_n/pred))
# P_n(1,1) vs 2 P_n(1,0) for odd n
print("\nP_n(1,1) and 2P_n(1,0):")
for n in range(1,9):
    en = sp.expand(e[n])
    deg_n = sum(cf*C.lv**a*C.lw**b for (a,b),cf in sp.Poly(en, C.lv, C.lw).terms() if a+b==n)
    p11 = deg_n.subs({C.lv:1, C.lw:1}); p10 = deg_n.subs({C.lv:1, C.lw:0})
    print(f"n={n}: P(1,1)={p11}  2P(1,0)={2*p10}  equal={sp.simplify(p11-2*p10)==0}")

# duality with careful quadrature
mp.mp.dps = 80
def Vnum(etaval):
    q = 1+mp.mpf(etaval)
    f = lambda v: (2/q)*(v*(1-v))**(1/q-1)*(v**(q-1)+(1-v)**(q-1))**(1/q)
    return mp.quad(f, [mp.mpf('1e-12'), mp.mpf('1e-6'), mp.mpf('1e-3'), mp.mpf('0.5')])
for ev in ['0.25','0.1','0.05']:
    a = Vnum(ev); b = Vnum(mp.mpf(1)/mp.mpf(ev))
    print(f"eta={ev}: V={mp.nstr(a,20)} V(1/eta)={mp.nstr(b,20)} diff={mp.nstr(abs(a-b),3)}")
'''
open('fix2.py','w').write(fix)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','fix2.py'],capture_output=True,text=True,timeout=900)
print(r.stdout); print(r.stderr[-300:])
