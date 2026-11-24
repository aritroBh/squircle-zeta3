# Main-session planner check of the minimised D = 17 form (own code): residual of
# T_coef*T + z5*zeta(5) + sum c_x Li5(x) + sum d_p prod(letters) with T = L(chi17,5)/sqrt17.
import json, mpmath as mp, sympy as sp
mp.mp.dps = 110
d = json.load(open('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S/d17_final.json'))
s = mp.sqrt(17)
def val(expr):  # "a+b*sqrt17" with rational a, b
    e = sp.sympify(expr.replace('sqrt17', 'sqrt(17)'))
    return mp.mpf(sp.Rational(e.subs(sp.sqrt(17), 0))) + s*mp.mpf(sp.Rational(sp.expand(e).coeff(sp.sqrt(17))))
leg = {a: (1 if pow(a, 8, 17) == 1 else -1) for a in range(1, 17)}
L = sum(leg[a]*mp.zeta(5, mp.mpf(a)/17) for a in range(1, 17))/mp.mpf(17)**5
T = L/s
lv = {'lE17': mp.log(4 + s), 'lp': mp.log((3 + s)/2), "lp'": mp.log(abs((3 - s)/2)), 'pi2': mp.pi**2}
tot = mp.mpf(int(d['T_coef']))*T + mp.mpf(int(d['z5_coef']))*mp.zeta(5)
mx = 0
for x, c in d['li5'].items():
    xv = val(x); assert abs(xv) < 1, x
    q = sp.Rational(c); tot += mp.mpf(q.p)/q.q*mp.polylog(5, xv); mx = max(mx, abs(q.p))
for k, c in d['products'].items():
    p = mp.mpf(1)
    for f in k.split('*'): p *= lv[f]
    q = sp.Rational(c); tot += mp.mpf(q.p)/q.q*p; mx = max(mx, abs(q.p))
print("L(chi17,5) =", mp.nstr(L, 30)); print("terms:", len(d['li5']), "Li5 +", len(d['products']), "products; max |coefficient| =", mx)
print("residual at 110 digits:", mp.nstr(tot, 5))
