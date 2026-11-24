# 02_P2c.py -- P2c = int_0^1 arsinh(1/q)/sqrt(1+q^2) ds, q = sqrt(s^2+(1-s)^2).
# Exact substitutions (sympy), exact antiderivative check, boundary values by named classical identities
# with explicit branches (each instance checked to 80 digits), exact assembly.
import sympy as sp, mpmath as mp, json
from common import lam, l3, L, al, PI, Li2, li2, num, show, B as Bq
mp.mp.dps = 90
ok = lambda e: sp.simplify(e) == 0
s, q, x, be = sp.symbols('s q x beta', positive=True)
r2, r3 = sp.sqrt(2), sp.sqrt(3)
c = r2/r3                                  # cos(alpha_c)
ref = json.load(open("reference_100.json"))

print("== Step 1: substitutions")
# (i) s in (0,1/2) <-> q in (1/sqrt2, 1):  s = (1 - sqrt(2q^2-1))/2
s_of_q = (1 - sp.sqrt(2*q**2 - 1))/2
assert ok(s_of_q**2 + (1 - s_of_q)**2 - q**2)
assert ok(sp.diff(s_of_q, q) + q/sp.sqrt(2*q**2 - 1))
print("  s=(1-sqrt(2q^2-1))/2: s^2+(1-s)^2=q^2, ds/dq=-q/sqrt(2q^2-1)   (exact)")
# (ii) q <-> x in (1, sqrt3):  beta = sqrt(1+q^2) = (x^2+1)/(2 c x),  q = sqrt(beta^2-1)
bx = (x**2 + 1)/(2*c*x)
qx = sp.sqrt(bx**2 - 1)
# q(x) as an explicit positive expression: beta^2 - 1 = (3x^4 - 2x^2 + 3)/(8x^2) > 0
assert ok(bx**2 - 1 - (3*x**4 - 2*x**2 + 3)/(8*x**2))
qx = sp.sqrt(3*x**4 - 2*x**2 + 3)/(2*r2*x)
# 2 q^2 - 1 = 3 (x^2-1)^2/(4 x^2)  => sqrt(2q^2-1) = sqrt3 (x^2-1)/(2x) for x>1
assert ok(2*qx**2 - 1 - 3*(x**2 - 1)**2/(4*x**2))
sq = r3*(x**2 - 1)/(2*x)
# q dq = beta dbeta (q^2 = beta^2 - 1), so 2 q q'(x)/(beta sqrt(2q^2-1)) = 2 beta'(x)/sqrt(2q^2-1)
assert ok(sp.diff(qx**2, x) - sp.diff(bx**2, x))
jac = sp.simplify(2*sp.diff(bx, x)/sq)
assert ok(jac - r2/x)
print("  beta=(x^2+1)/(2cx), c=sqrt(2/3): 2 q q'(x)/(beta sqrt(2q^2-1)) = sqrt2/x   (exact)")
# arsinh(1/q) = log((1+beta)/q) and ((1+beta)/q)^2 = (beta+1)/(beta-1) = (x^2+2cx+1)/(x^2-2cx+1)
assert ok(((1 + bx)/qx)**2 - (x**2 + 2*c*x + 1)/(x**2 - 2*c*x + 1))
assert ok(bx.subs(x, 1) - sp.sqrt(sp.Rational(3, 2))) and ok(bx.subs(x, r3) - r2)
assert ok(qx.subs(x, 1) - 1/r2) and ok(qx.subs(x, r3) - 1)
print("  ((1+beta)/q)^2 = (x^2+2cx+1)/(x^2-2cx+1); x=1 <-> q=1/sqrt2, x=sqrt3 <-> q=1   (exact)")
print("  arsinh(1/q) = (1/2) log((x^2+2cx+1)/(x^2-2cx+1)), so")
print("  => P2c = (1/sqrt2) int_1^sqrt3 [log(x^2+2cx+1) - log(x^2-2cx+1)] dx/x")
num_x = mp.quad(lambda t: (mp.log(t*t + 2*mp.sqrt(mp.mpf(2)/3)*t + 1) - mp.log(t*t - 2*mp.sqrt(mp.mpf(2)/3)*t + 1))/t,
                [1, mp.sqrt(3)])/mp.sqrt(2)
show("P2c: x-form (quadrature) vs reference", num_x, mp.mpf(ref["P2c"]), 80)

print("== Step 2: antiderivative G(x) = sqrt2 [Re Li2(x w) - Re Li2(-x w)],  w = e^{i alpha_c}")
w = sp.symbols('w')
d = sp.expand_func(sp.diff(sp.polylog(2, x*w), x))   # polylog(1,z) = -log(1-z)
assert ok(d - (-sp.log(1 - x*w)/x))
print("  d/dx Li2(x w) = -log(1 - x w)/x  (exact; Li2(xw) analytic since xw notin [1,oo) for x>0, w = e^{i alpha})")
wv = (r2 + sp.I)/r3
assert ok(sp.cos(sp.atan(1/r2)) - c) and ok(sp.sin(sp.atan(1/r2)) - 1/r3) and ok(wv - (c + sp.I/r3))
assert ok(sp.expand((1 - x*wv)*(1 - x*sp.conjugate(wv))) - (x**2 - 2*c*x + 1))
assert ok(sp.expand((1 + x*wv)*(1 + x*sp.conjugate(wv))) - (x**2 + 2*c*x + 1))
print("  |1 -+ x w|^2 = x^2 -+ 2cx + 1  (exact); Re log(1-xw) = log|1-xw|  ==>  G'(x) = integrand")
for xx in [mp.mpf('1.1'), mp.mpf('1.5')]:
    wn = mp.expj(mp.atan(1/mp.sqrt(2)))
    G = lambda t: mp.sqrt(2)*(mp.re(mp.polylog(2, t*wn)) - mp.re(mp.polylog(2, -t*wn)))
    gd = mp.diff(G, xx)
    integ = (mp.log(xx*xx + 2*mp.sqrt(mp.mpf(2)/3)*xx + 1) - mp.log(xx*xx - 2*mp.sqrt(mp.mpf(2)/3)*xx + 1))/(mp.sqrt(2)*xx)
    show("  G'(x) numerically at x=%s" % mp.nstr(xx, 3), gd, integ, 60)
print("  => P2c = G(sqrt3) - G(1) = 2 sqrt2 [Re chi2(sqrt2+i) - Re chi2(e^{i alpha})],  chi2(z) = (Li2(z)-Li2(-z))/2")

print("== Step 3: boundary values (unknowns X1=Re chi2(sqrt2+i), X2=Re chi2((sqrt2-i)/3), X3=Re chi2(y), X4=Re Li2(y))")
X0, X1, X2, X3, X4 = sp.symbols('X0:5')
eqs = []
# (a) Re Li2(e^{i th}) = pi^2/6 - th(2pi-th)/4, 0<=th<=2pi  (Fourier series sum cos(n th)/n^2)
ReLiU = lambda th: PI**2/6 - th*(2*PI - th)/4
X0val = sp.expand((ReLiU(al) - ReLiU(al + PI))/2)          # -e^{i alpha} = e^{i(alpha+pi)}, alpha+pi in [0,2pi]
print("  (a) Re chi2(e^{i alpha}) =", X0val)
aN = mp.atan(1/mp.sqrt(2))
show("  (a) numeric", (mp.re(mp.polylog(2, mp.expj(aN))) - mp.re(mp.polylog(2, -mp.expj(aN))))/2,
     mp.pi**2/8 - mp.pi*aN/4, 80)
# (b) inversion: chi2(z) + chi2(1/z) = (log^2 z - log^2(-z))/4, valid for z notin R.  z = sqrt2+i, 1/z = (sqrt2-i)/3
z = r2 + sp.I
assert ok(1/z - (r2 - sp.I)/3)
# log z = (1/2)log3 + i alpha;  log(-z) = (1/2)log3 + i(alpha - pi)   (principal: alpha-pi in (-pi,pi])
re_b = sp.expand(((al - PI)**2 - al**2)/4)                  # Re[(p+i t1)^2 - (p+i t2)^2] = t2^2 - t1^2
eqs.append(sp.Eq(X1 + X2, re_b))                             # Re chi2(1/z) = Re chi2(conj) = X2 (conjugate symmetry)
zN = mp.sqrt(2) + 1j
chiN = lambda t: (mp.polylog(2, t) - mp.polylog(2, -t))/2
show("  (b) inversion chi2 at z=sqrt2+i (real part)", mp.re(chiN(zN) + chiN(1/zN)), re_b.subs(al, sp.atan(1/r2)).evalf(90), 80)
# (c) Landen for chi2: chi2(t) + chi2((1-t)/(1+t)) = pi^2/8 + (1/2) log t log((1+t)/(1-t)),
#     valid on C \ ((-oo,0] U [1,oo)) (identity theorem from 0<t<1).  t = (sqrt2-i)/3, y = (1-t)/(1+t) = b e^{i pi/4}
t = (r2 - sp.I)/3
y = (1 - t)/(1 + t)
bb = r2 - 1; aa = 1 + r2
assert ok(sp.expand(y*(1 + t)) - sp.expand(1 - t))
assert ok(sp.radsimp(y) - bb*(1 + sp.I)/r2)
print("  (c) y = (1-t)/(1+t) = b e^{i pi/4} (exact);  log t = -(1/2)log3 - i alpha;  log((1+t)/(1-t)) = log(a e^{-i pi/4}) = L - i pi/4")
re_c = sp.expand(PI**2/8 + sp.Rational(1, 2)*((-l3/2)*L - (-al)*(-PI/4)))
eqs.append(sp.Eq(X2 + X3, re_c))
tN = (mp.sqrt(2) - 1j)/3
show("  (c) Landen chi2 at t=(sqrt2-i)/3 (full complex identity)", abs(chiN(tN) + chiN((1 - tN)/(1 + tN)) -
     (mp.pi**2/8 + mp.log(tN)*mp.log((1 + tN)/(1 - tN))/2)), 0, 80)
show("  (c) real part of RHS", mp.re(mp.pi**2/8 + mp.log(tN)*mp.log((1 + tN)/(1 - tN))/2), num(re_c), 80)
# (d) duplication Li2(y)+Li2(-y) = Li2(y^2)/2 (|y|<1) => chi2(y) = Li2(y) - Li2(y^2)/4;  y^2 = i b^2,
#     Re Li2(i b^2) = sum_m (-1)^m b^{4m}/(4m^2) = Li2(-b^4)/4
assert ok(sp.radsimp(sp.expand((bb*(1 + sp.I)/r2)**2)) - sp.I*bb**2)          # y^2 = i b^2
eqs.append(sp.Eq(X3, X4 - li2(-(Bq**4))/16))
yN = (mp.sqrt(2) - 1)*mp.expjpi(mp.mpf(1)/4)
show("  (d) Re chi2(y) = Re Li2(y) - Li2(-b^4)/16", mp.re(chiN(yN)), mp.re(mp.polylog(2, yN)) - mp.polylog(2, -(mp.sqrt(2) - 1)**4)/16, 80)
# (e) Lemma A (supplement, Lemma g3:lem:A): Re Li2(b e^{i pi/4}) = Li2(b)/2 + pi^2/192
eqs.append(sp.Eq(X4, li2(Bq)/2 + PI**2/192))
show("  (e) Lemma A", mp.re(mp.polylog(2, yN)), mp.polylog(2, mp.sqrt(2) - 1)/2 + mp.pi**2/192, 80)
sol = sp.solve(eqs, [X1, X2, X3, X4], dict=True)[0]
P2c = sp.expand(2*r2*(sol[X1] - X0val))
print("== Result")
print("  P2c =", P2c)
target = r2*(li2(Bq) - li2(-(Bq**4))/8 + L*l3/2 + PI**2/96 - PI*al/4)
dif = sp.expand(P2c - sp.expand(target))
print("  P2c - [Lemma P closed form]  (exact) =", dif)
assert dif == 0
show("P2c exact value vs 100-digit quadrature", num(P2c), mp.mpf(ref["P2c"]), 90)
json.dump({"P2c": sp.srepr(P2c)}, open("P2c_exact.json", "w"))
