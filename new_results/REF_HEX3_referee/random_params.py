"""Table formula vs direct t-quadrature of the referee's own d^3 integrand (derived as in
integrand_check_fast.py), at random edge parameters not tied to any hexagon, incl. r up to 6.
Tests J-table + assembly + A-elimination as an identity in the free parameters."""
import sympy as sp, mpmath as mp, random, formula
t, X, Y, W, d, lA0, lA1, lB = sp.symbols('t X Y W d lA0 lA1 lB')
O = 4
def trunc(e):
    e = sp.expand(e); return sum(e.coeff(d, k)*d**k for k in range(O))
def mul(a_, b_): return trunc(a_*b_)
def ser_exp(z):
    out = 1; p = 1
    for k in range(1, O):
        p = mul(p, z); out += p/sp.factorial(k)
    return trunc(out)
def ser_log1p(z):
    out = 0; p = 1
    for k in range(1, O):
        p = mul(p, z); out += sp.Rational((-1)**(k+1), k)*p
    return trunc(out)
def own_E3(B, r, L, q0, V02):
    L2 = L**2; A1 = 1 - B; A0 = 1 - B*r; kap = r - 1
    a = [A0*t, A1*(1-t), B*(1+kap*t)]; ap = [A0, -A1, B*kap]; u = [lA0 + X, lA1 + Y, lB + W]
    E = [ser_exp(d*uj) for uj in u]
    F = trunc(sum(aj*Ej for aj, Ej in zip(a, E))); S = trunc(sum(apj*Ej for apj, Ej in zip(ap, E)))
    logF = ser_log1p(F - 1); Phi = mul(S, ser_exp(-logF))
    Fpow = ser_exp(mul(-2*logF, sum((-d)**k for k in range(O))))
    G = mul(Fpow, trunc(L2 - 2*(q0 + L2*t)*Phi + (V02 + 2*q0*t + L2*t**2)*mul(Phi, Phi)))
    sq = ser_exp(ser_log1p(trunc((G - L2)/L2))/2)
    return sp.expand(L*sq.coeff(d, 3))
mp.mp.dps = 45
random.seed(7); worst = 0
for i in range(10):
    B = sp.Rational(random.randint(1, 99), 100)
    rmax = min(600, int(99/B))
    r = sp.Rational(random.randint(1, rmax), 100)
    if i % 3 == 0: r = sp.Rational(1, 1)/B*sp.Rational(9, 10)   # large r near the limit A0 -> 0
    L = sp.Rational(random.randint(5, 40), random.randint(1, 9)); q0 = sp.Rational(random.randint(-30, 30), 7)
    V02 = (q0**2 + sp.Rational(random.randint(1, 300), 13))/L**2
    w = [sp.Rational(random.randint(1, 50), 10) for _ in range(3)]
    ed = dict(B=B, r=r, q0=q0, L2=L**2, V02=V02, A0=1-B*r, A1=1-B, wA0=w[0], wA1=w[1], wB=w[2])
    vals = {lA0: sp.log(ed['A0']/w[0]), lA1: sp.log(ed['A1']/w[1]), lB: sp.log(B/w[2])}
    e3 = own_E3(B, r, L, q0, V02)
    f = sp.lambdify([t, X, Y, W, lA0, lA1, lB], e3, 'mpmath')
    lv = [mp.mpf(sp.N(vals[s], 60)) for s in (lA0, lA1, lB)]
    kap = mp.mpf(sp.N(r - 1, 60))
    num = 6*mp.quad(lambda tt: f(tt, mp.log(tt), mp.log(1-tt), mp.log(1+kap*tt), *lv), [0, mp.mpf(1)/2, 1])
    form = formula.E_formula(ed)
    rel = abs(num-form)/max(1, abs(num)); worst = max(worst, rel)
    print(f'B={B} r={r} L={L} q0={q0}: quad {mp.nstr(num, 25)}  table-quad {mp.nstr(form-num, 3)}', flush=True)
print('worst relative', mp.nstr(worst, 3))
