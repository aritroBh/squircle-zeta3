"""g8c_structure.py -- EXACT sympy checks for Section 'The weight-three part as an iterated integral' (supplement, group g8c).
Own code (the notes' m2_structure.py was read, not imported).  All identities are in explicit rings:
  ring A: Q(c)(t) or Q(c)(g)  (c, t, g symbols);  ring B: Q(q)(y), c=(1-q^2)/(1+q^2), s=2q/(1+q^2), q=tan(alpha/4)=tau.
Steps
 [1] u^2-coefficient of the integrand, from the DEFINITION, with L=log g, m=log(1+g) as symbols (independent series in u);
     compare with the hand formula  RA=(1/g)[(R-gR'+g^2R'')/(2(1+g))-1/2], RB=(gR'-R)/(g(1+g)), RC=R/(2g(1+g)),
     R=sqrt(S), S=1+g^2-2kg, R'=(g-k)/R, R''=(1-k^2)/R^3.
 [2] Euler substitution t=g+R: g, R, 1+g, dg/dt in t;  RB^t, RC^t closed forms (hand);  RA^t.
 [3] KEY IDENTITY RA^t = V' + RA^s, V=(t-1)(t-k)/(2P), partial fractions of RA^s, RB^t, RC^t on the real letters, no polynomial part.
 [4] V*dlog g/dt = 1/(2(t+1)), V = 1/2 - s^2 (t+1)/P, g(t1)=1.
 [5] change of variable t=1+2 s^2 y: letters, g(y), 1+g(y), residue table in y   (ring B).
 [6] expansion of W3 integrand in the weight-one functions lambda_j (table of coefficients of  lambda_a lambda_b dlambda_j).
 [7] counts: 64 products of logarithms (10,11,13,15,15 by outer index), 23 diagonal / 41 off-diagonal, 105 ordered words, 21 ending in 0.
Needs no data file.
"""
import sympy as sp

t, c, g, L, m, u, y, q = sp.symbols('t c g L m u y q')
k = 2*c**2 - 1
s2 = 1 - c**2
P = t**2 - 2*k*t + 1
SQ = sp.Symbol('SQ', positive=True)

def red(e):                     # reduce SQ^2 -> S0 in a rational expression in SQ
    S0 = 1 + g**2 - 2*k*g
    e = sp.together(sp.expand(e))
    num, den = sp.fraction(e)
    def r(p):
        p = sp.Poly(sp.expand(p), SQ)
        out = 0
        for (n,), co in p.terms():
            out += co*S0**(n//2)*SQ**(n % 2)
        return sp.expand(out)
    return sp.cancel(r(num)/r(den))

# ---------------- [1] -----------------
S0 = 1 + g**2 - 2*k*g
Q = 1 + g**2*sp.exp(-2*u*L) - 2*k*g*sp.exp(-u*L)
d = sp.series((Q - S0)/S0, u, 0, 3).removeO()
sqrtQ = SQ*(1 + d/2 - d**2/8)
integ = sp.exp(-u*m)/(1 + g)*sqrtQ*sp.exp(u*L)/g - sp.exp(u*L)/g
c2 = sp.expand(sp.series(integ, u, 0, 3).removeO().coeff(u, 2))
poly = sp.Poly(c2, L, m)
print('[1] monomials (L,m) of the u^2 coefficient:', sorted(poly.monoms()))
RA, RB, RC = [poly.coeff_monomial(mon) for mon in (L**2, L*m, m**2)]
Rr = SQ
Rp = (g - k)/SQ
Rpp = (1 - k**2)/SQ**3
RA_h = (Rr - g*Rp + g**2*Rpp)/(2*(1 + g))/g - 1/(2*g)
RB_h = (g*Rp - Rr)/(g*(1 + g))
RC_h = Rr/(2*g*(1 + g))
for nm, a, b in (('RA', RA, RA_h), ('RB', RB, RB_h), ('RC', RC, RC_h)):
    print(f'    {nm}(series) - {nm}(hand) = ', sp.simplify(red(a - b)))
print("    check R'=(g-k)/R, R''=(1-k^2)/R^3 :",
      sp.simplify(sp.diff(sp.sqrt(S0), g) - (g - k)/sp.sqrt(S0)),
      sp.simplify(sp.diff(sp.sqrt(S0), g, 2) - (1 - k**2)/sp.sqrt(S0)**3))

# ---------------- [2] -----------------
g_t = (t**2 - 1)/(2*(t - k))
R_t = P/(2*(t - k))
print('[2] R_t^2 - S0(g_t) =', sp.simplify(R_t**2 - S0.subs(g, g_t)), '; g_t+R_t-t =', sp.simplify(g_t + R_t - t))
dg = sp.diff(g_t, t)
print('    dg/dt - P/(2(t-k)^2) =', sp.simplify(dg - P/(2*(t - k)**2)), '; dg/dt - R/(t-k) =', sp.simplify(dg - R_t/(t - k)))
print('    1+g_t - ((t+1)^2-4c^2)/(2(t-k)) =', sp.simplify(1 + g_t - ((t + 1)**2 - 4*c**2)/(2*(t - k))))
def to_t(e):
    return sp.factor(sp.cancel(red(e).subs(SQ, R_t).subs(g, g_t)*dg))
RAt, RBt, RCt = to_t(RA), to_t(RB), to_t(RC)
RBt_hand = 2*(k*t**2 - 2*t + k)/((t**2 - 1)*((t + 1)**2 - 4*c**2))
RCt_hand = P**2/(2*(t - k)*(t**2 - 1)*((t + 1)**2 - 4*c**2))
print('    RB^t - hand =', sp.simplify(RBt - RBt_hand), '; RC^t - hand =', sp.simplify(RCt - RCt_hand))
print('    RA^t denominator (factored):', sp.factor(sp.denom(sp.together(RAt))))

# ---------------- [3] -----------------
V = (t - 1)*(t - k)/(2*P)
Vp = sp.diff(V, t)
RAs = sp.cancel(RAt - Vp)
print('[3] V - (1/2 - s^2 (t+1)/P) =', sp.simplify(V - (sp.Rational(1, 2) - s2*(t + 1)/P)))
print('    V\' + s^2 d/dt[(t+1)/P] =', sp.simplify(Vp + s2*sp.diff((t + 1)/P, t)))
print('    RA^s denominator:', sp.factor(sp.denom(RAs)), ' (P is not a factor:',
      sp.rem(sp.Poly(sp.denom(RAs), t), sp.Poly(P, t)).is_zero is False, ')')
rho = {'1': 1, '-1': -1, 'k': k, 'r1': 2*c - 1, 'r2': -1 - 2*c}
for nm, R in (('RA^s', RAs), ('RB^t', RBt), ('RC^t', RCt)):
    res = {key: sp.factor(sp.cancel(sp.cancel(R*(t - r)).subs(t, r))) for key, r in rho.items()}
    rest = sp.cancel(R - sum(res[key]/(t - rho[key]) for key in rho))
    print(f'    {nm}: residues at t=1,-1,k,2c-1,-1-2c:', [str(res[key]) for key in rho], ' remainder after subtracting simple poles =', rest)

# ---------------- [4] -----------------
dlogg = sp.diff(sp.log(g_t), t)
print('[4] V*dlog g/dt - 1/(2(t+1)) =', sp.simplify(V*dlogg - 1/(2*(t + 1))))
s = sp.Symbol('s', positive=True)
t1 = 1 + 2*s
print('    g(t1) with s^2=1-c^2:', sp.simplify(g_t.subs(c**2, 1 - s**2).subs(t, t1).subs(c, sp.sqrt(1 - s**2))))
print('    V(1) =', V.subs(t, 1), ';  P(1) - 4 s^2 =', sp.simplify(P.subs(t, 1) - 4*(1 - c**2)))

# ---------------- [5] ring B -----------------
cq = (1 - q**2)/(1 + q**2)
sq = 2*q/(1 + q**2)
print('[5] c^2+s^2-1 =', sp.simplify(cq**2 + sq**2 - 1), '; (1-c)/s - q =', sp.simplify((1 - cq)/sq - q))
tq = 1 + 2*sq**2*y                       # t = 1 + 2 s^2 y
kq = 2*cq**2 - 1
gy = sp.cancel(g_t.subs(c, cq).subs(t, tq).subs(t, tq) if False else ((tq**2 - 1)/(2*(tq - kq))))
print('    g(y) - y(1+s^2 y)/(1+y) =', sp.simplify(gy - y*(1 + sq**2*y)/(1 + y)))
onepg = sp.cancel(1 + gy)
print('    1+g(y) - (1+(1+c)y)(1+(1-c)y)/(1+y) =', sp.simplify(onepg - (1 + (1 + cq)*y)*(1 + (1 - cq)*y)/(1 + y)))
for key, r in rho.items():
    rq = r.subs(c, cq) if hasattr(r, 'subs') else r
    a = sp.simplify((rq - 1)/(2*sq**2))
    print(f'    letter t={key}: y = (rho-1)/(2s^2) =', sp.factor(a))
print('    1/s^2 = ', sp.factor(1/sq**2), ';  -1/(1+c) =', sp.factor(-1/(1 + cq)), ';  -1/(1-c) =', sp.factor(-1/(1 - cq)))
print('    y(t1) = 1/s: g(1/s) =', sp.simplify(gy.subs(y, 1/sq)))
print('    dt = 2 s^2 dy, t - rho = 2 s^2 (y - a_rho):  checked by the line above.')

# ---------------- [6] -----------------
lam = sp.symbols('l0:5')          # lambda_0=log y, lambda_1=log(1+y), lambda_2=log(1+s^2 y), lambda_3=log(1+(1+c)y), lambda_4=log(1+(1-c)y)
Lg = lam[0] + lam[2] - lam[1]
Lm = lam[3] + lam[4] - lam[1]
qq = (c**2 + 1)/(4*c)
tab = {0: (0, -1, sp.Rational(1, 2)), 2: (-1, 1, -sp.Rational(1, 2)), 1: (sp.Rational(1, 2), 0, sp.Rational(1, 2)),
       3: (-qq, c, -c), 4: (qq, -c, c)}
print('[6] letters a_j: a_0=0, a_1=-1, a_2=-1/(1-c^2), a_3=-1/(1+c), a_4=-1/(1-c); table (A_j,B_j,C_j):')
for j in range(5):
    print('    j =', j, tab[j])
coef = {}
for j, (A, B, C) in tab.items():
    ex = sp.expand(A*Lg**2 + B*Lg*Lm + C*Lm**2)
    for mon, co in sp.Poly(ex, *lam).terms():
        if co != 0:
            idx = []
            for i, e in enumerate(mon):
                idx += [i]*e
            coef[(tuple(idx), j)] = sp.factor(co)
print('    number of terms  coeff * int_0^Y lambda_a lambda_b dlambda_j :', len(coef))
for (ab, j), co in sorted(coef.items()):
    print('      (a,b)=', ab, ' j=', j, ' coeff =', co)
print('    letters used as the outer form (j):', sorted({j for (_, j) in coef}))

# ---------------- [7] counts: products of logarithms and ordered words -----------------
# lambda_a lambda_b = I(a,b)+I(b,a) (a != b), = 2 I(a,a) (a = b); so  int lambda_a lambda_b dlambda_j = I(j,a,b)+I(j,b,a)  or  2 I(j,a,a).
from collections import Counter
print('[7] counts')
print('    products of two logarithms against one form (nonzero coefficient):', len(coef))
print('    by outer index j = 0..4:', [sum(1 for (_, jj) in coef if jj == j) for j in range(5)])
diag = [(ab, j) for (ab, j) in coef if ab[0] == ab[1]]
off = [(ab, j) for (ab, j) in coef if ab[0] != ab[1]]
print('    diagonal (a=b):', len(diag), '; off-diagonal (a<b):', len(off))
words = {}
for (ab, j), co in coef.items():
    a, b = ab
    if a == b:
        words[(j, a, a)] = words.get((j, a, a), 0) + 2*co
    else:
        words[(j, a, b)] = words.get((j, a, b), 0) + co
        words[(j, b, a)] = words.get((j, b, a), 0) + co
words = {w: sp.factor(v) for w, v in words.items() if sp.simplify(v) != 0}
print('    distinct ordered words (j,a,b) with nonzero coefficient:', len(words), '  = diagonal + 2*off-diagonal =', len(diag) + 2*len(off))
print('    uncollected summands (diagonal counted twice):', len(diag)*2 + 2*len(off))
print('    words ending in 0 (innermost letter 0):', sum(1 for w in words if w[2] == 0))
print('    words (j,0,0):', sorted(w for w in words if w[1] == 0 and w[2] == 0))
print('    the pure word (0,0,0) occurs:', (0, 0, 0) in words)
print('    words (0,a,b) (j = 0) that contain the letter 0:', sorted(w for w in words if w[0] == 0 and 0 in w[1:]))
print('    every coefficient times c is a polynomial in c of degree <= 2 (so the coefficients lie in Q + Qc + Q/c):',
      all(sp.Poly(sp.cancel(v*c), c).degree() <= 2 for v in words.values()))
