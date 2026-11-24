"""step_07.py -- the pi^4 part of mu_{4,n} (z4 and z2^2 parts): EXACT rational functions, PROVED.

From the symbolic formula (mu4_formula.pkl): the z4-coefficient and the z2^2-coefficient are pure
Q(n) elements (z4 and z2^2 arise only from kappa_4 and kappa_2^2, never multiplied by harmonic
factors: kappa_4*kappa_1 would have total degree 5 > 4).
Extract them and identify:
  [z4]  mu_{4,n} = (n+3)/((2n+3)(4n^2-1))         (PROVED, rational-function identity)
  [z2^2]mu_{4,n} = -n/((2n+3)(4n^2-1))            (PROVED)
Generating functions (artanh c = (H_1 + H_{-1})/2):
  F1(c) = sum_{n>=1} (n+3)/((2n+3)(4n^2-1)) c^{2n},  F2(c) = sum_{n>=1} -n/((2n+3)(4n^2-1)) c^{2n}
are elementary by partial fractions in n: sum_{n>=1} c^{2n}/(2n-1) = c artanh c,
sum c^{2n}/(2n+1) = (artanh c - c)/c, sum c^{2n}/(2n+3) = (artanh c - c - c^3/3)/c^3  (all termwise exact).
Check: the resulting series agree with the exact table through c^1040 (EXACT-1040; PROVED for all n
by the partial-fraction argument).
"""
import json, pickle
from fractions import Fraction as Fr
import sympy as sp
import sys
sys.path.insert(0, "<project>/research/10_frontier/W1_corner_family/B_taylor_pi4")

n = sp.symbols('n')
terms = pickle.load(open('mu4_formula.pkl', 'rb'))
c4 = sp.Rational(0)   # coefficient of z4
c22 = sp.Rational(0)  # coefficient of z2^2
rest = 0
for mon, cf in terms:
    e2, e3, e4 = mon[8], mon[9], mon[10]
    if e4 == 1:
        assert all(e == 0 for e in mon[:8]) and e2 == 0 and e3 == 0
        c4 += sp.sympify(cf)
    elif e2 == 2:
        assert all(e == 0 for e in mon[:8]) and e3 == 0 and e4 == 0
        c22 += sp.sympify(cf)
    else:
        rest += 1
t1 = (n + 3)/((2*n + 3)*(4*n**2 - 1))
t2 = -n/((2*n + 3)*(4*n**2 - 1))
print('monomials not containing z4 or z2^2:', rest)
print('[z4]  coefficient extracted:', sp.factor(c4))
print('[z2^2] coefficient extracted:', sp.factor(c22))
print('PROVED [z4]  == (n+3)/((2n+3)(4n^2-1)) :', sp.simplify(c4 - t1) == 0)
print('PROVED [z2^2] == -n/((2n+3)(4n^2-1))  :', sp.simplify(c22 - t2) == 0)

# partial fractions
a1 = sp.apart(t1, n)
a2 = sp.apart(t2, n)
print('apart [z4]:', a1)
print('apart [z2^2]:', a2)

# closed forms  F = A G_{-1} + B G_1 + C G_3,  G_m as in the docstring
c = sp.symbols('c')
art = sp.atanh(c)
G = {-1: c*art, 1: (art - c)/c, 3: (art - c - c**3/3)/c**3}
def closed(apart_expr):
    F = 0
    for term in sp.Add.make_args(apart_expr):
        num, den = sp.fraction(term)
        dp = sp.Poly(sp.expand(den), n)
        assert dp.degree() == 1
        alpha = sp.Rational(dp.coeff_monomial(n))/2     # den = alpha*(2n + m)
        m = int(sp.Rational(dp.coeff_monomial(1))/alpha)
        A = sp.Rational(num)/alpha
        F += A*G[m]
    return sp.expand(F)
F1 = closed(a1)
F2 = closed(a2)
print('F1(c) =', F1)
print('F2(c) =', F2)

# exact check through c^1040 vs table
ORDER = 1040
s1 = sp.series(F1, c, 0, ORDER + 1).removeO().expand()
s2 = sp.series(F2, c, 0, ORDER + 1).removeO().expand()
tab = json.load(open('flat_exact_N520_K4_formula.json'))['4']
ok = True
for nn in range(1, 521):
    v1 = Fr(tab[str(nn)].get('(0, 0, 1)', '0'))
    v2 = Fr(tab[str(nn)].get('(2, 0, 0)', '0'))
    if Fr(s1.coeff(c, 2*nn)) != v1 or Fr(s2.coeff(c, 2*nn)) != v2:
        ok = False; print('mismatch at n=', nn)
print('EXACT through c^1040 (z4 and z2^2 parts):', ok)

# HPL form: artanh c = (H_1 + H_{-1})/2. Save as c^e H_w terms.
def to_hpl_terms(F):
    F = sp.expand(F)
    out = []
    for term in sp.Add.make_args(F):
        if term.has(art):
            co = sp.simplify(term/art)
            out.append((str(sp.factor(co)), 'artanh'))
        else:
            out.append((str(sp.factor(term)), '1'))
    return out
json.dump({'F1': to_hpl_terms(F1), 'F2': to_hpl_terms(F2)}, open('b4_pi4_part.json', 'w'), indent=1)
print('saved b4_pi4_part.json')
