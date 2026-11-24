"""check_symbols.py -- (1) d,m,p are genuine sympy Symbols, (2) free symbols / size of the identities,
(3) pointwise exact re-evaluation of the identities at rational points, (4) negative controls (perturbed
e-family must FAIL the certificates)."""
import sys, time, random
import sympy as sp
from fractions import Fraction
sys.path.insert(0, '.')
import tow3lib as T3
import prove_V56 as PV
import prove_adj56 as PA
import prove_E3b56 as PE

d, m = T3.d, T3.m
p = PA.p
print('--- symbols ---')
for name, s in (('tow3lib.d', d), ('tow3lib.m', m), ('prove_adj56.p', p), ('prove_E3b56.p', PE.p)):
    print(f'{name}: type={type(s).__name__}, is_Symbol={s.is_Symbol}, assumptions0={dict(s.assumptions0)}, free_symbols={s.free_symbols}')
print('prove_adj56.d is the expression', PA.d, '; free symbols', PA.d.free_symbols)
print('prove_V56 uses d,m from tow3lib:', PV.d is T3.d, PV.m is T3.m)

print('--- identity sizes / free symbols (non-adjacent) ---')
for (dl, s) in [(5,0),(5,1),(5,2),(6,0),(6,1),(6,2),(6,3)]:
    lhs = sp.expand(T3.e(dl, s) + T3.e(dl, s+1))
    rhs = T3.Lam(dl, s)
    print(f'({dl},{s}): lhs free={sorted(map(str,lhs.free_symbols))} ops={sp.count_ops(lhs)} | rhs free={sorted(map(str,rhs.free_symbols))} ops={sp.count_ops(rhs)}')
print('--- adjacent: free symbols of feed terms ---')
for (dl, s) in [(5,2),(6,3)]:
    terms = PA.feed_terms(dl, s)
    fs = set()
    for t in terms: fs |= t.free_symbols
    lhs = sp.expand(PA.eR(dl, s, PA.d, p) + PA.eR(dl, s+1, PA.d, p))
    print(f'adj ({dl},{s}): {len(terms)} feed terms, free symbols of terms={sorted(map(str,fs))}, lhs free={sorted(map(str,lhs.free_symbols))}')

print('--- pointwise exact evaluation, non-adjacent, random rational (d,m) incl. non-integer ---')
random.seed(7)
pts = [(Fraction(random.randint(9, 60), random.randint(1, 7)), Fraction(random.randint(1, 40), random.randint(1, 7))) for _ in range(4)] + [(sp.Integer(11), sp.Integer(6))]
for (dl, s) in [(5,0),(5,1),(5,2),(6,0),(6,1),(6,2),(6,3)]:
    lhs = sp.expand(T3.e(dl, s) + T3.e(dl, s+1))
    rhs = T3.Lam(dl, s)
    resid = []
    for (dv, mv) in pts:
        dv = sp.Rational(dv); mv = sp.Rational(mv)
        r = sp.simplify((lhs - rhs).subs({d: dv, m: mv}))
        resid.append(r)
    print(f'({dl},{s}) residual at {len(pts)} points:', resid)

print('--- negative controls (perturbed e-family must break the identity) ---')
def ident(dl, s):
    return PV.identity(dl, s)
# perturb a_5
base5 = T3.arho(5); T3._a[5] = sp.expand(base5 + 1)
print('perturb a_5 -> V(5,0) numerator==0 ?', ident(5, 0)); T3._a[5] = base5
base6 = T3.arho(6); T3._a[6] = sp.expand(base6 + 1)
print('perturb a_6 -> V(6,0) numerator==0 ?', ident(6, 0)); T3._a[6] = base6
# perturb e51
e51 = T3._EF[(5,1)]; T3._EF[(5,1)] = e51 + d
print('perturb e(5,1) -> V(5,1) numerator==0 ?', ident(5, 1), '| V(5,0):', ident(5, 0)); T3._EF[(5,1)] = e51
e61 = T3._EF[(6,1)]; T3._EF[(6,1)] = e61 + m
print('perturb e(6,1) -> V(6,1) numerator==0 ?', ident(6, 1), '| V(6,0):', ident(6, 0)); T3._EF[(6,1)] = e61
e62 = T3._EF[(6,2)]; T3._EF[(6,2)] = e62 + 1
print('perturb e(6,2) -> V(6,2) numerator==0 ?', ident(6, 2), '| V(6,1):', ident(6, 1)); T3._EF[(6,2)] = e62
A4 = T3.A[4]; T3.A[4] = A4 + 1
print('perturb zigzag A_4 5->6 -> V(6,3) numerator==0 ?', ident(6, 3)); T3.A[4] = A4
A3 = T3.A[3]; T3.A[3] = A3 + 1
print('perturb zigzag A_3 2->3 -> V(5,2) numerator==0 ?', ident(5, 2)); T3.A[3] = A3
print('restored: V(5,1), V(6,3) numerator==0 ?', ident(5, 1), ident(6, 3))

print('--- pointwise adjacent identities at integer p (unperturbed) and with perturbed e(5,1), a_5 ---')
def adj_resid(dl, s, pv):
    lhs = sp.expand(PA.eR(dl, s, PA.d, p) + PA.eR(dl, s+1, PA.d, p))
    terms = PA.feed_terms(dl, s)
    tot = lhs - sum(terms)
    return sp.simplify(tot.subs(p, pv))
t0 = time.time()
for (dl, s) in [(5,2),(6,3),(5,1),(6,2)]:
    print(f'adj ({dl},{s}) residual at p=3,4,5,8,13/2:', [adj_resid(dl, s, pv) for pv in (3, 4, 5, 8, sp.Rational(13, 2))], f'({round(time.time()-t0,1)}s)', flush=True)
e51 = T3._EF[(5,1)]; T3._EF[(5,1)] = e51 + d
print('perturbed e(5,1): adj (5,1) residual at p=3,5:', [adj_resid(5, 1, pv) for pv in (3, 5)], flush=True)
T3._EF[(5,1)] = e51
b5 = T3.arho(5); T3._a[5] = sp.expand(b5 + 1)
print('perturbed a_5: adj (5,2) residual at p=3,5 (top, uses s=2 only):', [adj_resid(5, 2, pv) for pv in (3, 5)], flush=True)
T3._a[5] = b5
A3 = T3.A[3]; T3.A[3] = A3 + 1
print('perturbed A_3: adj (5,2) residual at p=3,5:', [adj_resid(5, 2, pv) for pv in (3, 5)], flush=True)
T3.A[3] = A3
print('--- E3 negative control (perturbed e(5,1)): feed(5) numerator zero? ---')
e51 = T3._EF[(5,1)]; T3._EF[(5,1)] = e51 + d
f = PE.feed(5); num = sp.expand(sp.fraction(sp.together(f))[0])
print('perturbed e(5,1): E3 delta=5 numerator == 0 ?', num == 0)
T3._EF[(5,1)] = e51
print('done', round(time.time()-t0, 1))
