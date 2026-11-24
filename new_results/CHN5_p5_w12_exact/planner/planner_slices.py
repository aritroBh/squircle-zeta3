# Planner (W4): independent check of CHN5, by elimination instead of curve tracing.
# (1) A = (1/100) L P8_A, CLASS = (1/432) L P8_CLASS (L = a2 a3 a4 (Sum a)) at fresh rational points, CHN3's original engine.
# (2) Ideal memberships: expand U*C4 + V*Q7' and compare with Q_B - (75/14) Q_A and Q_Z + (25/3) Q_A (exact).
# (3) COMPLETE enumeration of Gamma = {C4 = Q7' = P7 = 0} at fixed rational a1 (a3 = 1): flint resultants
#     R1 = Res_a5(C4, Q7'), R2 = Res_a5(C4, P7) in (a2, a4); factor both; for each factor pair, S(a2) = Res_a4(f, g);
#     all real positive roots a2 of S (flint complex_roots); lift to a4 and a5 as common roots; polish (C4, Q7', P7)
#     with findroot at 80 digits; keep points with a2, a4, a5 > 0; report the signs of P8_A and P8_CLASS there.
import sys, os, time, random
from fractions import Fraction as F
import sympy as sp, mpmath as mp, flint
HERE = '<project>/research/08_theory/CHN5_p5_w12_exact'
CHN3 = HERE + '/../CHN3_p5_locus'
OUT = HERE + '/planner/planner_slices.out'
def log(s):
    with open(OUT, 'a') as f: f.write(time.strftime('%H:%M:%S ') + s + '\n')
a = sp.symbols('a1:6'); a1, a2, a3, a4, a5 = a
rd = lambda fn: sp.sympify(open(fn).read())
C4, Q7, P7 = rd(HERE + '/C4.txt'), rd(HERE + '/Q7p.txt'), rd(HERE + '/P7.txt')
PA, PB, PC, PZ = (rd(f'{HERE}/P8_{k}.txt') for k in ('A', 'B', 'CLASS', 'Z'))
t0 = time.time()
# (2)
import pickle
for fn, lhs in (('ideal_B-75_14A.pkl', PB / 224 - sp.Rational(75, 14) * PA / 100), ('ideal_Z+25_3A.pkl', PZ / 3 + sp.Rational(25, 3) * PA / 100)):
    d = pickle.load(open(f'{HERE}/{fn}', 'rb'))
    rem = sp.expand(lhs - sp.sympify(d['U']) * C4 - sp.sympify(d['V']) * Q7)
    log(f'(2) {fn}: remainder == 0 ? {rem == 0} (if False, the quotient normalisation differs; nterms {len(sp.Add.make_args(rem))})')
# (1)
os.chdir(CHN3); sys.path.insert(0, CHN3)
import point_engine as PE
def coords(pt, W):
    vw = PE.volreg_words_point(pt, W); vwr = PE.volreg_words_point(pt[::-1], W)
    T = PE.table(W); acc = [F(0)] * T.nullity
    for src in (vw[W], vwr[W]):
        for w, pp in src.items():
            rw = T.reduce_word(w)
            for b in range(T.nullity):
                if rw[b]: acc[b] += pp * rw[b]
    return acc
random.seed(31)
for k in range(2):
    pt = tuple(F(random.randint(1, 50), random.randint(1, 17)) for _ in range(5))
    c = coords(pt, 12); L = pt[1] * pt[2] * pt[3] * sum(pt); sub = dict(zip(a, pt))
    okA = c[7] == L * F(str(PA.subs(sub))) / 100
    okC = c[9] + F(14, 9) * c[10] == L * F(str(PC.subs(sub))) / 432
    log(f'(1) fresh pt {k}: A formula {okA}, CLASS formula {okC} ({time.time()-t0:.0f}s)')
# (3)
ctx = flint.fmpq_mpoly_ctx.get(('a2', 'a4', 'a5'), 'lex')
def tofl(expr, a1v):
    p = sp.Poly(sp.expand(expr.subs({a1: a1v, a3: 1})), a2, a4, a5)
    return ctx.from_dict({m: flint.fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q)) for m, c in p.terms()})
def nontriv(fl):
    c, fac = fl.factor()
    return [f for f, e in fac if f.degrees()[0] + f.degrees()[1] > 0]
mp.mp.dps = 80
def run_slice(a1v):
    t = time.time()
    c4, q7, p7 = (tofl(e, a1v) for e in (C4, Q7, P7))
    R1, R2 = nontriv(c4.resultant(q7, 'a5')), nontriv(p7.resultant(c4, 'a5'))
    f = [sp.lambdify((a2, a4, a5), e.subs({a1: a1v, a3: 1}), 'mpmath') for e in (C4, Q7, P7, PA, PC)]
    cands = []
    for g1 in R1:
        for g2 in R2:
            S = g1.resultant(g2, 'a4')
            if S.is_zero(): log(f'  a1={a1v}: common factor between R1 and R2 factors (skipped pair)'); continue
            # univariate in a2 -> fmpz_poly
            dd = S.to_dict(); deg = max(m[0] for m in dd)
            coeffs = [flint.fmpq(0)] * (deg + 1)
            for m, v in dd.items(): coeffs[m[0]] += v
            up = flint.fmpq_poly(coeffs)
            for fac, e in up.factor()[1]:
                if fac.degree() < 1: continue
                for r, mult in fac.numer().complex_roots():
                    if abs(float(r.imag.mid())) > 1e-20 * max(1, abs(float(r.real.mid()))) or float(r.real.mid()) <= 0: continue
                    cands.append(mp.mpf(r.real.mid().str(70, radius=False)))
    pts = []
    for x2 in sorted(set(mp.nstr(x, 30) for x in cands)):
        x2 = mp.mpf(x2)
        g4 = sp.Poly(sp.expand(C4.subs({a1: a1v, a3: 1, a2: sp.Float(mp.nstr(x2, 70), 70)})), a4, a5)
        # a4 candidates: positive real roots of R1 factors at this a2
        for g1 in R1:
            poly4 = sp.Poly(sp.sympify(str(g1).replace('^', '**')).subs({sp.Symbol('a2'): sp.Float(mp.nstr(x2, 70), 70)}), sp.Symbol('a4'))
            if poly4.degree() < 1: continue
            for x4 in mp.polyroots([mp.mpf(str(c)) for c in poly4.all_coeffs()], maxsteps=400, extraprec=400):
                if abs(mp.im(x4)) > mp.mpf('1e-25') or mp.re(x4) <= 0: continue
                x4 = mp.re(x4)
                p5 = sp.Poly(sp.expand(C4.subs({a1: a1v, a3: 1, a2: sp.Float(mp.nstr(x2, 70), 70), a4: sp.Float(mp.nstr(x4, 70), 70)})), a5)
                for x5 in mp.polyroots([mp.mpf(str(c)) for c in p5.all_coeffs()], maxsteps=400, extraprec=400):
                    if abs(mp.im(x5)) > mp.mpf('1e-25') or mp.re(x5) <= 0: continue
                    x5 = mp.re(x5)
                    try:
                        sol = mp.findroot([lambda u, v, w: f[0](u, v, w), lambda u, v, w: f[1](u, v, w), lambda u, v, w: f[2](u, v, w)], (x2, x4, x5))
                    except Exception: continue
                    u, v, w = sol[0], sol[1], sol[2]
                    if min(u, v, w) <= 0: continue
                    sc = (1 + a1v + u + v + w) ** 8
                    if max(abs(f[0](u, v, w)), abs(f[1](u, v, w)), abs(f[2](u, v, w))) / sc > mp.mpf('1e-60'): continue
                    if abs(u - x2) + abs(v - x4) + abs(w - x5) > mp.mpf('1e-10') * (1 + abs(x2) + abs(x4) + abs(x5)): continue
                    key = (mp.nstr(u, 20), mp.nstr(v, 20), mp.nstr(w, 20))
                    if any(k == key for k, _ in pts): continue
                    pts.append((key, (mp.sign(f[3](u, v, w)), mp.sign(f[4](u, v, w)), f[3](u, v, w) / sc, f[4](u, v, w) / sc)))
    log(f'a1={a1v}: R1 factors {len(R1)}, R2 factors {len(R2)}, positive a2 candidates {len(cands)}, Gamma points {len(pts)} ({time.time()-t:.0f}s)')
    for key, (sa, sc_, na, nc) in pts:
        log(f'   (a2,a4,a5)=({key[0]}, {key[1]}, {key[2]})  sign P8_A {int(sa)}, sign P8_CLASS {int(sc_)}  P8_A/S^8 {mp.nstr(na, 4)}  P8_C/S^8 {mp.nstr(nc, 4)}')
    return pts
allpts = []
for a1v in (sp.Rational(1, 10**4), sp.Rational(1, 100), sp.Rational(1, 10), sp.Rational(1, 3), sp.Rational(2, 3), sp.Rational(3, 2),
            sp.Rational(2), sp.Rational(5), sp.Rational(20), sp.Rational(100), sp.Rational(10**3), sp.Rational(10**4)):
    try: allpts += run_slice(a1v)
    except Exception as e: log(f'a1={a1v}: FAILED {type(e).__name__}: {e}')
neg = sum(1 for _, s in allpts if s[0] < 0); pos = sum(1 for _, s in allpts if s[1] > 0)
log(f'SUMMARY: {len(allpts)} Gamma points; P8_A < 0 at {neg}; P8_CLASS > 0 at {pos} ({time.time()-t0:.0f}s)')
log('DONE')
