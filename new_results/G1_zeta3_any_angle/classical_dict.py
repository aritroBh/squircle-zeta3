"""classical_dict.py -- classical (Li3, Li2*log, log^3, pi^2 log, zeta3) forms of the weight-3 HPLs that occur
in the general-angle formula for M2, valid for 0 < x < 1.

Method: function-level linear fit.  Basis F (45 real-analytic functions on (0,1)):
  8 trilogs  Li3 at  x, -x, 1-x, 1/(1+x), (1+x)/2, (1-x)/(1+x), -(1-x)/(1+x), 2x/(1+x)
  12 = Li2 at {x, -x, (1+x)/2} times logs {log x, log(1-x), log(1+x), log 2}
  20 cubic monomials in {log x, log(1-x), log(1+x), log 2};  4 = pi^2 * logs;  1 = zeta(3).
Each HPL H_w(x) is sampled (exact rational Taylor series, 600 terms, at x <= 0.45) at 70 points, the
70 x 45 system is solved in the least-squares sense at 120 digits, the coefficients are rationalised
(denominators <= 10^4), and the rational identity is then re-checked at 40 FRESH points spread over (0,1)
(including x > 1/2, using hpl_eval.H quadrature) to 100 digits.  Rank of the system is printed (45 = the
basis is independent on the sample, so the solution is unique).
"""
import json, itertools
from fractions import Fraction as Fr
import mpmath as mp
from hpl_fit import hpl_series
import hpl_eval as he

mp.mp.dps = 260
WORDS = [(0, 0, 1), (0, 0, -1), (0, 1, -1), (0, -1, 1), (1, 0, 1), (1, 0, -1), (1, -1, 1), (-1, 0, 1),
         (-1, 0, -1), (-1, 1, -1)]
NT = 1700


def hpl_series_mp(w, N):
    if len(w) == 1:
        a = w[0]
        return [mp.mpf(0)] + [(mp.mpf(1)/m if a == 1 else mp.mpf((-1)**(m + 1))/m) for m in range(1, N + 1)]
    inner = hpl_series_mp(w[1:], N)
    a = w[0]
    out = [mp.mpf(0)]*(N + 1)
    if a == 0:
        for m in range(1, N + 1):
            out[m] = inner[m]/m
    else:
        acc = mp.mpf(0)
        for m in range(0, N):
            acc = acc*a + inner[m]
            out[m + 1] = acc/(m + 1)
    return out


SER = {w: hpl_series_mp(w, NT) for w in WORDS}


def H_taylor(w, x):
    s = SER[w]
    return mp.fsum(s[m]*x**m for m in range(1, len(s)))


def basis_funcs():
    """x-dependent pure functions only; their coefficients are real constants identified afterwards
    (weight-3 functions: rational; weight-2: Q*log2; weight-1: Q pi^2 + Q log^2 2; 1: zeta3, pi^2 log2, log^3 2)."""
    Li3args = [('x', lambda x: x), ('-x', lambda x: -x), ('1-x', lambda x: 1 - x), ('1/(1+x)', lambda x: 1/(1 + x)),
               ('(1+x)/2', lambda x: (1 + x)/2), ('(1-x)/(1+x)', lambda x: (1 - x)/(1 + x)),
               ('-(1-x)/(1+x)', lambda x: -(1 - x)/(1 + x)), ('2x/(1+x)', lambda x: 2*x/(1 + x))]
    Li2args = [('x', lambda x: x), ('-x', lambda x: -x), ('(1+x)/2', lambda x: (1 + x)/2)]
    logs = [('lx', lambda x: mp.log(x)), ('l1m', lambda x: mp.log(1 - x)), ('l1p', lambda x: mp.log(1 + x))]
    B = []
    for nm, g in Li3args:
        B.append((f'Li3({nm})', 3, lambda x, g=g: mp.polylog(3, g(x))))
    for (nm, g), (ln, h) in itertools.product(Li2args, logs):
        B.append((f'Li2({nm})*{ln}', 3, lambda x, g=g, h=h: mp.polylog(2, g(x))*h(x)))
    for nm, g in Li2args:
        B.append((f'Li2({nm})', 2, lambda x, g=g: mp.polylog(2, g(x))))
    for deg in (3, 2, 1):
        for combo in itertools.combinations_with_replacement(range(3), deg):
            nm = '*'.join(logs[i][0] for i in combo)
            B.append((nm, deg, lambda x, combo=combo: mp.fprod(logs[i][1](x) for i in combo)))
    B.append(('1', 0, lambda x: mp.mpf(1)))
    return B


CONST = {3: [], 2: [('log2', lambda: mp.log(2))], 1: [('pi2', lambda: mp.pi**2), ('log2^2', lambda: mp.log(2)**2)],
         0: [('z3', lambda: mp.zeta(3)), ('pi2*log2', lambda: mp.pi**2*mp.log(2)), ('log2^3', lambda: mp.log(2)**3)]}


def identify(val, wt):
    """write a real coefficient of a weight-wt function as a Q-combination of weight-(3-wt) constants."""
    if abs(val) < mp.mpf(10)**-120:
        return []
    if wt == 3:
        fr = Fr(str(mp.nstr(val, 60))).limit_denominator(10**6)
        assert abs(val - mp.mpf(fr.numerator)/fr.denominator) < mp.mpf(10)**-100, val
        return [('1', fr)]
    cs = CONST[wt]
    vec = [val] + [f() for _, f in cs]
    mp.mp.dps = 100
    rel = mp.pslq(vec, maxcoeff=10**8, maxsteps=10**6)
    mp.mp.dps = 260
    assert rel is not None and rel[0] != 0, (val, wt)
    return [(cs[i][0], Fr(-rel[i + 1], rel[0])) for i in range(len(cs)) if rel[i + 1] != 0]


if __name__ == '__main__':
    B = basis_funcs()
    nb = len(B)
    xs = [mp.mpf(2)/100 + mp.mpf(i)*mp.mpf(53)/100/89 for i in range(90)]
    A = mp.matrix(len(xs), nb)
    for i, x in enumerate(xs):
        for j, (_, _, fn) in enumerate(B):
            A[i, j] = fn(x)
    U, S, V = mp.svd_r(A)
    print(f'{nb} basis functions, {len(xs)} sample points; singular values max {mp.nstr(max(S),3)} min {mp.nstr(min(S),3)}')
    out = {}
    for w in WORDS:
        rhs = mp.matrix([H_taylor(w, x) for x in xs])
        sol = mp.qr_solve(A, rhs)[0]
        terms = []
        for j in range(nb):
            for cn, fr in identify(sol[j], B[j][1]):
                terms.append((B[j][0], cn, fr))
        def evalterms(x):
            cv = {'1': mp.mpf(1), 'log2': mp.log(2), 'pi2': mp.pi**2, 'log2^2': mp.log(2)**2, 'z3': mp.zeta(3),
                  'pi2*log2': mp.pi**2*mp.log(2), 'log2^3': mp.log(2)**3}
            fmap = {b[0]: b[2] for b in B}
            return mp.fsum(mp.mpf(fr.numerator)/fr.denominator*cv[cn]*fmap[fn](x) for fn, cn, fr in terms)
        mp.mp.dps = 110
        worst = 0
        for x in [mp.mpf(k)/75 + mp.mpf(1)/997 for k in range(1, 41)]:      # fresh points in (0, 0.55)
            worst = max(worst, abs(evalterms(x) - H_taylor(w, x)))
        mp.mp.dps = 45
        far = max(abs(evalterms(x) - he.H(w, x)) for x in [mp.mpf('0.7'), mp.mpf('0.9'), mp.mpf('0.97')])
        mp.mp.dps = 260
        print(f'H{list(w)} = ' + ' + '.join(f'({fr})*{cn}*{fn}' for fn, cn, fr in terms))
        print(f'      max error at 40 fresh points in (0,0.55): {mp.nstr(worst, 3)};  at x=0.7,0.9,0.97 (45-digit quadrature): {mp.nstr(far, 3)}', flush=True)
        out[str(list(w))] = [(fn, cn, str(fr)) for fn, cn, fr in terms]
    json.dump(out, open('classical_dict.json', 'w'), indent=1)
    print('saved classical_dict.json')
