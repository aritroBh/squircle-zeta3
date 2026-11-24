"""check_data.py -- referee check of TOW3 against TOWN1 exact tower data (own QSI arithmetic).
Part A: rho(5), rho(6) closed forms (parsed straight out of proof_depth6.md) vs initial data u^{(k=d)}(delta).
Part B: tower law at delta=5,6: exact finite differences in k (degree floor(delta/2) exactly).
Part C: Newton coefficients Delta^s u(k0)/rho(0) vs tow3lib e-family (off-diagonal) and ehat (diagonal).
Run in the scratch folder (tow3lib, tow2lib, pickles copied there).
"""
import pickle, re, sys, time, os
from fractions import Fraction as Fr
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PKL = '<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl'
MD = '<project>/research/08_theory/TOW3_tower_depth6/proof_depth6.md'

# ---------------- own exact arithmetic in Q(sqrt2, i): (a + b r2) + i (c + d r2) ----------------
class Q:
    __slots__ = ('a', 'b', 'c', 'd')
    def __init__(s, a=0, b=0, c=0, d=0):
        s.a, s.b, s.c, s.d = Fr(a), Fr(b), Fr(c), Fr(d)
    @staticmethod
    def of(t):
        return Q(*t)
    def __add__(s, o): return Q(s.a+o.a, s.b+o.b, s.c+o.c, s.d+o.d)
    def __sub__(s, o): return Q(s.a-o.a, s.b-o.b, s.c-o.c, s.d-o.d)
    def __neg__(s): return Q(-s.a, -s.b, -s.c, -s.d)
    def sc(s, q): return Q(s.a*q, s.b*q, s.c*q, s.d*q)
    def __mul__(s, o):
        # real part: (a+b r)(e+f r) - (c+d r)(g+h r); imag part: (a+b r)(g+h r) + (c+d r)(e+f r)
        a, b, c, d = s.a, s.b, s.c, s.d
        e, f, g, h = o.a, o.b, o.c, o.d
        re0 = a*e + 2*b*f - c*g - 2*d*h
        re1 = a*f + b*e - c*h - d*g
        im0 = a*g + 2*b*h + c*e + 2*d*f
        im1 = a*h + b*g + c*f + d*e
        return Q(re0, re1, im0, im1)
    def conj_r2(s): return Q(s.a, -s.b, s.c, -s.d)
    def inv(s):
        t = s.conj_r2()
        n = s*t                         # in Q(i)
        assert n.b == 0 and n.d == 0
        den = n.a*n.a + n.c*n.c
        assert den != 0
        g = Q(n.a/den, 0, -n.c/den, 0)  # 1/n
        return t*g
    def iszero(s): return s.a == 0 and s.b == 0 and s.c == 0 and s.d == 0
    def in_Qi(s): return s.b == 0 and s.d == 0
    def __repr__(s): return f'Q({s.a},{s.b},{s.c},{s.d})'

# sanity of own arithmetic
x = Q(1, 2, 3, 4); y = Q(-2, 1, 5, -1)
assert (x*y.inv()*y - x).iszero()
assert (x*x.inv() - Q(1)).iszero()

def gauss_of(expr):
    expr = sp.simplify(expr)
    assert expr.is_number and not expr.has(sp.Float)
    re_, im_ = sp.re(expr), sp.im(expr)
    return Fr(int(sp.numer(re_)), int(sp.denom(re_))), Fr(int(sp.numer(im_)), int(sp.denom(im_)))

def make_eval(expr, syms):
    """exact evaluator of a rational function with Gaussian-rational coefficients in syms"""
    num, den = sp.fraction(sp.together(expr))
    def prep(poly_expr):
        P = sp.Poly(sp.expand(poly_expr), *syms)
        terms = []
        for mon, co in P.terms():
            cre, cim = sp.re(co), sp.im(co)
            terms.append((mon, Fr(int(sp.numer(cre)), int(sp.denom(cre))), Fr(int(sp.numer(cim)), int(sp.denom(cim)))))
        return terms
    tn, td = prep(num), prep(den)
    def ev(terms, vals):
        R = Fr(0); I = Fr(0)
        for mon, cr, ci in terms:
            f = Fr(1)
            for v, e_ in zip(vals, mon): f *= Fr(v)**e_
            R += cr*f; I += ci*f
        return R, I
    def f(*vals):
        nr, ni = ev(tn, vals); dr, di = ev(td, vals)
        dd = dr*dr + di*di
        assert dd != 0
        # (nr + i ni)(dr - i di)/dd
        return ((nr*dr + ni*di)/dd, (ni*dr - nr*di)/dd)
    return f

d_, m_ = sp.symbols('d m')

def parse_poly(s):
    s = s.replace('\n', ' ')
    s = s.replace('^', '**')
    for _ in range(4):
        s = re.sub(r'(\d)\s+([a-zA-Z(])', r'\1*\2', s)
        s = re.sub(r'(\d)([a-zA-Z(])', r'\1*\2', s)
        s = re.sub(r'\)\s*\(', ')*(', s)
        s = re.sub(r'\)\s*([a-zA-Z])', r')*\1', s)
        s = re.sub(r'([a-zA-Z0-9])\s+([a-zA-Z(])', r'\1*\2', s)
    return sp.sympify(s, locals={'i': sp.I, 'd': d_, 'm': m_})

def load_md_closed_forms():
    txt = open(MD).read()
    m5 = re.search(r'a_5 = \((.*?)\)\s*/ \(480 \(2d-5\)\(2d-3\)\)', txt, re.S)
    m6 = re.search(r'a_6 = \((.*?)\)\s*/ \(5760 \(2d-7\)\(2d-5\)\(2d-3\)\)', txt, re.S)
    N5 = parse_poly(m5.group(1)); N6 = parse_poly(m6.group(1))
    a5 = N5/(480*(2*d_-5)*(2*d_-3))
    a6 = N6/(5760*(2*d_-7)*(2*d_-5)*(2*d_-3))
    return a5, a6

def main():
    t00 = time.time()
    H = pickle.load(open(PKL, 'rb'))
    H.pop('_meta', None)
    a5, a6 = load_md_closed_forms()
    print('parsed md closed forms; free symbols:', a5.free_symbols, a6.free_symbols, flush=True)
    import tow3lib as T3
    # md text vs library closed forms (symbolic)
    d5 = sp.cancel(sp.together(a5 - T3.arho(5))); d6 = sp.cancel(sp.together(a6 - T3.arho(6)))
    print('md a_5 - tow3lib.arho(5) == 0 :', d5 == 0, '; md a_6 - tow3lib.arho(6) == 0 :', d6 == 0, flush=True)
    ev5 = make_eval(a5, (d_, m_)); ev6 = make_eval(a6, (d_, m_))
    evs = {5: ev5, 6: ev6}

    # ---------------- Part A ----------------
    cells = {5: [0, 0, 0], 6: [0, 0, 0]}   # total, bad, d<=40 total
    bad_list = []
    stats = {}
    for dlt in (5, 6):
        for dd in range(5, 53):
            for key in sorted(k for k in H if sum(k) == dd):
                p, q = key
                if 2*dd-4 < dlt: continue
                rho0 = Q.of(H[key][0][dd]); rhod = Q.of(H[key][dlt][dd])
                assert not rho0.iszero()
                ratio = rhod*rho0.inv()
                er, ei = evs[dlt](dd, p)
                ok = ratio.in_Qi() and ratio.a == er and ratio.c == ei
                cells[dlt][0] += 1
                if dd <= 40: cells[dlt][2] += 1
                if not ok:
                    cells[dlt][1] += 1; bad_list.append((dlt, key, ratio, (er, ei)))
    for dlt in (5, 6):
        n6_40 = sum(1 for key in H if 6 <= sum(key) <= 40 and 2*sum(key)-4 >= dlt)
        print(f'PART A delta={dlt}: cells(d=5..52) = {cells[dlt][0]}, bad = {cells[dlt][1]}; cells with d in 6..40: {sum(1 for key in H if 6<=sum(key)<=40)} (bad among them counted below)', flush=True)
    bad640 = [b for b in bad_list if 6 <= sum(b[1]) <= 40]
    print('PART A bad cells total:', len(bad_list), '; bad with d in 6..40:', len(bad640), flush=True)
    for b in bad_list[:5]: print('   BAD', b)
    # count of cells in d 6..40 per delta
    for dlt in (5, 6):
        c = sum(1 for dd in range(6, 41) for key in H if sum(key) == dd and 2*dd-4 >= dlt)
        print(f'   delta={dlt}: d in 6..40 -> {c} cells checked (ratio in Q(i) asserted, exact equality of real & imag parts)', flush=True)
    # d = 5 only
    c5 = sum(1 for key in H if sum(key) == 5)
    print('   d=5 keys:', c5, '(included in the d=5..52 totals above)', flush=True)

    # ---------------- Part B ----------------
    def fdiff(vals, order):
        cur = list(vals)
        for _ in range(order):
            cur = [cur[t+1] - cur[t] for t in range(len(cur)-1)]
        return cur
    res = {}
    for dlt in (5, 6):
        P = dlt//2
        tot = 0; genuine = 0; bad = 0; nondeg = 0; nongen = 0
        tot40 = 0; gen40 = 0; bad40 = 0; nondeg40 = 0
        minpts = 10**9; minchecks = 10**9
        wrong_diag_k0 = 0; diag_cells = 0
        for key in sorted(H):
            p, q = key; dd = p+q
            if dd < 5 or 2*dd-4 < dlt: continue
            diag = (p == q)
            k0 = dd+1 if diag else dd
            ks = sorted(k for k in H[key][dlt] if k >= k0)
            assert ks == list(range(k0, 53)), (key, ks[:3], ks[-3:])
            ys = [Q.of(H[key][dlt][k]) for k in ks]
            tot += 1
            if dd <= 40: tot40 += 1
            if len(ys) < P+2:
                nongen += 1; continue
            genuine += 1
            if dd <= 40: gen40 += 1
            dP1 = fdiff(ys, P+1); dP = fdiff(ys, P)
            minpts = min(minpts, len(ys)); minchecks = min(minchecks, len(dP1))
            if not all(x.iszero() for x in dP1):
                bad += 1
                if dd <= 40: bad40 += 1
            if dP[0].iszero() or all(x.iszero() for x in dP):
                nondeg += 1
                if dd <= 40: nondeg40 += 1
            # diagonal exceptional level: is the initial level k=d on the polynomial?
            if diag:
                diag_cells += 1
                y0 = Q.of(H[key][dlt][dd])
                # extrapolate polynomial back one step: use P+1 forward diffs from k0
                # p(k0-1) = sum_{j} (-1)^j Delta^j y(k0)
                fw = []
                cur = list(ys)
                for j in range(P+1):
                    fw.append(cur[0]); cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
                ext = Q(0)
                for j, v in enumerate(fw):
                    ext = ext + (v if j % 2 == 0 else -v)
                if not (ext - y0).iszero(): wrong_diag_k0 += 1
        res[dlt] = (tot, genuine, bad, nondeg, nongen, tot40, gen40, bad40, nondeg40, minpts, minchecks, diag_cells, wrong_diag_k0)
        print(f'PART B delta={dlt} (poly degree {P}): keys d=5..52: {tot}, genuine (>= {P+2} levels): {genuine}, '
              f'non-genuine skipped: {nongen}, finite-diff order {P+1} nonzero: {bad}, order-{P} identically zero (degree < {P}): {nondeg}; '
              f'min points {minpts}, min #checks {minchecks}', flush=True)
        print(f'          restricted to d<=40: keys {tot40}, genuine {gen40}, bad {bad40}, deg-deficient {nondeg40}', flush=True)
        print(f'          diagonal keys: {diag_cells}; initial level k=2p OFF the stable polynomial in {wrong_diag_k0}', flush=True)

    # ---------------- Part C: Newton coefficients vs tow3lib e-family ----------------
    sys.setrecursionlimit(10000)
    evf = {}
    for dlt in (5, 6):
        for s in range(0, dlt//2+1):
            evf[(dlt, s)] = make_eval(T3.e(dlt, s), (T3.d, T3.m))
    badC = 0; totC = 0
    for dlt in (5, 6):
        P = dlt//2
        for key in sorted(H):
            p, q = key; dd = p+q
            if dd < 5 or 2*dd-4 < dlt or p == q: continue
            ks = sorted(k for k in H[key][dlt] if k >= dd)
            ys = [Q.of(H[key][dlt][k]) for k in ks]
            if len(ys) < P+1: continue
            rho0 = Q.of(H[key][0][dd]); r0i = rho0.inv()
            cur = list(ys)
            for s in range(P+1):
                Es = cur[0]; cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
                got = Es*r0i
                er, ei = evf[(dlt, s)](dd, p)
                totC += 1
                if not (got.in_Qi() and got.a == er and got.c == ei):
                    badC += 1
                    if badC <= 5: print('   C-BAD', dlt, key, s, got, (er, ei))
    print(f'PART C off-diagonal Newton coefficients: {totC} coefficients compared to tow3lib.e(delta,s)(d,m)*rho(0); bad = {badC}', flush=True)

    # diagonal: two ehat conventions.  PE.ehat is normalised by rho_pp(0); PA.ehat is normalised by the CHILD's
    # e00 = rho_{(p,p-1)}(0) (level k = 2p-1).  Both are compared with the data's shifted-Newton coefficients.
    import prove_adj56 as PA
    import prove_E3b56 as PE
    p = PA.p
    assert PE.sp.Symbol('p') == p   # same symbol object semantics
    badE = 0; totE = 0; badA = 0; totA = 0; badDelta = 0; totDelta = 0; badR = 0; totR = 0
    ehatE = {}; ehatA = {}
    for dlt in (5, 6):
        for s in range(dlt//2+1):
            ehatE[(dlt, s)] = sp.cancel(sp.together(PE.ehat(dlt, s)))
            ehatA[(dlt, s)] = sp.cancel(sp.together(PA.ehat(dlt, s)))
    Dl_cache = {dlt: sp.cancel(PE.Delta(dlt)) for dlt in (5, 6)}
    # symbolic link between the two conventions
    ratio_sym = sp.cancel(PA.rB0 / PE.rB0d)
    for (dlt, s) in ehatE:
        lk = sp.cancel(ehatA[(dlt, s)] - ratio_sym*ehatE[(dlt, s)])
        print(f'   symbolic link ehat_adj({dlt},{s}) - (rB0/rB0d) * ehat_E3({dlt},{s}) == 0 : {lk == 0}')
    print('   rB0_adj/rB0d_E3 =', ratio_sym)
    for dlt in (5, 6):
        P = dlt//2
        for pp in range(3, 27):
            key = (pp, pp); dd = 2*pp
            if 2*dd-4 < dlt: continue
            ks = sorted(k for k in H[key][dlt] if k >= dd+1)
            ys = [Q.of(H[key][dlt][k]) for k in ks]
            if len(ys) < P+1: continue
            rho0 = Q.of(H[key][0][dd]); r0i = rho0.inv()
            child = (pp, pp-1)
            if 2*(2*pp-1)-4 >= dlt:
                e00c = Q.of(H[child][0][2*pp-1]); e00ci = e00c.inv()
                rr = rho0*e00ci          # rho_pp(0)/e00_child  should equal rB0/rB0d(p)
                totR += 1
                exr = gauss_of(ratio_sym.subs(p, pp))
                if not (rr.in_Qi() and rr.a == exr[0] and rr.c == exr[1]):
                    badR += 1; print('   R-BAD', pp, rr, exr)
            else:
                e00ci = None
            cur = list(ys); fw = []
            for s in range(P+1):
                Es = cur[0]; fw.append(Es); cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
                got = Es*r0i
                ex = gauss_of(ehatE[(dlt, s)].subs(p, pp))
                totE += 1
                if not (got.in_Qi() and got.a == ex[0] and got.c == ex[1]):
                    badE += 1
                    if badE <= 5: print('   E-BAD', dlt, key, s, got, ex)
                if e00ci is not None:
                    gotA = Es*e00ci
                    exA = gauss_of(ehatA[(dlt, s)].subs(p, pp))
                    totA += 1
                    if not (gotA.in_Qi() and gotA.a == exA[0] and gotA.c == exA[1]):
                        badA += 1
                        if badA <= 5: print('   A-BAD', dlt, key, s, gotA, exA)
            poly_at_2p = Q(0)
            for j, v in enumerate(fw): poly_at_2p = poly_at_2p + (v if j % 2 == 0 else -v)
            Dd = (Q.of(H[key][dlt][dd]) - poly_at_2p)*r0i
            ex = gauss_of(Dl_cache[dlt].subs(p, pp))
            totDelta += 1
            if not (Dd.in_Qi() and Dd.a == ex[0] and Dd.c == ex[1]):
                badDelta += 1
                if badDelta <= 5: print('   Delta-BAD', dlt, key, Dd, ex)
    print(f'PART C diagonal ehat (E3 convention, /rho_pp(0)) vs data shifted-Newton coeffs, p=3..26: {totE} compared, bad = {badE}', flush=True)
    print(f'PART C diagonal ehat (adjacent-script convention, /e00_child) vs data: {totA} compared, bad = {badA}', flush=True)
    print(f'PART C data ratio rho_pp(0)/e00_(p,p-1) vs rB0/rB0d: {totR} compared, bad = {badR}', flush=True)
    print(f'PART C diagonal Delta_pp(delta) (E3 input, symbolic) vs data: {totDelta} compared, bad = {badDelta}', flush=True)
    print('total time', round(time.time()-t00, 1), 's')

if __name__ == '__main__':
    main()
