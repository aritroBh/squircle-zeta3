
import os
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
verify = r'''"""WT6S consolidated verification: verify_wt6s.py  (< 10 min; PASS/FAIL per item).

Reruns, in one process (shared g_symbol cache):
  A   binomial rewrite of deg5 residues == polynomials at k=6,7,10,21,52 (EXACT)
  B-i CONTROL: W6(k=5) pairings from kern_5.pkl == (449 s2/4, -30, 15 s2, -9)
  B-ii CONTROL: weight-5 E1,E2,E3 (3,2) pairings == (-17/2,0), (-2,0), (0,0)
  B   12 N_{(i,j),r} pairing vectors; E_{5,m} pairings (E_{5,4}=E_{5,5}=0 classical)
  C   W6(k) pairing polynomial; independent full-symbol check at k=7,10
  D   atom pairings == M5a table; det M = 9/512; N(E_{5,m}) solved+verified; M5b N6 control
  E   numerics: N values (sub method) vs stored; c(6)=E55/24 vs LEAD1 (>=40d); W6(6) direct quad
"""
import json, pickle, time
from fractions import Fraction as F
import sympy as sp
import flint
import qsqrt2
from w4_symbol import Q, R, A, B, ONE as ONEQ, ban_float, self_check
import wt6s_core as C6
from wt6s_core import (ZEROQ, POLE, symbol_of_N, symbol_of_kernel_terms, pairing_vector6,
                       pair_ab, B3, L2, depth2_symbol, add, scale, clean)

k = sp.Symbol('k')
KEYS5, POLES5 = ["50", "41", "32"], ["0", "a", "-a", "-b"]
KEYS4, POLES4 = ["40", "31", "22"], ["0", "a", "-a", "-b"]
RESULTS = []
T0 = time.time()

def report(name, ok, detail=""):
    RESULTS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}  {detail}  ({time.time()-T0:.0f}s)", flush=True)

def parse_poly(s):
    return sp.expand(sp.sympify(s)) if s is not None else sp.Integer(0)

def binom_coeffs(poly, base, deg):
    poly = sp.expand(poly)
    out, rem = [], poly
    for m in range(deg, -1, -1):
        bp = sp.expand(sp.prod([k - base - j for j in range(m)]) / sp.factorial(m))
        lc = sp.LC(sp.Poly(bp, k))
        cd = sp.Poly(rem, k).degree() if rem != 0 else -1
        if cd < m:
            out.append(sp.Integer(0)); continue
        c = sp.Rational(sp.Poly(rem, k).nth(m)) / lc
        out.append(c); rem = sp.expand(rem - c * bp)
    assert rem == 0
    out.reverse(); return out

def kern_pkl_to_omega(K, keys):
    P2 = {(0, 0): '0', (1, 1): 'a', (-1, -1): '-a', (-1, 1): 'b', (1, -1): '-b'}
    out = {}
    for key in keys:
        rec = []
        for (p0, p1), (c0, c1, c2, c3) in K[key]['real'].items():
            assert c2 == 0 and c3 == 0
            if c0 or c1:
                rec.append([P2[(int(p0), int(p1))], [str(c0), str(c1)]])
        out[f"{key[0]},{key[1]}"] = rec
    return out

def main():
    ban_float(); self_check()
    dd = json.load(open('deg5_residue_polynomials.json'))
    # ---------- A ----------
    okA = True
    for key in KEYS5:
        for pole in POLES5:
            cell = dd[key].get(pole, {})
            for part in ['rat', 'sqrt2']:
                poly = parse_poly(cell.get(part))
                co = binom_coeffs(poly, 6, 4)
                for kv in [6, 7, 10, 21, 52]:
                    val = sum(sp.Rational(co[m]) * sp.binomial(kv - 6, m) for m in range(5))
                    okA &= sp.expand(val - poly.subs(k, kv)) == 0
    report('A: binomial rewrite of deg5 residues exact at k=6,7,10,21,52', okA)
    e_bin = {f"{a}|{b}": {'rat': [str(c) for c in binom_coeffs(parse_poly(dd[a].get(b, {}).get('rat')), 6, 4)],
                          'sqrt2': [str(c) for c in binom_coeffs(parse_poly(dd[a].get(b, {}).get('sqrt2')), 6, 4)]}
             for a in KEYS5 for b in POLES5}

    # ---------- B-i ----------
    K5 = pickle.load(open('kern_5.pkl', 'rb'))
    om5 = kern_pkl_to_omega(K5, [(5, 0), (4, 1), (3, 2)])
    S5 = {}
    for ij in [(5, 0), (4, 1), (3, 2)]:
        S5 = add(S5, symbol_of_kernel_terms(om5[f"{ij[0]},{ij[1]}"], *ij))
    pv5 = pairing_vector6(S5)
    exp5 = [Q(0, F(449, 4)), Q(-30), Q(0, 15), Q(-9)]
    report('B-i: CONTROL W6(k=5) pairings == (449 s2/4, -30, 15 s2, -9)', pv5 == exp5,
           str([str(v) for v in pv5]))

    # ---------- B-ii ----------
    d4 = json.load(open('deg4_residue_polynomials.json'))
    pw5 = {}
    for which, mi in [('E1', 1), ('E2', 2), ('E3', 3)]:
        Sw = {}
        for key in KEYS4:
            i, j = int(key[0]), int(key[1])
            for pole in POLES4:
                cell = d4[key].get(pole, {})
                rat = binom_coeffs(parse_poly(cell.get('rat')), 5, 3)[mi]
                sq = binom_coeffs(parse_poly(cell.get('sqrt2')), 5, 3)[mi]
                c = Q(F(sp.Rational(rat)), F(sp.Rational(sq)))
                if c:
                    Sw = add(Sw, scale(symbol_of_N(i, j, pole), c))
        pw5[which] = [pair_ab(Sw, B3[0], L2), pair_ab(Sw, B3[1], L2)]
    okii = (pw5['E1'] == [Q(F(-17, 2)), Q(0)] and pw5['E2'] == [Q(-2), Q(0)]
            and pw5['E3'] == [Q(0), Q(0)])
    report('B-ii: CONTROL weight-5 E1,E2,E3 (3,2) pairings', okii,
           str({w: [str(v) for v in pw5[w]] for w in pw5}))

    # ---------- B ----------
    n_pv = {}
    for key in KEYS5:
        i, j = int(key[0]), int(key[1])
        for pole in POLES5:
            n_pv[f"{key}|{pole}"] = pairing_vector6(symbol_of_N(i, j, pole))
    ePv = {}
    ePv_pq = {}
    for m in range(1, 6):
        tot = [ZEROQ] * 4
        for key in KEYS5:
            for pole in POLES5:
                cell = e_bin[f"{key}|{pole}"]
                c = Q(F(cell['rat'][m - 1]), F(cell['sqrt2'][m - 1]))
                if c:
                    for t in range(4):
                        tot[t] += c * n_pv[f"{key}|{pole}"][t]
        ePv[f"E5,{m}"] = tot
        ePv_pq[f"E5,{m}"] = [[str(v.p), str(v.q)] for v in tot]
    expE = {'E5,1': [Q(0, 46), Q(F(-75, 2)), Q(0, 15), Q(0)],
            'E5,2': [Q(0, F(-21, 2)), Q(F(-35, 2)), Q(0), Q(21)],
            'E5,3': [Q(0, -6), Q(-10), Q(0), Q(12)],
            'E5,4': [Q(0)] * 4, 'E5,5': [Q(0)] * 4}
    okB = all(ePv[m] == expE[m] for m in expE)
    report('B: 12 N pairing vectors; E_{5,m} pairings (E_{5,4},E_{5,5} classical)', okB,
           str({m: [str(v) for v in ePv[m]] for m in ePv}))

    # ---------- C ----------
    def w6_pairings_at(kv):
        om = {}
        for key in KEYS5:
            i, j = int(key[0]), int(key[1])
            terms = []
            for pole in POLES5:
                cell = dd[key].get(pole, {})
                rat = parse_poly(cell.get('rat')).subs(k, kv)
                sq = parse_poly(cell.get('sqrt2')).subs(k, kv)
                if rat or sq:
                    terms.append([pole, [str(sp.Rational(rat)), str(sp.Rational(sq))]])
            om[f"{i},{j}"] = terms
        S = {}
        for ij in [(5, 0), (4, 1), (3, 2)]:
            S = add(S, symbol_of_kernel_terms(om[f"{ij[0]},{ij[1]}"], *ij))
        return pairing_vector6(S)
    okC = True
    for kv in [7, 10]:
        direct = w6_pairings_at(kv)
        poly = [sum((ePv[f"E5,{m}"][t] * int(sp.binomial(kv - 6, m - 1)) for m in range(1, 6)), ZEROQ)
                for t in range(4)]
        okC &= all(direct[t] == poly[t] for t in range(4))
    report('C: W6(k) pairing polynomial == independent full symbol at k=7,10', okC)

    # ---------- D ----------
    ATOMS = [("Li_{3,3}(b,1)", 3, 3, B, ONEQ), ("Li_{3,3}(1/sqrt2,1)", 3, 3, R / 2, ONEQ),
             ("Li_{3,3}(b,2)", 3, 3, B, Q(2)), ("Li_{4,2}(b,1/sqrt2)", 4, 2, B, R / 2)]
    TAB = {"Li_{3,3}(b,1)": ['0', '0', '0', '3/4'], "Li_{3,3}(1/sqrt2,1)": ['0', '3/16', '0', '0'],
           "Li_{3,3}(b,2)": ['1', '0', '0', '-3'], "Li_{4,2}(b,1/sqrt2)": ['1/4', '-1/8', '1/8', '3/4']}
    avec = {}
    okD = True
    for name, m, n, x, y in ATOMS:
        pv = pairing_vector6(depth2_symbol(m, n, x, y))
        avec[name] = pv
        okD &= [str(v) for v in pv] == TAB[name]
    M = flint.fmpq_mat([[str(avec[nm][t].p) for nm, *_ in ATOMS] for t in range(4)])
    okD &= M.det() == F(9, 512)
    Mi = M.inv()
    combo = [sum((Q(-12, 209) * avec['Li_{3,3}(b,1)'][t],), ZEROQ) for t in range(1)]  # placeholder
    # M5b N6 control
    N6 = {"Li_{3,3}(b,1)": Q(-12, 209), "Li_{3,3}(1/sqrt2,1)": Q(-160, 80),
          "Li_{3,3}(b,2)": Q(0, F(329, 4)), "Li_{4,2}(b,1/sqrt2)": Q(0, 120)}
    c5 = [sum((N6[nm] * avec[nm][t] for nm, *_ in ATOMS), ZEROQ) for t in range(4)]
    okD &= c5 == exp5
    sols = {}
    for m in range(1, 6):
        tgt = ePv[f"E5,{m}"]
        cr = Mi * flint.fmpq_mat([[str(v.p)] for v in tgt])
        cq = Mi * flint.fmpq_mat([[str(v.q)] for v in tgt])
        coeffs = [Q(F(str(cr[i, 0])), F(str(cq[i, 0]))) for i in range(4)]
        S = {}
        for (name, mm, nn, x, y), c in zip(ATOMS, coeffs):
            if c:
                S = add(S, scale(depth2_symbol(mm, nn, x, y), c))
        chk = pairing_vector6(S) if S else [ZEROQ] * 4
        okm = all(chk[t] == tgt[t] for t in range(4))
        okD &= okm
        sols[f"E5,{m}"] = {nm: str(c) for (nm, *_), c in zip(ATOMS, coeffs) if c}
    report('D: atom pv==M5a; det M=9/512; M5b N6 control; N(E_{5,m}) solved+verified', okD,
           str(sols))

    # ---------- E ----------
    import mpmath as mp
    mp.mp.dps = 120
    r2, aN, bN = mp.sqrt(2), 1 + mp.sqrt(2), mp.sqrt(2) - 1
    POLEN = {"0": mp.mpf(0), "a": aN, "-a": -aN, "-b": -bN}
    P = 12
    def N_sub(i, j, rname):
        r = POLEN[rname]
        def f1(u):
            if u == 0: return mp.mpf(0)
            vb = (1 - bN) * u**P; v = bN + vb
            l1 = mp.log(vb * (v + aN) / (4 * v)); l2 = mp.log((aN - v) * (v + bN) / (4 * v))
            return l1**i * l2**j / (v - r) * ((1 - bN) * P * u**(P - 1))
        def f2(u):
            if u == 0: return mp.mpf(0)
            da = (aN - 1) * u**P; v = aN - da
            l1 = mp.log((v - bN) * (v + aN) / (4 * v)); l2 = mp.log(da * (v + bN) / (4 * v))
            if r == aN: return -l1**i * l2**j * (P / u)
            return l1**i * l2**j / (v - r) * ((aN - 1) * P * u**(P - 1))
        return mp.quad(f1, [0, 1], maxdegree=12) + mp.quad(f2, [0, 1], maxdegree=12)
    stored = json.load(open('n6_values.json'))['values']
    okE = True
    Nv = {}
    for (i, j) in [(5, 0), (4, 1), (3, 2)]:
        for rname in POLES5:
            v = N_sub(i, j, rname)
            Nv[f"{i},{j}|{rname}"] = v
            okE &= abs(v - mp.mpf(stored[f"{i},{j}|{rname}"])) < mp.mpf(10) ** -100
    Evals = {}
    for m in range(1, 6):
        tot = mp.mpf(0)
        for key in KEYS5:
            i, j = int(key[0]), int(key[1])
            for pole in POLES5:
                cell = e_bin[f"{key}|{pole}"]
                cr, cq = F(cell['rat'][m - 1]), F(cell['sqrt2'][m - 1])
                if cr or cq:
                    c = mp.mpf(cr.numerator)/cr.denominator + r2*mp.mpf(cq.numerator)/cq.denominator
                    tot += c * Nv[f"{i},{j}|{pole}"]
        Evals[m] = tot
    c6ref = mp.mpf(json.load(open('c6_lead_reference.json'))['value'])
    dc6 = abs(Evals[5] / 24 - c6ref)
    okE &= dc6 < mp.mpf(10) ** -40
    # W6(6) direct quadrature from kern_6.pkl full kernel
    K6 = pickle.load(open('kern_6.pkl', 'rb'))
    om6 = kern_pkl_to_omega(K6, [(5, 0), (4, 1), (3, 2)])
    w6d = mp.mpf(0)
    for ij in [(5, 0), (4, 1), (3, 2)]:
        i, j = ij
        ts = [(POLEN[nm], mp.mpf(F(pq[0]).numerator)/F(pq[0]).denominator
               + r2*mp.mpf(F(pq[1]).numerator)/F(pq[1]).denominator) for nm, pq in om6[f"{i},{j}"]]
        def f(v, ts=ts, i=i, j=j):
            l1 = mp.log((v - bN) * (v + aN) / (4 * v)); l2 = mp.log((aN - v) * (v + bN) / (4 * v))
            return l1**i * l2**j * mp.fsum(c / (v - r) for r, c in ts)
        w6d += mp.quad(f, [bN, mp.mpf(1)/2, 1, 2, aN], maxdegree=12)
    okE &= abs(w6d - Evals[1]) < mp.mpf(10) ** -100
    report('E: N values (sub) vs stored (100d); c(6)=E55/24 vs LEAD1 (>=40d); W6(6)==E_{5,1} direct',
           okE, f"c6 diff {mp.nstr(dc6,3)}, W6(6) diff {mp.nstr(abs(w6d-Evals[1]),3)}")

    # ---------- summary ----------
    print("\n===== SUMMARY =====")
    nfail = sum(1 for _, ok, _ in RESULTS if not ok)
    for name, ok, _ in RESULTS:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    print(f"verify_wt6s: {'ALL PASS' if nfail == 0 else f'{nfail} FAILURES'}  ({time.time()-T0:.0f}s)")
    return nfail == 0

if __name__ == '__main__':
    import sys
    sys.exit(0 if main() else 1)
'''
open(os.path.join(WD, 'verify_wt6s.py'), 'w').write(verify)
print("written verify_wt6s.py")
import subprocess, time
PY = '<project>/.venv/bin/python'
t0 = time.time()
r = subprocess.run(['nice', '-n', '10', PY, '-B', 'verify_wt6s.py'], capture_output=True, text=True, cwd=WD, timeout=900)
print(f"runtime: {time.time()-t0:.0f}s rc={r.returncode}")
print(r.stdout)
print("STDERR:", r.stderr[-2000:])
