"""WT6S step B: weight-6 symbols of the 12 N_{(i,j),r} and (3,3)/(4,2) pairings of E_{5,m}.

CONTROLS (must pass first):
 (i)  W6 at k=5 from kern_5.pkl -> (33: 449 sqrt2/4, 42: -30, 15 sqrt2, -9)  [M5a]
 (ii) weight five: E1,E2,E3 (3,2) pairings -> (-17/2,0), (-2,0), (0,0)        [WT5]
Then the 12 weight-6 N-symbols and pairing vectors, and pairings of E_{5,m}, m=1..5.
"""
import json, pickle
from fractions import Fraction as F
import sympy as sp
import qsqrt2
from w4_symbol import Q, R, A, B, ban_float, self_check
import wt6s_core as C6
from wt6s_core import (ZEROQ, ONEQ, POLE, symbol_of_N, symbol_of_kernel_terms,
                       pairing_vector6, pair_ab, B3, L2, add, scale, clean, log)

k = sp.Symbol('k')
KEYS5 = ["50", "41", "32"]
POLES5 = ["0", "a", "-a", "-b"]
KEYS4 = ["40", "31", "22"]
POLES4 = ["0", "a", "-a", "-b"]

def parse_poly(s):
    return sp.expand(sp.sympify(s)) if s is not None else sp.Integer(0)

def binom_poly(m, base):
    return sp.expand(sp.prod([k - base - j for j in range(m)]) / sp.factorial(m))

def binomial_coeffs(poly, base, deg):
    poly = sp.expand(poly)
    assert sp.Poly(poly, k).degree() <= deg
    out, rem = [], poly
    for m in range(deg, -1, -1):
        binom = binom_poly(m, base)
        lc = sp.LC(sp.Poly(binom, k))
        cdeg = sp.Poly(rem, k).degree() if rem != 0 else -1
        if cdeg < m:
            out.append(sp.Integer(0)); continue
        coeff = sp.Rational(sp.Poly(rem, k).nth(m)) / lc
        out.append(coeff)
        rem = sp.expand(rem - coeff * binom)
    assert rem == 0
    out.reverse()
    return out

def kern_pkl_to_omega(K, keys):
    PAIR_TO_POLE = {(0, 0): '0', (1, 1): 'a', (-1, -1): '-a', (-1, 1): 'b', (1, -1): '-b'}
    out = {}
    for key in keys:
        rec = []
        for (p0, p1), (c0, c1, c2, c3) in K[key]['real'].items():
            assert c2 == 0 and c3 == 0
            if c0 == 0 and c1 == 0:
                continue
            rec.append([PAIR_TO_POLE[(int(p0), int(p1))], [str(c0), str(c1)]])
        out[f"{key[0]},{key[1]}"] = rec
    return out

def pv_str(pv):
    return [str(v) for v in pv]

def main():
    ban_float(); self_check()
    # ---------------- CONTROL (i): k=5 from kern_5.pkl
    K5 = pickle.load(open('kern_5.pkl', 'rb'))
    om5 = kern_pkl_to_omega(K5, [(5, 0), (4, 1), (3, 2)])
    S6k5 = {}
    for ij in [(5, 0), (4, 1), (3, 2)]:
        S6k5 = add(S6k5, symbol_of_kernel_terms(om5[f"{ij[0]},{ij[1]}"], *ij))
    pv5 = pairing_vector6(S6k5)
    exp5 = [Q(0, F(449, 4)), Q(-30), Q(0, 15), Q(-9)]
    log(f"CONTROL (i) W6(k=5) pairings {pv_str(pv5)} expect {pv_str(exp5)}")
    assert pv5 == exp5, (pv5, exp5)
    log("CONTROL (i): PASS")

    # ---------------- CONTROL (ii): weight five E1,E2,E3 (3,2) pairings
    d4 = json.load(open('deg4_residue_polynomials.json'))
    deltas = json.load(open('deltas_deg4.json'))   # WT5 finite differences D1,D2,D3
    # extract e1,e2,e3 (base 5) from the deg4 polynomials; verify vs deltas_deg4.json
    pairings_w5 = {}
    for which, mi in [('E1', 1), ('E2', 2), ('E3', 3)]:
        Sw = {}
        for key in KEYS4:
            i, j = int(key[0]), int(key[1])
            for pole in POLES4:
                cell = d4[key].get(pole, {})
                rat = binomial_coeffs(parse_poly(cell.get('rat')), 5, 3)[mi]
                sq = binomial_coeffs(parse_poly(cell.get('sqrt2')), 5, 3)[mi]
                c = Q(sp.Rational(rat).p, sp.Rational(rat).q) if False else None
                c = Q(F(sp.Rational(rat)), F(sp.Rational(sq)))
                if c == ZEROQ:
                    continue
                Sw = add(Sw, scale(symbol_of_N(i, j, pole), c))
        pv = [pair_ab(Sw, B3[0], L2), pair_ab(Sw, B3[1], L2)]
        pairings_w5[which] = pv_str(pv)
        log(f"CONTROL (ii) {which} (3,2) pairings {pv_str(pv)}")
    # verify the deg4-polynomial e's reproduce WT5 deltas_deg4.json (exact)
    okd = True
    for key in KEYS4:
        i, j = int(key[0]), int(key[1])
        dkey = f"{i},{j}"
        for mi, dn in [(1, 'D1'), (2, 'D2'), (3, 'D3')]:
            ref = {t[0]: (F(t[1][0]), F(t[1][1])) for t in deltas[dkey][dn]}
            got = {}
            for pole in POLES4:
                cell = d4[key].get(pole, {})
                rat = binomial_coeffs(parse_poly(cell.get('rat')), 5, 3)[mi]
                sq = binomial_coeffs(parse_poly(cell.get('sqrt2')), 5, 3)[mi]
                if rat or sq:
                    got[pole] = (F(sp.Rational(rat)), F(sp.Rational(sq)))
            okd &= (got == ref)
    log(f"CONTROL (ii) deg4-polynomial binomial e's == WT5 deltas_deg4.json: {okd}")
    assert okd
    assert pairings_w5['E1'] == ['-17/2', '0'], pairings_w5['E1']
    assert pairings_w5['E2'] == ['-2', '0'], pairings_w5['E2']
    assert pairings_w5['E3'] == ['0', '0'], pairings_w5['E3']
    log("CONTROL (ii): PASS  (E1=(-17/2,0), E2=(-2,0), E3=(0,0); deltas match)")

    # ---------------- the 12 weight-6 N-symbols and their pairing vectors
    n_pv = {}
    n_terms = {}
    for key in KEYS5:
        i, j = int(key[0]), int(key[1])
        for pole in POLES5:
            S = symbol_of_N(i, j, pole)
            pv = pairing_vector6(S)
            n_terms[f"{key}|{pole}"] = len(S)
            n_pv[f"{key}|{pole}"] = pv
            log(f"N_{{({i},{j}),{pole}}}: symbol words={len(S)} pv={pv_str(pv)}")
    json.dump({'symbol_words': n_terms,
               'pairing_vectors': {kk: pv_str(vv) for kk, vv in n_pv.items()},
               'pairing_vectors_pq': {kk: [[str(v.p), str(v.q)] for v in vv] for kk, vv in n_pv.items()}},
              open('n6_symbols_pairings.json', 'w'), indent=1)

    # ---------------- pairings of E_{5,m}, m=1..5
    e = json.load(open('e5m_binomial_coeffs.json'))
    ePv = {}
    ePv_pq = {}
    for m in range(1, 6):
        tot = [ZEROQ] * 4
        for key in KEYS5:
            for pole in POLES5:
                cell = e[f"{key}|{pole}"]
                c = Q(F(cell['rat'][m - 1]), F(cell['sqrt2'][m - 1]))
                if c == ZEROQ:
                    continue
                pvn = n_pv[f"{key}|{pole}"]
                for t in range(4):
                    tot[t] = tot[t] + c * pvn[t]
        ePv[f"E5,{m}"] = pv_str(tot)
        ePv_pq[f"E5,{m}"] = [[str(v.p), str(v.q)] for v in tot]
        log(f"E_{{5,{m}}} pairings (33,42a,42b,42c): {pv_str(tot)}  classical={all(v==ZEROQ for v in tot)}")
    json.dump({'pairings': ePv, 'pairings_pq': ePv_pq}, open('e5m_pairings.json', 'w'), indent=1)
    with open('progress.log', 'a') as f:
        f.write('step_B: controls (i) k=5 449s2/4,(-30,15s2,-9) PASS; (ii) w5 E1,E2,E3 PASS; '
                '12 N-symbols; E5,m pairings saved\n')
    log("saved n6_symbols_pairings.json, e5m_pairings.json")

if __name__ == '__main__':
    main()
