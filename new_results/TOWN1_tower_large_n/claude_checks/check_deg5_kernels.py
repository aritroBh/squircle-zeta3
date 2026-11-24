# M1 session: degree-5 kernels (weight six) from TOWN1's exact degree-6 towers (depths 0..8), k = 6..52.
# (1) compare key by key with the INDEPENDENT full-reducer kernels (WTS kern_k.pkl, k = 5..18; TWOEND1 red/kern_k.pkl, k = 19..21);
# (2) test that every residue is a quartic polynomial in k on 6 <= k <= 52 (exact finite differences), no simple poles at +-i;
# (3) print the residue polynomials. Triangular system (main.tex eq:triangular) at degree 5, k > 5:
#     W50 = -6 G60 w1 - G51 w2,  W41 = -5 G51 w1 - 2 G42 w2,  W32 = -4 G42 w1 - 3 G33 w2 - 3 F(G33 w1).
import os, sys, pickle, json
from fractions import Fraction as F
R = '<project>/research/08_theory'
sys.path.insert(0, f'{R}/TOWN1_tower_large_n'); sys.path.insert(0, f'{R}/WTS_weight_k_structure/theory')
import kernel_from_tower as kft
from fastred import QSI, ZERO, AQ, BQ, pf_add, pf_scale, pf_mul_DL, pf_fold, DL1, DL2

H = pickle.load(open(f'{R}/TOWN1_tower_large_n/data/towers_D8_K52.pkl', 'rb'))
meta = H.pop('_meta'); print('towers meta', meta, flush=True)
KMAX = meta['k_done']
def h_orders(key, k):
    d = key[0] + key[1]
    return {2*d - 2 - delta: kft.unq4(kd[k]) for delta, kd in H[key].items() if k in kd}
def kernel_deg5_level(k):
    g60, g51, g42, g33 = (kft.gm_from_tower(h_orders(key, k)) for key in [(6, 0), (5, 1), (4, 2), (3, 3)])
    W50 = pf_add(pf_scale(pf_mul_DL(g60, DL1), F(-6)), pf_scale(pf_mul_DL(g51, DL2), F(-1)))
    W41 = pf_add(pf_scale(pf_mul_DL(g51, DL1), F(-5)), pf_scale(pf_mul_DL(g42, DL2), F(-2)))
    W32 = pf_add(pf_add(pf_scale(pf_mul_DL(g42, DL1), F(-4)), pf_scale(pf_mul_DL(g33, DL2), F(-3))),
                 pf_scale(pf_fold(pf_mul_DL(g33, DL1)), F(-3)))
    return {key: kft.simple_part_of(W)[0] for key, W in [((5, 0), W50), ((4, 1), W41), ((3, 2), W32)]}
def as_dict(P):  # TOWN1 PF -> {(r0, r1): (c0, c1)}, plus flag for any i-part / simple poles at +-i
    got = {(r.p0, r.p1): (c.p0, c.p1) for r, c in P[0].items() if c.p0 or c.p1}
    imag = any(c.p2 or c.p3 for c in P[0].values()) or any(not c.is_zero() for c in P[1][:2]) or any(not c.is_zero() for c in P[2][:2])
    return got, imag
def ref_dict(P):  # reducer kernel (either format)
    if isinstance(P, dict) and 'real' in P:
        got = {p: (c[0], c[1]) for p, c in P['real'].items() if c[0] or c[1]}
        imag = any(c[2] or c[3] for c in P['real'].values()) or any(any(e) for e in P['ipos']) or any(any(e) for e in P['ineg'])
    else:
        got = {(r.p0, r.p1): (c.p0, c.p1) for r, c in P[0].items() if c.p0 or c.p1}
        imag = any(c.p2 or c.p3 for c in P[0].values())
    return got, imag

KS = list(range(6, KMAX + 1))
ker = {}
for k in KS:
    ker[k] = {key: as_dict(P) for key, P in kernel_deg5_level(k).items()}
print(f'degree-5 kernels from towers built for k = 6..{KMAX}', flush=True)

# (1) against the independent reducers
n = bad = 0
for k in range(6, 22):
    fn = f'{R}/WTS_weight_k_structure/theory/kern_{k}.pkl' if k <= 18 else f'{R}/TWOEND1_weight_poles/red/kern_{k}.pkl'
    Kref = pickle.load(open(fn, 'rb'))
    for key in [(5, 0), (4, 1), (3, 2)]:
        g, gi = ker[k][key]; r, ri = ref_dict(Kref[key])
        ok = (g == r) and not gi and not ri
        n += 1; bad += (not ok)
        if not ok: print(f'  MISMATCH k={k} key={key}: towers={g} imag={gi} | reducer={r} imag={ri}')
print(f'(1) towers vs reducers, k = 6..21, 3 keys: {n - bad}/{n} match', 'ALL MATCH' if not bad else 'MISMATCHES', flush=True)

# (2) cubic law on 5..KMAX, exact finite differences, per key/pole/(rational, sqrt2) part
def fd_poly(ys, P):
    cols = [list(ys)]
    for _ in range(P + 1):
        cols.append([cols[-1][t + 1] - cols[-1][t] for t in range(len(cols[-1]) - 1)])
    return all(c == 0 for c in cols[P + 1]), cols[P][0] != 0, len(cols[P + 1])
imag_any = any(ker[k][key][1] for k in KS for key in ker[k])
seqs = {}
for k in KS:
    for key, (g, _) in ker[k].items():
        for p in g: seqs.setdefault((key, p), None)
cells = ok_cells = 0; minchecks = 10**9; top = {}
for (key, p) in sorted(seqs):
    for part in (0, 1):
        ys = [F(ker[k][key][0].get(p, (0, 0))[part]) for k in KS]
        if not any(ys): continue
        is_poly, deg3, nch = fd_poly(ys, 4)
        cells += 1; ok_cells += is_poly; minchecks = min(minchecks, nch)
        if deg3: top[(key, p, part)] = True
        if not is_poly: print(f'  NOT QUARTIC: key={key} pole={p} part={"rat" if part == 0 else "sqrt2"}')
print(f'(2) quartic law k = 6..{KMAX}: {ok_cells}/{cells} residue sequences quartic (min {minchecks} spare checks); '
      f'i-parts/simple +-i poles present: {imag_any}', flush=True)

# (3) residue polynomials y(k) = sum c_j k^j from 4 points (exact)
def interp(ks, ys):
    import sympy as sp
    k = sp.Symbol('k'); return sp.expand(sp.interpolate(list(zip(ks, [sp.Rational(y.numerator, y.denominator) for y in ys])), k))
PNAME = {(0, 0): '0', (1, 1): 'a', (-1, -1): '-a', (-1, 1): 'b', (1, -1): '-b'}
out = {}
for (key, p) in sorted(seqs):
    for part in (0, 1):
        ys = [F(ker[k][key][0].get(p, (0, 0))[part]) for k in KS[:5]]
        if not any(F(ker[k][key][0].get(p, (0, 0))[part]) for k in KS): continue
        poly = interp(KS[:5], ys)
        out.setdefault(f'{key[0]}{key[1]}', {}).setdefault(PNAME.get((int(p[0]), int(p[1])), str(p)), {})['rat' if part == 0 else 'sqrt2'] = str(poly)
        print(f'   K_{key[0]}{key[1]} pole {PNAME.get((int(p[0]), int(p[1])), p):>3} {"rat  " if part == 0 else "sqrt2"}: {poly}')
json.dump(out, open(f'{R}/TOWN1_tower_large_n/claude_checks/deg5_residue_polynomials.json', 'w'), indent=1)
print('RESULT:', 'PASS' if (not bad and ok_cells == cells and not imag_any) else 'FAIL')
