# M1 session: settle LEAD2R's d = 9, 10 spare-point count and do the end-to-end c(w) check it left unfinished (w = 2..11).
# Own code, exact Fractions for the kernels: for each d, key (i,j), pole, Q(sqrt2)-component, take the residue sequence on
# k = d+1..21 (WTS kern_k.pkl for k <= 18, TWOEND1 red/kern_k.pkl for k = 19..21), test that the d-th forward difference vanishes
# (polynomial of degree <= d-1), count the checks, extract the top coefficient [k^{d-1}] = (d-1)-st difference/(d-1)!, then
# integrate the top coefficients and compare with the law c(d+1) = sigma_{d-1} A_d/(2 (d-1)!) [y^{d+1}] E(y).
import pickle, sys
from fractions import Fraction as F
from math import factorial
import mpmath as mp
mp.mp.dps = 50
R = '<project>/research/08_theory'
sys.path.insert(0, f'{R}/WTS_weight_k_structure/theory')
import fastred  # noqa: F401  (needed to unpickle QSI objects in the TWOEND1 files)
def ref_dict(P):
    if isinstance(P, dict) and 'real' in P:
        got = {(F(p[0]), F(p[1])): (F(c[0]), F(c[1])) for p, c in P['real'].items() if c[0] or c[1]}
        imag = any(c[2] or c[3] for c in P['real'].values()) or any(any(e) for e in P['ipos']) or any(any(e) for e in P['ineg'])
    else:
        got = {(F(r.p0), F(r.p1)): (F(c.p0), F(c.p1)) for r, c in P[0].items() if c.p0 or c.p1}
        imag = any(c.p2 or c.p3 for c in P[0].values())
    return got, imag
KER = {}
for k in range(2, 22):
    fn = f'{R}/WTS_weight_k_structure/theory/kern_{k}.pkl' if k <= 18 else f'{R}/TWOEND1_weight_poles/red/kern_{k}.pkl'
    KER[k] = pickle.load(open(fn, 'rb'))
def diffs(ys, n):
    cols = [list(ys)]
    for _ in range(n):
        cols.append([cols[-1][t + 1] - cols[-1][t] for t in range(len(cols[-1]) - 1)])
    return cols
s2 = mp.sqrt(2); A = 1 + s2; B = s2 - 1
tt = lambda v: (v*v + 2*v - 1)/(4*v)
L1 = lambda v: mp.log(tt(v)); L2 = lambda v: mp.log(1 - tt(v))
E = lambda y: mp.gamma(1 + y/2)*mp.gamma(mp.mpf(1)/2 - y/2)/mp.sqrt(mp.pi)
Ecoef = mp.taylor(E, 0, 12)
def zigzag(n):  # A_n via the Seidel-Entringer triangle
    row = [1]
    for m in range(1, n + 1):
        new = [0]
        for x in reversed(row): new.append(new[-1] + x)
        row = new
    return row[-1]
def law(w):
    m = w - 2; sig = (-1)**((m + 1)//2)
    return sig*zigzag(w - 1)/(2*mp.factorial(w - 2))*Ecoef[w]
allok = True
for d in range(1, 11):
    ks = list(range(d + 1, 22))
    minchk = 10**9; ok = True; top = {}
    for key in [(d - j, j) for j in range(d//2 + 1)]:
        seq = {k: ref_dict(KER[k][key]) for k in ks}
        if any(im for _, im in seq.values()): ok = False; print(f'  d={d} key={key}: i-part present')
        poles = set().union(*[set(g) for g, _ in seq.values()])
        for p in poles:
            for part in (0, 1):
                ys = [seq[k][0].get(p, (F(0), F(0)))[part] for k in ks]
                cols = diffs(ys, d)
                if not all(c == 0 for c in cols[d]): ok = False; print(f'  d={d} key={key} pole={p} part={part}: NOT degree <= {d-1}')
                minchk = min(minchk, len(cols[d]))
                tc = cols[d - 1][0]/factorial(d - 1)
                if tc: top[(key, p, part)] = tc
    # integrate the top coefficients: sum over keys of int top_kernel(v) L1^i L2^j dv
    def integrand(v):
        tot = 0
        for (key, p, part), c in top.items():
            pole = mp.mpf(p[0].numerator)/p[0].denominator + mp.mpf(p[1].numerator)/p[1].denominator*s2
            coef = mp.mpf(c.numerator)/c.denominator*(s2 if part else 1)
            tot += coef/(v - pole)*L1(v)**key[0]*L2(v)**key[1]
        return tot
    cw = mp.re(mp.quad(integrand, [B, mp.mpf('0.6'), 1, mp.mpf('1.8'), A]))
    sq2_top = any(part for (_, _, part) in top)
    res = cw - law(d + 1)
    good = ok and abs(res) < mp.mpf(10)**-30
    allok &= good
    print(f'd={d} (w={d+1}): points k={d+1}..21 ({len(ks)}), spare checks = {minchk}, degree<=d-1: {ok}, '
          f'sqrt2 in top: {sq2_top}, c(w) = {mp.nstr(cw, 25)}, law = {mp.nstr(law(d+1), 25)}, resid {mp.nstr(res, 3)}', flush=True)
print('RESULT:', 'PASS' if allok else 'FAIL')
