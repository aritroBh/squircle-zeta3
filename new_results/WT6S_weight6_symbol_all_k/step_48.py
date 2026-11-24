
import os
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
step = r'''"""WT6S step D: depth-two parts N(E_{5,m}) as Q(sqrt2)-combinations of Li33/Li42.

Atoms (M5b's N6 argument set): Li33(b,1), Li33(1/sqrt2,1), Li33(b,2), Li42(b,1/sqrt2).
Compute each atom's pairing vector with wt6s_core (independent of M5a table; asserted equal),
verify the M5b N6 combination reproduces (449s2/4, -30, 15s2, -9), then solve M c = pv(E_{5,m})
exactly (flint fmpq), report det M (uniqueness), and re-verify by direct pairing recomputation.
"""
import json
from fractions import Fraction as F
import flint
import qsqrt2
from w4_symbol import Q, R, B, ONE as ONEQ, ban_float, self_check
import wt6s_core as C6
from wt6s_core import depth2_symbol, pairing_vector6, ZEROQ, log

# atoms: (name, m, n, x, y)
ATOMS = [
    ("Li_{3,3}(b,1)", 3, 3, B, ONEQ),
    ("Li_{3,3}(1/sqrt2,1)", 3, 3, R / 2, ONEQ),
    ("Li_{3,3}(b,2)", 3, 3, B, Q(2)),
    ("Li_{4,2}(b,1/sqrt2)", 4, 2, B, R / 2),
]
M5A_TABLE = {
    "Li_{3,3}(b,1)": ['0', '0', '0', '3/4'],
    "Li_{3,3}(1/sqrt2,1)": ['0', '3/16', '0', '0'],
    "Li_{3,3}(b,2)": ['1', '0', '0', '-3'],
    "Li_{4,2}(b,1/sqrt2)": ['1/4', '-1/8', '1/8', '3/4'],
}
M5B_N6 = {   # M5b n6_solution.json coefficients
    "Li_{3,3}(b,1)": Q(-12, 209),
    "Li_{3,3}(1/sqrt2,1)": Q(-160, 80),
    "Li_{3,3}(b,2)": Q(0, F(329, 4)),
    "Li_{4,2}(b,1/sqrt2)": Q(0, 120),
}
TARGET_K5 = [Q(0, F(449, 4)), Q(-30), Q(0, 15), Q(-9)]
LABELS = ['33', '42a', '42b', '42c']

def main():
    ban_float(); self_check()
    # atom pairing vectors (own computation)
    avec = {}
    for name, m, n, x, y in ATOMS:
        pv = pairing_vector6(depth2_symbol(m, n, x, y))
        avec[name] = pv
        mine = [str(v) for v in pv]
        assert mine == M5A_TABLE[name], (name, mine, M5A_TABLE[name])
        log(f"atom {name}: pv={mine}  (matches M5a table)")
    # all rational?
    assert all(v.q == 0 for pv in avec.values() for v in pv)
    # matrix M (pairing components x atoms), rational
    M = flint.fmpq_mat([[str(avec[name][t].p) for name, *_ in ATOMS] for t in range(4)])
    det = M.det()
    log(f"pairing matrix det = {det}  (nonzero => unique solution for any target)")
    assert det != 0
    Mi = M.inv()
    # CONTROL: M5b N6 combination reproduces the k=5 target
    combo = [sum((M5B_N6[name] * avec[name][t] for name, *_ in ATOMS), ZEROQ) for t in range(4)]
    log(f"CONTROL M5b N6 pairings {[str(v) for v in combo]} expect {[str(v) for v in TARGET_K5]}")
    assert combo == TARGET_K5
    log("CONTROL M5b N6: PASS")
    # solve for each E_{5,m}
    ePv = json.load(open('e5m_pairings.json'))['pairings_pq']
    out = {'det': str(det), 'atoms': [a[0] for a in ATOMS], 'solutions': {}}
    for m in range(1, 6):
        tgt = [Q(F(p), F(q)) for p, q in ePv[f"E5,{m}"]]
        # split rational/sqrt2 parts; solve M cr = tr, M cq = tq
        tr = flint.fmpq_mat([[str(v.p)] for v in tgt])
        tq = flint.fmpq_mat([[str(v.q)] for v in tgt])
        cr = Mi * tr
        cq = Mi * tq
        coeffs = [Q(F(str(cr[i, 0])), F(str(cq[i, 0]))) for i in range(4)]
        # verify by direct pairing recomputation
        S = {}
        from wt6s_core import add, scale
        for (name, mm, nn, x, y), c in zip(ATOMS, coeffs):
            if c == ZEROQ:
                continue
            S = add(S, scale(depth2_symbol(mm, nn, x, y), c))
        chk = pairing_vector6(S) if S else [ZEROQ] * 4
        ok = all(chk[t] == tgt[t] for t in range(4))
        log(f"N(E_{{5,{m}}}): target={[str(v) for v in tgt]}")
        for (name, *_), c in zip(ATOMS, coeffs):
            if c != ZEROQ:
                log(f"     ({c}) * {name}")
        log(f"     direct-pairing verification: {'PASS' if ok else 'FAIL'}")
        assert ok
        out['solutions'][f"E5,{m}"] = {
            'target': [str(v) for v in tgt],
            'coeffs': {name: [str(c.p), str(c.q)] for (name, *_), c in zip(ATOMS, coeffs) if c != ZEROQ}}
    json.dump(out, open('n_e5m_solution.json', 'w'), indent=1)
    with open('progress.log', 'a') as f:
        f.write('step_D: atom pv match M5a table; M5b N6 control PASS; det(M)=%s; N(E5,m) solved+verified m=1..5\n' % det)
    log("saved n_e5m_solution.json")

if __name__ == '__main__':
    main()
'''
open(os.path.join(WD, 'step_06.py'), 'w').write(step)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice', '-n', '10', PY, '-B', 'step_06.py'], capture_output=True, text=True, cwd=WD)
print(r.stdout)
print("STDERR:", r.stderr[-2000:])
