code = r'''"""step_04_audit_w4.py: weight-4 audit.
1) numeric evaluation of ALL 17625 relation rows at ~60 digits (acb, prec=200)
2) mod-p rank with a second prime 1000000021 (Codex used only 1000000009), plus D residue check
"""
import sys, pickle, time, json
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
from fractions import Fraction as F

def numeric():
    import check_conversion as CC
    from check_conversion import evaluate, approx, mp
    with open('Q/weight4.pickle', 'rb') as f:
        cols, rows = pickle.load(f)
    t0 = time.monotonic(); worst = mp.mpf(0); bad = 0
    for k, row in enumerate(rows):
        v = abs(approx(evaluate(dict(row))))
        worst = max(worst, v)
        if v > mp.mpf('1e-40'):
            bad += 1
            if bad <= 3: print('BAD', k, v, row[:5])
        if k % 4000 == 3999:
            print(f'  checked {k+1}/{len(rows)} worst so far {mp.nstr(worst,8)} ({time.monotonic()-t0:.0f}s)', flush=True)
    print(f'weight 4 numeric: rows={len(rows)} bad={bad} worst={mp.nstr(worst,8)} ({time.monotonic()-t0:.0f}s)', flush=True)
    return bad == 0

def modular():
    from flint import nmod_mat, fmpz
    with open('Q/weight4.pickle', 'rb') as f:
        cols, rows = pickle.load(f)
    p = 1000000021
    assert fmpz(p).is_prime()
    ids = {w: i for i, w in enumerate(cols)}
    t0 = time.monotonic()
    M = nmod_mat(len(rows), len(cols), p)
    for i, row in enumerate(rows):
        for w, c in row:
            M[i, ids[w]] = c % p
    print(f'matrix built ({time.monotonic()-t0:.0f}s), rref mod {p}...', flush=True)
    rr, rank = M.rref(inplace=True)
    print(f'weight 4 mod {p}: rank={rank} quotient={len(cols)-rank} ({time.monotonic()-t0:.0f}s)', flush=True)
    # pivots
    pp = []; start = 0
    for i in range(rank):
        while not rr[i, start]: start += 1
        pp.append(start); start += 1
    free = [j for j in range(len(cols)) if j not in set(pp)]
    # D residues
    data = json.loads(open('Q/vectors.json').read())
    res = []
    for part in data['D']:
        v = [0]*len(cols)
        for w, c in part:
            c = F(c); v[ids[tuple(w)]] = c.numerator*pow(c.denominator, -1, p) % p
        rem = [(v[j] - sum(v[q]*int(rr[i, j]) for i, q in enumerate(pp))) % p for j in free]
        res.append(sum(bool(x) for x in rem))
    print(f'D residue nonzero coordinates mod {p}: {res} of {len(free)}', flush=True)
    return rank, len(cols)-rank

ok = numeric()
rank, quot = modular()
print(f'AUDIT w4: numeric_ok={ok} rank={rank} quotient={quot}')
assert ok and quot == 142
print('AUDIT w4 PASS (second prime confirms 142; all rows numerically zero at 60 digits)')
'''
open('step_04_audit_w4.py','w').write(code)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','-B','step_04_audit_w4.py'],
                   capture_output=True, text=True, timeout=3000)
print(r.stdout); print('STDERR:', r.stderr[-2000:])
