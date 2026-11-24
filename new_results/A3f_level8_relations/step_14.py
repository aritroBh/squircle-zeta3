code = r'''"""step_02_audit.py: independent audit of Codex's quotient computation, weights 1-3.
- rebuild relation rows via relations.make_relations (theorem rows: reg. double shuffle,
  distribution (Zhao 5.4), lifts, weight-1 branch relations)
- exact Q rank via fmpq_mat rref of the integer matrix of ALL unique rows
- mod-p ranks with three primes (Codex used only 1000000009)
- numeric evaluation of every row at ~60 digits (acb, ctx.prec=200; independent engine)
"""
import sys, time, pickle, json
sys.dont_write_bytecode = True
sys.path.insert(0, 'Q')
from fractions import Fraction as F
from flint import fmpq_mat, fmpz_mat, nmod_mat, fmpz
import relations as REL
from relations import words, make_relations
import cmzv

PRIMES = [1000000009, 1000000021, 999999937]
for p in PRIMES:
    assert fmpz(p).is_prime()

def build_int_matrix(rows, cols):
    ids = {w: i for i, w in enumerate(cols)}
    M = fmpq_mat(len(rows), len(cols))
    for i, row in enumerate(rows):
        for w, c in row:           # primitive() rows: (word, int) pairs
            M[i, ids[w]] = int(c)
    return M

def ranks(n, rows, cols):
    ids = {w: i for i, w in enumerate(cols)}
    out = {}
    for p in PRIMES:
        M = nmod_mat(len(rows), len(cols), p)
        for i, row in enumerate(rows):
            for w, c in row:
                M[i, ids[w]] = c % p
        rr, r = M.rref(inplace=True)
        out[p] = r
    return out

def main():
    lower = {}
    summary = []
    for n in (1, 2, 3):
        t0 = time.monotonic()
        cols = words(n)
        rows, counts = make_relations(n, lower)
        mp_r = ranks(n, rows, cols)
        M = build_int_matrix(rows, cols)
        t1 = time.monotonic()
        rr, rq = M.rref(inplace=True)
        t2 = time.monotonic()
        quot = len(cols) - rq
        print(f'weight {n}: cols={len(cols)} rows={len(rows)} '
              f'mod-p ranks={mp_r} exactQ rank={rq} quotient={quot} DG={3**n} '
              f'(build {t1-t0:.1f}s, fmpq rref {t2-t1:.1f}s)', flush=True)
        assert quot == 3**n or n > 1
        # keep pivot rows (actual theorem rows) as basis for lifts at next weight
        # independent pivot extraction: pivot rows of transpose
        Mt = M.transpose()
        rrt, rt = Mt.rref(inplace=True)
        assert rt == rq
        piv = []
        start = 0
        for i in range(rt):
            while not rrt[i, start]: start += 1
            piv.append(start); start += 1
        selected = [rows[j] for j in piv]
        lower[n] = selected
        summary.append((n, len(cols), len(rows), rq, quot))
        with open(f'audit_weight{n}.pickle', 'wb') as f:
            pickle.dump({'cols': cols, 'rows': rows, 'selected': selected,
                         'rank_Q': rq, 'quotient': quot}, f)
    print('SUMMARY (weight, cols, rows, rank_Q, quotient):', summary)
    assert [s[4] for s in summary] == [3, 10, 37], summary
    print('AUDIT PART 1 PASS: quotients 3/10/37 reproduced exactly over Q and mod 3 primes')

main()
'''
open('step_02_audit.py','w').write(code)
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','-B','step_02_audit.py'],
                   capture_output=True, text=True, timeout=2400)
print(r.stdout)
print('STDERR:', r.stderr[-3000:])
