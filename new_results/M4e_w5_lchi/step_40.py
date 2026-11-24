script = r'''"""M4e: saturated integer kernel of the S0 symbol matrix (65 atoms, 32 words).

Adapted from M4d nice_sat.py, for nice_search.S0 (no extras). HNF of the
transposed integer symbol matrix; last (atoms-rank) rows of U give a saturated
Z-basis of the symbol kernel. Prints the index of the RREF kernel lattice
(nice_symbol.json, primitive) in the saturated one.
"""
import json, sys, time
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path
import flint
import nice_search as NS

HERE = Path(__file__).resolve().parent
def log(s): print(s, flush=True)
def prim(vec):
    D = lcm(*[c.denominator for c in vec]); ints=[int(c*D) for c in vec]; g=0
    for i in ints: g=gcd(g,abs(i))
    return [i//g for i in ints]

def main():
    S = list(NS.S0)
    sol = NS.solve(S, 'satS0')
    assert sol['ok']
    names = sol['names']; basis = sol['basis']; words = NS.WORDS
    n = len(names)
    cols = []
    for name, s, _ in basis:
        assert all(c.q == 0 for c in s.values()), name
        cols.append([s.get(w, NS.ZERO).p for w in words])
    D2 = lcm(*[c.denominator for col in cols for c in col])
    M2 = [[int(c * D2) for c in col] for col in cols]   # atoms x words
    rows = [[M2[ai][wi] for ai in range(n)] for wi in range(32)]  # words x atoms
    A = flint.fmpz_mat(rows)
    assert A.nrows() == 32 and A.ncols() == n
    t0 = time.monotonic()
    H, U = A.transpose().hnf(transform=True)
    log(f'HNF of transpose: {time.monotonic()-t0:.1f}s, U is {U.nrows()}x{U.ncols()}')
    rank = sum(1 for i in range(H.nrows()) if any(H[i, j] for j in range(H.ncols())))
    log(f'rank = {rank} (expected 30)')
    Ksat = [[int(U[i, j]) for j in range(n)] for i in range(rank, n)]
    for krow in Ksat:
        for wi in range(32):
            assert sum(F(M2[ai][wi]) * F(krow[ai]) for ai in range(n)) == 0
    log('saturated kernel rows have zero symbol: PASS')
    snf = flint.fmpz_mat(Ksat).snf()
    diag = [int(snf[i, i]) for i in range(min(snf.nrows(), snf.ncols()))]
    log(f'saturated kernel SNF invariants: {sorted(set(diag))} (all 1 => saturated)')
    # index of old RREF kernel lattice in saturated: write old rows as T*Ksat
    Kold = [prim([F(c) for c in v]) for v in sol['kernel']]
    assert len(Kold) == len(Ksat)
    # solve T over Q: Kold = T Ksat  => T = Kold Ksat^T (Ksat Ksat^T)^{-1}
    import sympy as sp
    Ko = sp.Matrix(Kold); Ks = sp.Matrix(Ksat)
    G = (Ks * Ks.T).inv()
    T = Ko * Ks.T * G
    detT = sp.factor(T.det())
    log(f'det(T) = {detT}  (Kold = T*Ksat); index [Ksat:Kold] = {sp.factor(1/T.det()) if T.det()!=0 else "?"}')
    # LLL-reduce saturated kernel
    t0 = time.monotonic()
    red = flint.fmpz_mat(Ksat).lll(gram='exact')
    K = [[int(red[i, j]) for j in range(red.ncols())] for i in range(red.nrows())]
    log(f'LLL-reduce saturated kernel: {time.monotonic()-t0:.1f}s; max entries {sorted(max(abs(i) for i in r) for r in K)}')
    (HERE / 'nice_ksat_s0.json').write_text(json.dumps(dict(names=names, K=K, index=str(1/T.det()))) + '\n')
    log('saved nice_ksat_s0.json')

if __name__ == '__main__':
    main()
'''
M4E = '<project>/research/08_theory/M4e_w5_lchi'
open(M4E+'/nice_sat_s0.py','w').write(script)
import subprocess, time
PY = '<project>/.venv/bin/python'
t0=time.monotonic()
r = subprocess.run(['nice', PY, 'nice_sat_s0.py'], cwd=M4E, capture_output=True, text=True, timeout=1500)
print(r.stdout); print("STDERR:", r.stderr[-2000:]); print(f"elapsed {time.monotonic()-t0:.0f}s")
open(M4E+'/progress.log','a').write("step_40: nice_sat_s0.py (saturated S0 kernel + index)\n")
