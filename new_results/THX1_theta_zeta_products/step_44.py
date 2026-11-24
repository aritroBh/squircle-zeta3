
cert_module = '''"""thx1_cert.py -- certificates of MZV identities from proved relation instances.

Pipeline (per weight W):
  relations R_j (DS+DUAL instances, exact Fractions)  [proved identities]
  columns = sorted convergent words
  flint rref(R) mod p  -> pivot columns
  flint rref(R^T) mod p -> independent original row indices B
  for target V (exact Fraction word-combo): solve h_B . R_B = V over Q:
      solve mod primes, CRT, rational reconstruction, EXACT verification (sparse Fraction)
  certificate: V = sum_j h_j R_j  (exact); then value(V)=0 since each R_j is proved.
"""
from fractions import Fraction as Fr
import flint

def build_nmod(rels, cols, p):
    colpos = {c: i for i, c in enumerate(cols)}
    A = flint.nmod_mat(len(rels), len(cols), p)
    for i, r in enumerate(rels):
        for k, v in r.items():
            A[i, colpos[k]] = (v.numerator % p) * pow(v.denominator % p, -1, p) % p
    return A, colpos

def rref_pivots(A):
    """return (rank, pivot_col_indices) from flint rref."""
    R = A.rref()
    nrows, ncols = R.nrows(), R.ncols()
    piv = []
    for i in range(nrows):
        for j in range(ncols):
            if R[i, j] != 0:
                piv.append(j)
                break
    return len(piv), piv, R

def structure(rels, cols, p):
    """returns dict(rank, pivcols, indep_rows)"""
    A, colpos = build_nmod(rels, cols, p)
    rank, pivcols, _ = rref_pivots(A)
    At = A.transpose()
    rank2, indep_rows, _ = rref_pivots(At)
    assert rank == rank2
    return {'rank': rank, 'pivcols': pivcols, 'indep_rows': indep_rows, 'colpos': colpos}

def vec_to_nmod_row(V, cols, p, colpos):
    M = flint.nmod_mat(1, len(cols), p)
    for k, v in V.items():
        M[0, colpos[k]] = (v.numerator % p) * pow(v.denominator % p, -1, p) % p
    return M

def solve_h_mod(V, rels, B, pivcols, cols, colpos, p):
    """solve h (len r) with sum_i h_i R_{B[i]} = V  (mod p), using pivot columns."""
    r = len(B)
    S = flint.nmod_mat(r, r, p)
    for i, bi in enumerate(B):
        for k, v in rels[bi].items():
            j = colpos[k]
            if j in pivcols:
                S[pivcols.index(j), i] = (v.numerator % p) * pow(v.denominator % p, -1, p) % p
    b = flint.nmod_mat(r, 1, p)
    for j, pc in enumerate(pivcols):
        # V coefficient at column pc
        pass
    Vv = {colpos[k]: v for k, v in V.items()}
    for j, pc in enumerate(pivcols):
        v = Vv.get(pc, Fr(0))
        b[j, 0] = (v.numerator % p) * pow(v.denominator % p, -1, p) % p
    x = S.solve(b)
    h = [int(x[i, 0]) % p for i in range(r)]
    return h
'''
open('thx1_cert.py','w').write(cert_module)
print("written thx1_cert.py")
