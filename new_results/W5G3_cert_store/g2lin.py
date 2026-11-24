"""bases() and exact_combination(): copied verbatim from W5G1_frame_galois/common.py (read-only source)."""
from g2common import *
def bases():
    """base_n rows (list of dict word->Fraction) used by W5C1 lift families, n=1..4."""
    out = {}
    for n in (1, 2, 3):
        sel = pickle.load(open(T + f'/A3c_claude_elementary/Q/weight{n}.pickle', 'rb'))[2]
        out[n] = [{tuple(w): F(c) for w, c in row} for row in sel]
    T0 = pickle.load(open(T + '/A3h_weight4_target/w4_systems_T.pkl', 'rb'))['T0']
    S = pickle.load(open(T + '/A3h5_exact_certificate/stage2_rows_S.pkl', 'rb'))['S']
    b4 = []
    for i in S:
        d = {}
        for w, (q, c) in T0[i][1]:
            assert q == 0
            d[tuple(w)] = F(c)
        b4.append(d)
    out[4] = b4
    return out
FAM_N = {'lift1': 1, 'lift2old': 2, 'lift3old': 3, 'lift4old': 4}

def exact_combination(gens, rhss, cols, p=1000000007, tag=''):
    """Exact: for each rhs (dict word->Fraction) find x (list of Fraction, len(gens)) with sum x_i gens_i == rhs,
    or None if rhs not in the Q-span. gens: list of dicts with integer (or Fraction) coefficients.
    Method: independent rows + pivot columns mod p, square fmpz solve, exact Fraction verification on all cols."""
    from flint import nmod_mat, fmpz_mat
    from math import lcm
    ids = {w: i for i, w in enumerate(cols)}
    def modv(c): c = F(c); return c.numerator * pow(c.denominator, -1, p) % p
    M = nmod_mat(len(gens), len(cols), p)
    for i, g in enumerate(gens):
        for w, c in g.items(): M[i, ids[w]] = modv(c)
    tr, rk = M.transpose().rref()
    I = []; j = 0
    for i in range(rk):
        while tr[i, j] == 0: j += 1
        I.append(j); j += 1
    MI = nmod_mat(len(I), len(cols), p)
    for a, i in enumerate(I):
        for w, c in gens[i].items(): MI[a, ids[w]] = modv(c)
    rr, rk2 = MI.rref(); assert rk2 == rk
    C = []; j = 0
    for i in range(rk):
        while rr[i, j] == 0: j += 1
        C.append(j); j += 1
    Cpos = {c: k for k, c in enumerate(C)}
    # integer scaling of gens
    gl = [lcm(*[F(c).denominator for c in gens[i].values()]) if gens[i] else 1 for i in I]
    A = fmpz_mat(rk, rk)     # A[k, a] = gens[I[a]][C[k]] * gl[a]  (columns = generators)
    for a, i in enumerate(I):
        for w, c in gens[i].items():
            k = Cpos.get(ids[w])
            if k is not None: A[k, a] = int(F(c) * gl[a])
    rl = [lcm(*[F(c).denominator for c in r.values()]) if r else 1 for r in rhss]
    Bm = fmpz_mat(rk, len(rhss))
    for t, r in enumerate(rhss):
        for w, c in r.items():
            k = Cpos.get(ids[w])
            if k is not None: Bm[k, t] = int(F(c) * rl[t])
    X = A.solve(Bm)
    out = []
    for t, r in enumerate(rhss):
        x = [F(0)] * len(gens)
        for a, i in enumerate(I):
            v = X[a, t]
            x[i] = F(int(v.numerator), int(v.denominator)) * gl[a] / rl[t]
        acc = {}
        for i, xi in enumerate(x):
            if xi:
                for w, c in gens[i].items(): acc[w] = acc.get(w, F(0)) + xi * c
        ok = all(acc.get(w, F(0)) == F(r.get(w, 0)) for w in set(acc) | set(r))
        out.append(x if ok else None)
    return out, rk
