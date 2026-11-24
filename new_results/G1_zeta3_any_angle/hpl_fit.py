"""hpl_fit.py -- find the generating function of the exact flat-corner Taylor coefficients as a combination
of harmonic polylogarithms (HPLs) of c with letters {0, 1, -1} and prefactors c^e, e in {-1, 0, 1}.

Input: flat_exact_N{N}_K{K}.json (exact mu_{k,n} in Q[z2,z3]).
For a given k and zeta-monomial ('1', 'z2', 'z3'), the target series is  T(c) = sum_n coef_n c^(2n).
Basis: c^e H_w(c), w nonempty word over {0,1,-1} of length <= 3 with last letter != 0, plus 1 and c.
HPLs (Remiddi-Vermaseren): H_1 = -log(1-c), H_{-1} = log(1+c), H_0 = log c,
  H_{a,w}(c) = int_0^c f_a(t) H_w(t) dt,  f_0 = 1/t, f_1 = 1/(1-t), f_{-1} = 1/(1+t).
Exact rational linear algebra (python-flint fmpq_mat) on Taylor coefficients 0..ORDER (odd ones = 0).
The system is overdetermined (ORDER+1 equations); a solution is accepted only if it satisfies ALL equations.
"""
import sys, json, itertools
from fractions import Fraction as Fr
import flint

FILE = sys.argv[1] if len(sys.argv) > 1 else 'flat_exact_N100_K3.json'
K = int(sys.argv[2]) if len(sys.argv) > 2 else 2
MONO = sys.argv[3] if len(sys.argv) > 3 else '(0, 0, 0)'
MAXLEN = int(sys.argv[4]) if len(sys.argv) > 4 else 3
EXPS = [int(e) for e in (sys.argv[5].split(',') if len(sys.argv) > 5 else ['-1', '0', '1'])]

data = json.load(open(FILE))[str(K)]
NMAX = max(int(n) for n in data)
ORDER = 2*NMAX

def target():
    T = [Fr(0)]*(ORDER + 1)
    for nn, P in data.items():
        T[2*int(nn)] = Fr(P.get(MONO, '0'))
    return T

def hpl_series(w, order):
    """Taylor coefficients of H_w(c), list length order+3 (index = power of c)."""
    N = order + 3
    if len(w) == 1:
        a = w[0]
        s = [Fr(0)]*(N + 1)
        for m in range(1, N + 1):
            s[m] = Fr(1, m) if a == 1 else Fr((-1)**(m + 1), m)
        return s
    inner = hpl_series(w[1:], order)
    a = w[0]
    out = [Fr(0)]*(N + 1)
    if a == 0:
        for m in range(1, N + 1):
            out[m] = inner[m]/m
    else:
        # f_a(t) * inner(t) = sum_m t^m sum_{i<=m} inner_i a^{m-i}  (1/(1-a t) expansion)
        acc = Fr(0)
        for m in range(0, N):
            acc = acc*a + inner[m]          # coefficient of t^m in inner/(1 - a t)
            out[m + 1] = acc/(m + 1)
    return out

def words(maxlen):
    W = []
    for L in range(1, maxlen + 1):
        for w in itertools.product([0, 1, -1], repeat=L):
            if w[-1] != 0:
                W.append(w)
    return W

if __name__ == '__main__':
    T = target()
    basis = [('1', [Fr(1)] + [Fr(0)]*(ORDER + 2)), ('c', [Fr(0), Fr(1)] + [Fr(0)]*(ORDER + 1))]
    for w in words(MAXLEN):
        s = hpl_series(w, ORDER)
        for e in EXPS:
            if e >= 0:
                ser = [Fr(0)]*e + s[:len(s) - e]
            else:
                ser = s[-e:] + [Fr(0)]*(-e)
            basis.append((f'c^{e} H{list(w)}', ser))
    nb = len(basis)
    print(f'k={K} monomial {MONO}: {nb} basis functions, {ORDER+1} equations (n <= {NMAX})')
    Mx = flint.fmpq_mat(ORDER + 1, nb)
    for j, (_, ser) in enumerate(basis):
        for i in range(ORDER + 1):
            v = ser[i]
            Mx[i, j] = flint.fmpq(v.numerator, v.denominator)
    rhs = flint.fmpq_mat(ORDER + 1, 1)
    for i in range(ORDER + 1):
        rhs[i, 0] = flint.fmpq(T[i].numerator, T[i].denominator)
    # least-squares-free exact solve: augment and row reduce
    aug = flint.fmpq_mat(ORDER + 1, nb + 1)
    for i in range(ORDER + 1):
        for j in range(nb):
            aug[i, j] = Mx[i, j]
        aug[i, nb] = rhs[i, 0]
    R, rank = aug.rref()
    rankA = Mx.rref()[1]
    print('rank(A) =', rankA, ' rank([A|b]) =', rank)
    if rank != rankA:
        print('NO SOLUTION in this basis')
        sys.exit(0)
    if rankA < nb:
        print('WARNING: basis not independent on these coefficients (rank deficit', nb - rankA, ')')
    # read particular solution from rref (pivot columns)
    sol = [flint.fmpq(0)]*nb
    for i in range(rank):
        row = [R[i, j] for j in range(nb + 1)]
        piv = next(j for j in range(nb + 1) if row[j] != 0)
        if piv == nb:
            print('inconsistent'); sys.exit(0)
        sol[piv] = row[nb]
    terms = [(basis[j][0], sol[j]) for j in range(nb) if sol[j] != 0]
    print(f'SOLUTION ({len(terms)} terms):')
    for name, v in terms:
        print(f'   {v}  *  {name}')
    json.dump([(nm, str(v)) for nm, v in terms], open(f'hpl_fit_k{K}_{MONO.replace(" ", "")}.json', 'w'))
