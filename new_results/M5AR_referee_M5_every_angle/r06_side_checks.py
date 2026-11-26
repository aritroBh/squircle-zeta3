"""r06 (referee, own code): side claims.
 (a) [z2^2]mu_{5,n} and [z4]mu_{5,n} formulas printed by s03 (with Z_{-1}(2n)) vs own exact mu, n<=150.
 (b) top part [z2^3]mu_{6,n} = n(5n-1)/(3(2n-1)(2n+1)(2n+3)(2n+5)) vs own exact mu at k=6, n<=80.
 (c) uniqueness: rank (mod a prime) of the coefficient map of the fit bases {1, c^-2, c^e He[w]: -4<=e<=1, |w|<=W, last letter 1}.
"""
import os, sys, pickle, itertools
from fractions import Fraction as Fr
os.environ['KMAX'] = '6'
import r01_mu_exact as R
MU = pickle.load(open("r01_mu_k0to5_n150.pkl", "rb"))
Z = lambda **kw: tuple(kw.get(f"z{w}", 0) for w in range(2, 6))
def Zm1(M): return sum((Fr((-1)**i, i) for i in range(1, M+1)), Fr(0))
okA = okB = True
for n in range(1, 151):
    d = 4*(n+1)*(n+2)*(2*n-1)*(2*n+1)**2*(2*n+3)**2
    z22 = Fr(-n, (2*n-1)*(2*n+1)*(2*n+3))*Zm1(2*n) - Fr(84*n**5 + 296*n**4 + 275*n**3 - 44*n**2 - 137*n - 36, d)
    z4 = Fr(n+3, (2*n-1)*(2*n+1)*(2*n+3))*Zm1(2*n) - Fr(4*n**6 - 128*n**5 - 457*n**4 - 346*n**3 + 231*n**2 + 312*n + 72, n*d)
    okA &= MU[5][n].get(Z(z2=2)) == z22
    okB &= MU[5][n].get(Z(z4=1)) == z4
print("(a) [z2^2]mu5 formula (s03) == own mu, n<=150:", okA, "; [z4]mu5 formula == own mu:", okB)
m6 = tuple([3, 0, 0, 0, 0])
okC = all(R.mu(n, 6).get(m6, 0) == Fr(n*(5*n-1), 3*(2*n-1)*(2*n+1)*(2*n+3)*(2*n+5)) for n in range(1, 81))
print("(b) [z2^3]mu_(6,n) = n(5n-1)/(3(2n-1)(2n+1)(2n+3)(2n+5)), own mu at k=6, n<=80:", okC)
mons6 = sorted({R.name(m) for n in (1, 2, 3, 7) for m in R.mu(n, 6)})
print("    monomials of mu_(6,n):", mons6)
# (c) ranks
import flint
p = 2305843009213693951
def words_upto(W):
    return [tuple(pre) + (1,) for L in range(1, W+1) for pre in itertools.product((0, 1, -1), repeat=L-1)]
def coeffs_mod(words, NT):
    inv = [0] + [pow(m, -1, p) for m in range(1, NT+1)]
    memo = {}
    def get(w):
        if w in memo: return memo[w]
        if len(w) == 1:
            r = [0] + [inv[m] if (w[0] == 1 or m % 2) else (-inv[m]) % p for m in range(1, NT+1)]
        else:
            inner = get(w[1:]); a = w[0]; r = [0]*(NT+1)
            if a == 0:
                for m in range(1, NT+1): r[m] = inner[m]*inv[m] % p
            else:
                s = 0
                for m in range(1, NT+1):
                    s = (s*a + inner[m-1]) % p; r[m] = s*inv[m] % p
        memo[w] = r; return r
    return {w: get(w) for w in words}, {w: get(tuple(-x for x in w)) for w in words}
for W, label in ((1, "[z5],[z2z3]"), (2, "G, [z4]"), (3, "[z3]"), (4, "Psi3"), (6, "R5")):
    words = words_upto(W); ncol = 2 + 6*len(words)
    nrow_even = ncol + 40
    NT = 2*nrow_even + 10
    A, Am = coeffs_mod(words, NT)
    rows = []
    Ns = list(range(-6, 2*nrow_even + 1))       # all Laurent orders (even and odd)
    for N in Ns:
        row = [1 if N == 0 else 0, 1 if N == -2 else 0]
        for e in range(-4, 2):
            for w in words:
                m = N - e
                if m < 1: row.append(0); continue
                v = (A[w][m] + (1 if (e % 2 == 0) == (m % 2 == 0) else -1)*A[w][m]) % p   # c^e H_w(c) + (-c)^e H_w(-c)
                row.append(v)
        rows.append(row)
    M = flint.nmod_mat(len(rows), ncol, [x for r in rows for x in r], p)
    print(f"(c) basis W={W} ({label}): {ncol} unknowns, {len(rows)} equations (c^N, -6<=N<={Ns[-1]}), rank mod p = {M.rank()}", flush=True)
