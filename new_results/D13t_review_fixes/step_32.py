import sympy as sp
k = sp.symbols('k', integer=True, positive=True)
ck = (-1)**(k+1)*sp.binomial(sp.Rational(1,2), k)
S = sp.summation(ck, (k, 1, sp.oo))
print("sympy summation sum_{k>=1} (-1)^{k+1} binom(1/2,k) =", S)
# direct partial sum at K=2e6 via recursion
import mpmath as mp
mp.mp.dps = 50
c = mp.mpf('0.5'); s = c; K = 2000000
for kk in range(2, K+1):
    c *= -(mp.mpf('0.5')-kk+1)/kk   # binom(1/2,k) recursion: b_k = b_{k-1}*(0.5-k+1)/k; c_k = (-1)^{k+1} b_k
    # careful: recursion for b_k then sign
print("(recursion quick check skipped; using closed tail model)")
# partial sum trend S(K) vs 1 - 1/sqrt(pi K)
def partial(K):
    b = mp.mpf('0.5'); tot = b
    for kk in range(2, K+1):
        b *= (mp.mpf(kk)-mp.mpf('1.5'))/kk   # |binom(1/2,k)| = |b_{k-1}|*(k-1.5)/k, and c_k = |b_k|
        tot += b
    return tot
for K in [10**4, 10**5, 10**6]:
    SK = partial(K)
    print(f"  S({K}) = {mp.nstr(SK,12)} ; 1 - 1/sqrt(pi K) = {mp.nstr(1-1/mp.sqrt(mp.pi*K),12)} ; diff = {mp.nstr(SK-(1-1/mp.sqrt(mp.pi*K)),6)}")
