import mpmath as mp
mp.mp.dps = 80
def chi8(n):
    r = n % 8
    return {1:1, 7:1, 3:-1, 5:-1}.get(r, 0)
h = mp.mpf(8)**-5 * (mp.zeta(5, mp.mpf(1)/8) + mp.zeta(5, mp.mpf(7)/8)
                     - mp.zeta(5, mp.mpf(3)/8) - mp.zeta(5, mp.mpf(5)/8))
# direct series with Richardson/Shanks via nsum with more tolerance, and plain partial sums
for N in [10**4, 10**6]:
    s = mp.fsum(chi8(n)/mp.mpf(n)**5 for n in range(1, N+1))
    # tail ~ sum_{n>N} chi8(n)/n^5 ~ chi-averaged ~ integral with zero mean -> O(1/N^5); but cancellation
    print(N, mp.nstr(s, 50), "diff vs hurwitz:", mp.nstr(abs(s-h), 6))
# nsum with better method
s2 = mp.nsum(lambda n: chi8(n)/mp.mpf(n)**5, [1, mp.inf], tol=mp.mpf('1e-70'), method='euler-maclaurin')
print("nsum EM:", mp.nstr(s2, 60), mp.nstr(abs(s2-h), 6))
