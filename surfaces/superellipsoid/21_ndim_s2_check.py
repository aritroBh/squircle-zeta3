# 21_ndim_s2_check.py -- check of the universal second coefficient in dimension n=4,5:
#   S_n(p) = 2^n n + 2^n n(n-1) M0/p + 2^n n(n-1) [M1 - (n-2)(pi^2/12 + Q0/6)]/p^2 + O(p^-3).
import mpmath as mp, importlib, sys
from multiprocessing import Pool
A_n = None
def work(args):
    n, p = args
    mp.mp.dps = 45
    mod = importlib.import_module("20_ndim_check")
    return (n, p), mp.nstr(2**n*mod.A_n(n, 1/mp.mpf(p)), 45)
if __name__ == "__main__":
    jobs = [(n, 10**4*2**j) for n in (4, 5) for j in range(5)]
    with Pool(len(jobs)) as pool:
        res = dict(pool.map(work, jobs))
    mp.mp.dps = 45
    from squircle_consts import consts
    M0, M1, M2 = consts()
    Q0 = mp.mpf(open("Q0_value.txt").read().strip())
    for n in (4, 5):
        ps = [10**4*2**j for j in range(5)]
        es = [1/mp.mpf(p) for p in ps]; Ss = [mp.mpf(res[(n, p)]) for p in ps]
        M = mp.matrix([[e**k for k in range(5)] for e in es]); c = mp.lu_solve(M, mp.matrix(Ss))
        s1 = 2**n*n*(n-1)*M0; s2 = 2**n*n*(n-1)*(M1 - (n-2)*(mp.pi**2/12 + Q0/6))
        print("n=%d: c0=%s (2^n n=%d)  c1-s1=%s  c2=%s  s2(formula)=%s  diff=%s" % (n, mp.nstr(c[0], 20), 2**n*n,
              mp.nstr(c[1]-s1, 3), mp.nstr(c[2], 15), mp.nstr(s2, 15), mp.nstr(c[2]-s2, 3)))
