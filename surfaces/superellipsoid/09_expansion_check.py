# 09_expansion_check.py -- S(p) for p = 10 .. 10^5 from the exact 1-D formula (superlib.py),
# and comparison with the asymptotic series 24 + s1/p + s2/p^2 + s3/p^3 (+ s4/p^4) using
#   s1 = 48 M0, s2 = 48 M1 - 4 pi^2 - 8 Q0, s3 = 48 M2 - 8 pi^2 M0 + 48 zeta(3) - 8 Q1.
# If a log p term were present, the scaled remainders below would drift like log p.
import mpmath as mp, sys, json
from multiprocessing import Pool
DPS = 60
PS = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 100000]

def work(p):
    mp.mp.dps = DPS
    from superlib import S_area
    return p, mp.nstr(S_area(p), DPS)

if __name__ == "__main__":
    with Pool(len(PS)) as pool:
        res = dict(pool.map(work, PS))
    json.dump(res, open("S_values.json", "w"), indent=1)
    mp.mp.dps = DPS
    r2 = mp.sqrt(2); L = mp.log(1+r2); lam = mp.log(2); pi = mp.pi
    from squircle_consts import consts
    M0, M1, M2 = consts()
    Q0 = mp.mpf(open("Q0_value.txt").read().strip())
    Q1 = mp.mpf(open("Q1_value.txt").read().strip())
    s1 = 48*M0; s2 = 48*M1 - 4*pi**2 - 8*Q0; s3 = 48*M2 - 8*pi**2*M0 + 48*mp.zeta(3) - 8*Q1
    print("s1 =", mp.nstr(s1, 40)); print("s2 =", mp.nstr(s2, 40)); print("s3 =", mp.nstr(s3, 40))
    print("%8s %42s %26s %26s %26s" % ("p", "S(p)", "(S-24)p", "(S-24-s1/p)p^2", "(S-..-s2/p^2)p^3", ))
    for p in PS:
        S = mp.mpf(res[p]); e = mp.mpf(1)/p
        r1 = (S-24)/e; r2_ = (S-24-s1*e)/e**2; r3 = (S-24-s1*e-s2*e**2)/e**3
        r4 = (S-24-s1*e-s2*e**2-s3*e**3)/e**4
        print("%8d %42s %26s %26s %26s   R4=%s" % (p, mp.nstr(S, 38), mp.nstr(r1, 20), mp.nstr(r2_, 20), mp.nstr(r3, 20), mp.nstr(r4, 15)))
