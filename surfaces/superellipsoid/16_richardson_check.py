# 16_richardson_check.py -- independent check of s1, s2, s3 from S(p) itself at huge p.
# S(p) is evaluated from the exact 1-D formula (superlib.py) at p = 10^8 * 2^j (j=0..7), dps 90,
# and the polynomial  S = c0 + c1 e + ... + c7 e^7  (e=1/p) is fitted exactly through the 8 points.
# Then c0 must be 24, c1 = 48 M0, c2 = s2 (closed form, Theorem B), c3 = s3 (numerical).
import mpmath as mp, json, sys
from multiprocessing import Pool
DPS = 90
PS = [10**8 * 2**j for j in range(8)]
def work(p):
    mp.mp.dps = DPS
    from superlib import S_area
    return p, mp.nstr(S_area(p), DPS)
if __name__ == "__main__":
    with Pool(len(PS)) as pool:
        res = dict(pool.map(work, PS))
    json.dump({str(k): v for k, v in res.items()}, open("S_values_huge.json", "w"), indent=1)
    mp.mp.dps = DPS
    es = [mp.mpf(1)/p for p in PS]; Ss = [mp.mpf(res[p]) for p in PS]
    n = len(PS)
    M = mp.matrix([[e**k for k in range(n)] for e in es]); c = mp.lu_solve(M, mp.matrix(Ss))
    from squircle_consts import consts
    M0, M1, M2 = consts()
    r2 = mp.sqrt(2); lam, l3, L = mp.log(2), mp.log(3), mp.log(1+r2); pi = mp.pi; b = r2-1; al = mp.atan(1/r2)
    Li2 = lambda x: mp.polylog(2, x)
    s2 = (-5*pi**2 + 37*r2/2*pi**2 - 12*r2*pi*al + 48*lam**2 + 96*r2*lam*L + 72*r2*l3*L
          - 96*L**2 - 168*r2*L**2 - 192*r2*Li2(b) - 18*r2*Li2(-b**4))
    try:
        Q1 = mp.mpf(json.load(open("Q1_hp_250.json"))["Q1"])
    except Exception:
        Q1 = mp.mpf(open("Q1_value_50.txt").read().strip())
    s3 = 48*M2 - 8*pi**2*M0 + 48*mp.zeta(3) - 8*Q1
    for k, ref, name in ((0, mp.mpf(24), "24"), (1, 48*M0, "s1=48 M0"), (2, s2, "s2 (closed form)"), (3, s3, "s3 = 48M2-8pi^2M0+48zeta3-8Q1")):
        print("c%d = %s\n     %s = %s   residual %s" % (k, mp.nstr(c[k], 45), name, mp.nstr(ref, 45), mp.nstr(c[k]-ref, 3)))
    print("c4 =", mp.nstr(c[4], 25), "  (s4 estimate)")
