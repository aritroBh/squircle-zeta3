# SM6 item 3 (weight 7): c7 numeric vs pure basis; free moment gamma33 vs 21-element alternating weight-7 basis.
import mpmath as mp, itertools, sys
from multiprocessing import Pool
DPS = 420
def integ(spec):
    a, b, c, den = spec
    mp.mp.dps = DPS
    # int_0^1 (-log x)^a log^b(1-x) log^c(1+x) / den dx ; den in {'x','1+x','1-x'... } encoded
    if den == 'x': w = lambda x: 1/x
    elif den == '1+x': w = lambda x: 1/(1+x)
    else: raise
    f = lambda x: (-mp.log(x))**a*mp.log(1-x)**b*mp.log(1+x)**c*w(x)
    return mp.quad(f, [0, mp.mpf(1)/2, 1], maxdegree=8)
CAND = [(4, 1, 1, 'x'), (3, 2, 1, 'x'), (3, 1, 2, 'x'), (2, 2, 2, 'x'), (5, 1, 0, '1+x'), (4, 2, 0, '1+x'), (4, 0, 2, 'x')]
def gam(spec):
    a, b = spec
    mp.mp.dps = DPS
    f = lambda x: (-mp.log(x))**a*mp.log(1+x)**b/x
    return mp.quad(f, [0, mp.mpf(1)/4, mp.mpf(1)/2, 1], maxdegree=8)*mp.mpf(2)**(-a-1)/(mp.factorial(a)*mp.factorial(b))
if __name__ == "__main__":
    with Pool(8) as P:
        cand_vals = P.map(integ, CAND + [(4, 1, 0, '1+x')])
        gams = P.map(gam, [(3, 3), (1, 5)])
    import pickle
    pickle.dump((CAND, cand_vals, gams), open("sm6_w7_raw.pkl", "wb"))
    print("done")
