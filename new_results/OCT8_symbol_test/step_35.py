import mpmath as mp
mp.mp.dps = 150
# difference from synthetic test
Fnum = 2*mp.polylog(3,-3) + mp.polylog(3,-5) - mp.polylog(2,-3)*mp.log(2)
# recovered combination
from fractions import Fraction as Fr
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator
G = [('Li3(-80)', Fr(-803,1128)), ('Li3(-125/3)', Fr(-3469,6768)), ('Li3(-24)', Fr(-88,141)),
     ('Li3(-15)', Fr(-3637,2256)), ('Li3(-25/2)', Fr(-328,141)), ('Li3(-9)', Fr(19213,13536)),
     ('Li3(-8)', Fr(76693,30456)), ('Li3(-5)', Fr(30743,1128)), ('Li3(-4)', Fr(18491,6768)),
     ('Li3(-3)', Fr(253,1692)), ('Li3(-16/9)', Fr(3325,27072)), ('Li3(-5/3)', Fr(3469,2256)),
     ('Li3(-3/2)', Fr(-160,47)), ('Li3(-5/4)', Fr(101,141)), ('Li3(-2/3)', Fr(-3397,1692))]
D = [('-Li2(-80)*log2', Fr(-5,47)), ('-Li2(-125/3)*log2', Fr(-2,47)), ('-Li2(-24)*log2', Fr(-7,47))]
RHS = sum(M(c)*mp.polylog(3, M(Fr(g[4:-1]))) for g,c in G)
RHS += sum(M(c)*(-mp.polylog(2, M(Fr(d[5:d.index(')')]) )))*mp.log(int(d[-1])) for d,c in D)
diff = Fnum - RHS
print("diff =", mp.nstr(diff, 30))
z3, z2, pi = mp.zeta(3), mp.zeta(2), mp.pi
basis = [z3, z2*mp.log(2), z2*mp.log(3), z2*mp.log(5), pi**3, mp.log(2)**3, mp.log(3)**3,
         z2*mp.log(2)**2, mp.log(2)*mp.log(3), mp.log(2)**2]
rel = mp.pslq([diff]+basis, tol=mp.mpf('1e-120'), maxcoeff=10**30, maxsteps=100000)
print("PSLQ vs {z3, z2log p, pi^3, ...}:", rel)
# simpler: is diff / z3 rational-ish?
print("diff/z3 =", mp.nstr(diff/z3, 30))
print("diff/pi^3 =", mp.nstr(diff/pi**3, 30))
