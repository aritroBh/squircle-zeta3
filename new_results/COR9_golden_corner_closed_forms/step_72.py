import itertools, math, sys
import sympy as sp
sys.path.insert(0, '<project>/research/08_theory/COR9_golden_corner_closed_forms')
g = sp.symbols('g')
from cor9_rank import true_rank, rank_int
import mpmath as mp
mp.mp.dps = 80
dist = {}
checked = 0
B = 6
for m, n, k in itertools.product(range(-B, B+1), repeat=3):
    if m == 0: continue
    if math.gcd(math.gcd(abs(m), abs(n)), abs(k)) != 1: continue
    if (n, k) == (0, 0): continue
    exps = [n+k, m, k, m+n, n, m+k]
    sh = -min(exps)
    for s, t, tp in itertools.product([1, -1], repeat=3):
        F = (2*t*tp*g**(n+k+sh) + 2*s*g**(m+sh) - tp*g**(k+sh) - s*t*g**(m+n+sh)
             - t*g**(n+sh) - s*tp*g**(m+k+sh))
        for f, e in sp.factor_list(sp.expand(F), g)[1]:
            if not f.has(g): continue
            fd = sp.Poly(f, g)
            if fd.degree() != 2: continue
            Tc = -fd.nth(1); N = fd.nth(0)
            if N not in (1, -1): continue
            Dd = int(Tc*Tc - 4*N)
            if Dd <= 0: continue
            sq = sp.factorint(Dd)
            Dsq = 1
            for p, ee in sq.items():
                if ee % 2: Dsq *= p
            if Dsq == 1: continue
            rD = mp.sqrt(Dsq)
            for sgn in (1, -1):
                gr = (mp.mpf(int(Tc)) + sgn*rD*mp.sqrt(mp.mpf(Dd)/Dsq**0)/mp.mpf(Dsq)**0 if False else (mp.mpf(int(Tc)) + sgn*mp.sqrt(Dd))/2)
                cp = (1 - t*gr**n)/(s*gr**m - t*gr**n)
                cv = s*gr**m*cp
                # identify (a + b sqrt Dsq)/e with small height
                hit = None
                for e2 in range(1, 500):
                    am = cv*e2
                    a_ = int(mp.nint(am - int(mp.nint(am)) + int(mp.nint(am))))
                    for b_ in range(-1000, 1001):
                        if abs((mp.mpf(a_) + b_*rD)/e2 - cv) < mp.mpf('1e-60'):
                            if math.gcd(math.gcd(abs(a_), abs(b_)), e2) == 1:
                                hit = (a_, b_, e2); break
                    if hit: break
                if not hit: continue
                r, M, pids, rel = true_rank(Dsq, *hit)
                oddcols = [j for j, pid in enumerate(pids) if pid[0] != 2]
                rv = rank_int(M); rvodd = rank_int([[row[j] for j in oddcols] for row in M])
                key = (r, rv, rvodd, r - rv)
                dist[key] = dist.get(key, 0) + 1
                checked += 1
print("reconstructed and rank-checked:", checked)
print("(rank, V, V_odd, delta) distribution:")
for kk in sorted(dist, key=lambda x: (x[0] is None, x)):
    print("  ", kk, "->", dist[kk])
