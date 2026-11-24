import itertools, math, sys
import sympy as sp
sys.path.insert(0, '<project>/research/08_theory/COR9_golden_corner_closed_forms')
g = sp.symbols('g')
sols = []
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
            disc = Tc*Tc - 4*N
            if disc <= 0 or sp.sqrt(disc, evaluate=False) == int(sp.sqrt(disc)): 
                if disc <= 0: continue
            if disc <= 0: continue
            sols.append((m, n, k, s, t, tp, Tc, N))
print("real-quadratic candidates:", len(sols))
# reconstruct c and check rank
from cor9_rank import true_rank
import mpmath as mp
mp.mp.dps = 60
found = []
for m, n, k, s, t, tp, Tc, N in sols:
    Dd = Tc*Tc - 4*N
    for sgn in (1, -1):
        gr = (mp.mpf(int(Tc)) + sgn*mp.sqrt(int(Dd)))/2
        if abs(gr) < 1e-50: continue
        cp = (1 - t*gr**n)/(s*gr**m - t*gr**n)
        cv = s*gr**m*cp
        if not (mp.almosteq(cv, cv) and abs(cv) < 1e6): continue
        # identify field of c: c should be in Q(sqrt Dd') with Dd' squarefree part of Dd
        sq = sp.factorint(int(Dd))
        Dsq = 1
        for p, ee in sq.items():
            if ee % 2: Dsq *= p
        rD = mp.sqrt(Dsq)
        # write cv = (a + b sqrt Dsq)/e via nsimplify-ish: try small e
        hit = None
        for e in range(1, 200):
            a_ = mp.nint(e*(cv + 0*rD)); 
            for b_ in range(-400, 401):
                if abs((mp.mpf(a_) + b_*rD)/e - cv) < mp.mpf('1e-40'):
                    a2 = int(a_); 
                    if math.gcd(math.gcd(abs(a2), abs(b_)), e) == 1:
                        hit = (a2, b_, e); break
            if hit: break
        if not hit: continue
        a2, b_, e = hit
        r, M, pids, rel = true_rank(Dsq, a2, b_, e)
        oddcols = [j for j, pid in enumerate(pids) if pid[0] != 2]
        rv = cor9_rank_rank = None
        from cor9_rank import rank_int
        rv = rank_int(M); rvodd = rank_int([[row[j] for j in oddcols] for row in M])
        key = (r, rv, rvodd, (r or 0) - (rv or 0))
        if r is not None and r <= 2:
            found.append((Dsq, a2, b_, e, key, (m, n, k, s, t, tp), cv))
            print(f"  RANK<=2! c = ({a2},{b_})/{e} in Q(sqrt{Dsq}) ~ {mp.nstr(cv,12)}: (rank,V,Vodd,delta) = {key}, from (m,n,k,s,t,t') = {(m,n,k,s,t,tp)}")
print("total rank<=2 reconstructed:", len(found))
