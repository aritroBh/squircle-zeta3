# Depth-2 alternating (level-2) Euler sums: exact reduction by proved relations.
#   Z(a,s;b,t) = sum_{m>l>=1} s^m t^l / (m^a l^b),  s,t in {+1,-1};  L(a,s) = sum_{m>=1} s^m/m^a.
# Relations used (all classical, each valid whenever every series involved converges):
#  (St) stuffle  L(a,s)L(b,t) = Z(a,s;b,t)+Z(b,t;a,s)+L(a+b,st)          [exact for truncated sums, then limit]
#  (Sh) shuffle  L(a,s)L(b,t) = sum_{j=1}^a C(w-j-1,b-1) Z(w-j,t;j,st) + sum_{j=1}^b C(w-j-1,a-1) Z(w-j,s;j,st)
#       [partial fractions of 1/(m^a n^b) in m+n; with weights x^{m+n}, 0<x<1, everything converges absolutely,
#        then x->1- by Abel's theorem]
#  (ShR) the same shuffle when a=1, s=+1 (or b=1, t=+1) and the other factor converges: the divergence is
#       linear; with weights x^{m+n}, -log(1-x) plays the role of T on both sides, and for a linearly divergent
#       series Abel regularisation (T=-log(1-x)) equals harmonic regularisation (T=H_N).  Hence the relation holds
#       with L(1,+) := T and Z(1,+;b,t) := T L(b,t) - Z(b,t;1,+) - L(b+1,t) (truncated stuffle), identically in T.
#  Euler: zeta(2)^2 = 5/2 zeta(4).
#  (Di) distribution  sum_{s,t} Z(a,s;b,t) = 2^{2-w} Z(a,+;b,+)   (a>=2)   [exact for truncated sums]
# Depth-1:  L(s,+) = zeta(s) (s>=2), L(s,-) = -(1-2^{1-s}) zeta(s) (s>=2), L(1,-) = -log 2.
import sympy as sp, itertools
l2 = sp.Symbol('log2'); ZS = {s: sp.Symbol('zeta%d' % s) for s in range(2, 8)}
def L(a, s):
    if a == 1: return -l2 if s == -1 else None
    return ZS[a] if s == 1 else -(1-sp.Rational(2)**(1-a))*ZS[a]
def zsym(a, s, b, t): return sp.Symbol('Z_%d%s_%d%s' % (a, 'p' if s == 1 else 'm', b, 'p' if t == 1 else 'm'))
def valid(a, s): return not (a == 1 and s == 1)
T = sp.Symbol('T')
def Lr(a, s): return T if (a == 1 and s == 1) else L(a, s)
def reduce_weight(w):
    unk = [(a, s, w-a, t) for a in range(1, w) for s in (1, -1) for t in (1, -1) if valid(a, s)]
    U = {u: zsym(*u) for u in unk}
    rels = []
    for a in range(1, w):
        b = w-a
        for s, t in itertools.product((1, -1), repeat=2):
            if valid(a, s) and valid(b, t):
                rels.append(L(a, s)*L(b, t) - U[(a, s, b, t)] - U[(b, t, a, s)] - L(w, s*t))
                sh = sum(sp.binomial(w-j-1, b-1)*U[(w-j, t, j, s*t)] for j in range(1, a+1)) \
                   + sum(sp.binomial(w-j-1, a-1)*U[(w-j, s, j, s*t)] for j in range(1, b+1))
                rels.append(L(a, s)*L(b, t) - sh)
            # regularised shuffle (ShR): exactly one factor is (1,+), the other convergent
            if (a == 1 and s == 1) != (b == 1 and t == 1):
                def Ur(A, S, B, Tt):
                    if valid(A, S): return U[(A, S, B, Tt)]
                    return T*L(B, Tt) - U[(B, Tt, 1, 1)] - L(B+1, Tt)
                sh = sum(sp.binomial(w-j-1, b-1)*Ur(w-j, t, j, s*t) for j in range(1, a+1)) \
                   + sum(sp.binomial(w-j-1, a-1)*Ur(w-j, s, j, s*t) for j in range(1, b+1))
                rr = sp.expand(Lr(a, s)*Lr(b, t) - sh)
                assert sp.expand(rr.coeff(T, 1)) == 0, rr   # divergent parts agree
                rels.append(rr.coeff(T, 0))
        if a >= 2:
            rels.append(sum(U[(a, s, b, t)] for s in (1, -1) for t in (1, -1)) - sp.Rational(2)**(2-w)*U[(a, 1, b, 1)])
    syms = [U[u] for u in unk]
    rels = [sp.expand(sp.expand(r).subs(ZS[2]**2, sp.Rational(5, 2)*ZS[4])) for r in rels]
    sol = sp.linsolve(rels, syms)
    (sol,) = list(sol)
    free = set().union(*[e.free_symbols for e in sol]) & set(syms)
    return {U[u]: sp.expand(v) for u, v in zip(unk, sol)}, free
if __name__ == '__main__':
    for w in range(2, 6):
        red, free = reduce_weight(w)
        print('weight', w, ': unknowns', len(red), ' free (irreducible modulo the relations):', free)
        for k_, v in red.items(): print('   ', k_, '=', v)
