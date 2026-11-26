"""s07: the 'top' parts [zeta(2)^e] P_{2e,n} (wt(z) = k), symbolic n, for e = 1,2,3 (mu_side of g8_allN, exact).
By the ladder (Lemma ladderprecise) every part [z]M_k with wt(z) = k is a rational combination of E_e = [zeta(2)^e]M_{2e}."""
import sys, time
sys.setrecursionlimit(20000)
import sympy as sp
import g8_allN as G
for e in (1, 2, 3):
    t0 = time.time()
    P = G.mu_side(2*e)
    zm = tuple([e] + [0]*(2*e - 2))
    E = P[zm]
    print(f"[z2^{e}]P_({2*e},n) = ", {w: sp.factor(c.as_expr()) for w, c in E.items()}, f"  [{time.time()-t0:.0f}s]")
    if list(E) == [()]:
        n = sp.Symbol('n')
        print("   partial fractions:", sp.apart(E[()].as_expr(), n), flush=True)
