import itertools, math
import sympy as sp
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
        Fp = sp.Poly(sp.expand(F), g)
        for f, e in sp.factor_list(Fp.as_expr(), g)[1]:
            fd = sp.Poly(f, g) if f.has(g) else None
            if fd is None or fd.degree() != 2: continue
            N = fd.nth(0)
            if N not in (1, -1): continue
            sols.append((m, n, k, s, t, tp, sp.expand(f)))
print("quadratic-factor candidates with N = +-1:", len(sols))
seen = set()
out = []
for m, n, k, s, t, tp, f in sols:
    if (f, s, t, tp) in seen: continue
    seen.add((f, s, t, tp))
    out.append(((m,n,k,s,t,tp), f))
    print(f"  (m,n,k,s,t,t') = {(m,n,k,s,t,tp)}: g solves {f}")
