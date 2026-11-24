import cmath
# Quick numeric scan: maps determined by images (a,b,c) of (0,1,inf), all in mu_N.
# f(z) = [c(b-a) z + a(c-b)] / [(b-a) z + (c-b)]
def f13(a,b,c):
    num = c*(b-a) + 3*a*(c-b); den = (b-a) + 3*(c-b)
    if abs(den) < 1e-12: return None if abs(num)>1e-9 else 'degenerate'
    return num/den
def fm1(a,b,c):
    num = -c*(b-a) + a*(c-b); den = -(b-a) + (c-b)
    if abs(den) < 1e-12: return None if abs(num)>1e-9 else 'degenerate'
    return num/den

def classify(z, N, tol=1e-9):
    if z is None: return 'inf'
    if z == 'degenerate': return 'deg'
    if abs(z) < tol: return 'zero'
    if abs(abs(z)-1) > tol: return None
    w = z**N
    if abs(w-1) < 1e-6: return 'root'
    return None

for N in [6, 12, 18, 24, 30, 36]:
    roots = [cmath.exp(2j*cmath.pi*k/N) for k in range(N)]
    cnt_strict = 0; cnt_relaxed = 0
    examples = []
    for i in range(N):
        for j in range(N):
            if j==i: continue
            for k in range(N):
                if k==i or k==j: continue
                a,b,c = roots[i], roots[j], roots[k]
                u = f13(a,b,c); v = fm1(a,b,c)
                cu = classify(u,N); cv = classify(v,N)
                if cu=='root' and cv=='root':
                    cnt_strict += 1
                    if len(examples)<3: examples.append((i,j,k,u,v))
                if cu in ('root','zero','inf') and cv in ('root','zero','inf') and not (cu=='deg' or cv=='deg'):
                    cnt_relaxed += 1
    print(f"N={N}: strict(into mu_N)={cnt_strict}, relaxed(mu_N+0+inf)={cnt_relaxed}, examples={examples[:1]}")
