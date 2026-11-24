import sympy as sp, pickle, os, time
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
t = sp.Symbol('t')
L0 = sp.log(t); L1 = sp.log(1-t)
Li2 = lambda z: sp.polylog(2, z)
Li3 = lambda z: sp.polylog(3, z)
def p1norm(e):
    # polylog(1, z) -> -log(1-z)
    return e.replace(lambda f: isinstance(f, sp.Function) and f.func == sp.polylog and f.args[0]==1,
                     lambda f: -sp.log(1-f.args[1]))
def check0(e):
    return sp.simplify(p1norm(sp.expand(e))) == 0
def int_t(n,a,b):
    """antiderivative of t^n L0^a L1^b (a+b<=3), certificate-verified."""
    if b==0:
        F = 0
        for i in range(a+1):
            F += (-1)**(a-i) * sp.factorial(a)/sp.factorial(i) * L0**i / (n+1)**(a-i+1)
        return sp.expand(t**(n+1)*F)
    if a==0:
        # int t^n L1^b = t^(n+1)L1^b/(n+1) + (b/(n+1))[ -L1^b/b - sum_{j=0}^n int t^j L1^(b-1) ]
        F = t**(n+1)*L1**b/(n+1) + sp.Rational(b,n+1)*(-L1**b/b)
        for j in range(n+1):
            F -= sp.Rational(b,n+1)*int_t(j,0,b-1)
        return sp.expand(F)
    # a>=1,b>=1: int t^n L0^a L1^b = B(n,a)(t)*L1^b + b*int B(n,a)(t) L1^(b-1)/(1-t) dt
    B = int_t(n,a,0)
    F = sp.expand(B*L1**b)
    # B = sum c_i t^(n+1) L0^i ; B/(1-t): sub s=1-t handled directly in t:
    # int t^(n+1) L0^i L1^(b-1)/(1-t) dt : expand t^(n+1) = (1-(1-t))^(n+1) = sum_j (-1)^j C(n+1,j) (1-t)^j
    s = 1-t
    acc = 0
    for i in range(a+1):
        ci = sp.expand(B).coeff(L0,i)
        # ci should be c_i * t^(n+1) (possibly 0)
        if ci == 0: continue
        ci_t = sp.simplify(ci/t**(n+1))
        assert not ci_t.has(t), (i, ci)
        for j in range(n+2):
            cj = (-1)**j * sp.binomial(n+1,j)
            # int (1-t)^j L0^i L1^(b-1)/(1-t) dt with sign ds=-dt and t->s=1-t:
            # int t^(n+1) L0^i L1^(b-1)/(1-t) dt = -int (1-s)^n+1... we already expanded:
            # term = cj * int (1-t)^j L0^i L1^(b-1)/(1-t) dt
            if j == 0:
                term = cj * int_over_t(i, b-1)   # int L0^i L1^(b-1)/(1-t) dt  [s-form]
            else:
                term = cj * int_polylogweight(j-1, i, b-1)  # int (1-t)^(j-1) L0^i L1^(b-1) dt
            acc += ci_t*term
    return sp.expand(F + b*acc)
def int_over_t(i, q):
    """int L0^i L1^q/(1-t) dt = -int (log(1-s))^i (log s)^q / s ds  (s=1-t) -> -J(q,i)."""
    return -int_s_over(q, i)
def int_s_over(p, q):
    """int (log s)^p (log(1-s))^q / s ds  (p+q<=2)  [s-form, s=1-t]."""
    s = 1-t
    ls = sp.log(s); l1s = sp.log(1-s)
    if (p,q)==(0,1): return -Li2(s)
    if (p,q)==(1,1): return -Li2(s)*ls + Li3(s)
    if (p,q)==(0,2): return ls*l1s**2 + 2*l1s*Li2(1-s) - 2*Li3(1-s)
    if (p,q)==(1,0): return ls**2/2
    if (p,q)==(2,0): return ls**3/3
    if (p,q)==(0,0): return ls
    raise ValueError((p,q))
def int_polylogweight(m, p, q):
    """int (1-t)^m L0^p L1^q dt, m>=0, p+q<=2 (and not the /t case)."""
    # expand (1-t)^m = sum_j (-1)^j C(m,j) t^j
    acc = 0
    for j in range(m+1):
        cj = (-1)**j*sp.binomial(m,j)
        acc += cj*int_t(j,p,q)
    return sp.expand(acc)
# build the 40 local antiderivatives + verify certificates
O5 = '<project>/research/08_theory/OCT6_octagon_closed/oct5'
Jloc = pickle.load(open(O5+'/Jloc.pkl','rb'))
keys = sorted(Jloc.keys())
print("keys:", len(keys))
t0=time.time()
bad = 0
anti = {}
for (n,a,b) in keys:
    F = int_t(n,a,b)
    f = t**n * L0**a * L1**b
    ok = check0(sp.diff(F,t) - f)
    if not ok:
        bad += 1
        print("CERT FAIL", (n,a,b))
    anti[(n,a,b)] = F
print("certificates verified:", len(keys)-bad, "of", len(keys), " in %.1fs"%(time.time()-t0))
pickle.dump(anti, open('Jloc_antiderivatives.pkl','wb'))
