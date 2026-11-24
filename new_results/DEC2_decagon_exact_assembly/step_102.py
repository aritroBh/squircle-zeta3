import sympy as sp, pickle
t = sp.Symbol('t'); s_ = sp.Symbol('s')
L0 = sp.log(t); L1s = sp.log(1-t)
Li2 = lambda z: sp.polylog(2,z); Li3 = lambda z: sp.polylog(3,z)
def int_s_over(p,q):
    s__ = 1-t; ls = sp.log(s__); l1s_ = sp.log(1-s__)
    if (p,q)==(0,1): return -Li2(s__)
    if (p,q)==(1,1): return -Li2(s__)*ls + Li3(s__)
    if (p,q)==(0,2): return ls*l1s_**2 + 2*l1s_*Li2(1-s__) - 2*Li3(1-s__)
    if (p,q)==(1,0): return ls**2/2
    if (p,q)==(2,0): return ls**3/3
    if (p,q)==(0,0): return ls
    raise ValueError((p,q))
def int_polylogweight(m,p,q):
    acc = 0
    for j in range(m+1):
        acc += (-1)**j*sp.binomial(m,j)*int_t(j,p,q)
    return sp.expand(acc)
def int_t(n,a,b):
    if b==0:
        Fv = 0
        for i in range(a+1):
            Fv += (-1)**(a-i) * sp.factorial(a)/sp.factorial(i) * L0**i / (n+1)**(a-i+1)
        return sp.expand(t**(n+1)*Fv)
    if a==0:
        Fv = t**(n+1)*L1s**b/(n+1) + sp.Rational(b,n+1)*(-L1s**b/b)
        for j in range(n+1):
            Fv -= sp.Rational(b,n+1)*int_t(j,0,b-1)
        return sp.expand(Fv)
    B = int_t(n,a,0)
    Fv = sp.expand(B*L1s**b)
    acc = 0
    for i in range(a+1):
        ci = sp.expand(B).coeff(L0,i)
        if ci == 0: continue
        ci_t = sp.simplify(ci/t**(n+1))
        for j in range(n+2):
            cj = (-1)**j * sp.binomial(n+1,j)
            if j == 0:
                term = cj * (-int_s_over(i, b-1))
            else:
                term = cj * int_polylogweight(j-1, i, b-1)
            acc += ci_t*term
    return sp.expand(Fv + b*acc)
z2s = sp.pi**2/6; z3s = sp.zeta(3)
def texpo(term, sym): return int(term.as_powers_dict().get(sym,0))
def limit1(Fv):
    e = Fv.subs(t, 1-s_)
    e = e.subs({sp.polylog(2,s_):0, sp.polylog(3,s_):0, sp.polylog(2,1-s_):z2s,
                sp.polylog(3,1-s_):z3s, sp.log(1-s_):0})
    e = sp.expand(e)
    keep = sp.expand(sum(term for term in sp.Add.make_args(e) if texpo(term,s_)==0))
    ls = sp.log(s_)
    bad = []
    for q in range(1,7):
        c = sp.simplify(keep.coeff(ls,q))
        if c != 0: bad.append((q,c))
    return keep, bad
# test on the keys
O5 = '<project>/research/08_theory/OCT6_octagon_closed/oct5'
JlocO = pickle.load(open(O5+'/Jloc.pkl','rb'))
for key in sorted(JlocO.keys()):
    n,a,b = key
    Fv = int_t(n,a,b)
    keep, bad = limit1(Fv)
    if bad:
        print("KEY", key, "bad coeffs:", bad[:2])
        break
else:
    print("all keys OK in this standalone replica")
