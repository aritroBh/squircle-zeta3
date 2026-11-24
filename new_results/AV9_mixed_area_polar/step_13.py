
# step_13.py : (i) c_7/c_8 PSLQ at dps 250/350; (ii) top-weight theorem check; (iii) V(eta)=V(1/eta) numeric
import mpmath as mp
import sympy as sp
import pickle
import av9core as C

c_exact = {int(n): sp.sympify(v) for n,v in pickle.load(open('scratchpad/c_exact_selfdual.pkl','rb')).items()}
def basis_elements(N):
    out = []
    for a in range(N+1):
        for b in range(N//2+1):
            for c in range(N//3+1):
                for dd in range(N//5+1):
                    for e in range(N//7+1):
                        if 1 <= a+2*b+3*c+5*dd+7*e <= N: out.append((a,b,c,dd,e))
    return out
for n, dps in [(6,250),(7,250),(8,350)]:
    mp.mp.dps = dps
    L2 = mp.log(2); z = {k: mp.zeta(k) for k in (2,3,5,7)}
    target = mp.mpf(sp.N(c_exact[n], dps+10))
    basis = basis_elements(n)
    vals = [L2**m[0]*z[2]**m[1]*z[3]**m[2]*z[5]**m[3]*z[7]**m[4] for m in basis]
    rel = mp.pslq([target]+vals, maxcoeff=10**9, maxsteps=20000)
    if rel is None:
        print(f"c_{n} dps={dps}: FAILED"); continue
    resid = abs(rel[0]*target + sum(r*v for r,v in zip(rel[1:],vals)))
    # exact vector for comparison
    Z2 = sp.Symbol('Z2'); dpi = sp.Symbol('dpi')
    e = sp.expand(sp.expand(c_exact[n]).subs(sp.pi, dpi).subs(dpi**2, 6*Z2))
    vec = []
    for m in basis:
        cf = sp.expand(e).coeff(Z2,m[1]).coeff(sp.log(2),m[0]).coeff(sp.zeta(3),m[2]).coeff(sp.zeta(5),m[3]).coeff(sp.zeta(7),m[4])
        vec.append(cf)
    i0 = next(i for i,v in enumerate(vec) if v)
    rat = sp.Rational(rel[1:][i0], vec[i0])
    match = all(sp.Rational(r) == rat*v for r,v in zip(rel[1:],vec))
    print(f"c_{n} dps={dps}: resid={mp.nstr(resid,3)} max|c|={max(abs(r) for r in rel)} dim={len(basis)+1} match={match}")

# (ii) top-weight theorem: degree-n part of e_n equals (lv^n+lw^n)/(2 n!) for self-dual
print("\n-- top-weight check (self-dual) --")
h, e = C.selfdual_coeffs(8)
for n in range(1,9):
    en = sp.expand(e[n])
    deg_n = sum(cf*C.lv**a*C.lw**b for (a,b),cf in sp.Poly(en, C.lv, C.lw).terms() if a+b==n)
    pred = (C.lv**n + C.lw**n)/(2*sp.factorial(n))
    print(f"n={n}: deg-n part of e_n == (lv^n+lw^n)/(2 n!) : {sp.simplify(deg_n-pred)==0}")

# (iii) duality numeric
mp.mp.dps = 60
def Vnum(etaval):
    q = 1+etaval
    f = lambda v: (2/q)*(v*(1-v))**(1/q-1)*(v**(q-1)+(1-v)**(q-1))**(1/q)
    return mp.quad(f, [0, mp.mpf('0.5'), 1])
print("\n-- V(eta) vs V(1/eta) --")
for ev in [mp.mpf('0.25'), mp.mpf('0.1'), mp.mpf('0.05'), mp.mpf('3')]:
    a, b = Vnum(ev), Vnum(1/ev)
    print(f"eta={ev}: V(eta)={mp.nstr(a,25)} V(1/eta)={mp.nstr(b,25)} diff={mp.nstr(abs(a-b),3)}")
