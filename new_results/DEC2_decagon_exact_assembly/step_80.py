import sympy as sp, pickle, os, sys, time
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, '<project>/research/08_theory/DEC1_decagon')
t = sp.Symbol('t'); r5 = sp.sqrt(5); phi = (1+r5)/2
kap1 = 1/phi; kap4 = -1/phi**2
u1 = 1+kap1*t; u4 = 1+kap4*t
w1 = sp.log(u1); w4 = sp.log(u4)
L0 = sp.log(t); L1t = sp.log(1-t)
import dec_engine as E
Jtab = pickle.load(open('Jtab_dec2.pkl','rb'))
Lf = sp.Symbol('Lf')
def K1v(n,a,b,c,side):
    """exact single-ratio value at side (1: x1/w1, 4: x4/w4), in engine's named atoms."""
    sub = E.sub1 if side==1 else E.sub4
    return sp.expand(Jtab[(n,a,b,c)].subs(sub, simultaneous=True))
D = pickle.load(open('D_all_proved.pkl','rb'))
def canon5(e):
    e = sp.expand(e)
    P = sp.Poly(e, *E.NAMED)
    out = sp.Integer(0)
    for monom, coeff in P.terms():
        c = sp.together(coeff)
        cc = c.subs(r5, -r5)
        p = sp.radsimp(sp.expand((c+cc)/2)); q = sp.radsimp(sp.expand((c-cc)/(2*r5)))
        mm = sp.Integer(1)
        for a_, e_ in zip(E.NAMED, monom): mm *= a_**e_
        out += (p + q*r5)*mm
    return sp.expand(out)
# ---- partial-fraction identities (machine check, all n, both kappas) ----
def pf_check(n, kap):
    lhs = t**(n+1)
    ssum = sum((-1/kap)**j * t**(n-j) for j in range(n+1))/kap
    rem = (-1/kap)**(n+1)
    return sp.simplify(lhs - (ssum + rem/(1+kap*t))*(1+kap*t))
def pf_check_minus(n, kap):
    lhs = t**(n+1) - 1
    ssum = sum((-1/kap)**j * t**(n-j) for j in range(n+1))/kap
    rem = (-1/kap)**(n+1) - 1
    return sp.simplify(lhs - (ssum + rem/(1+kap*t))*(1+kap*t))
pf_ok = all(pf_check(n,kap)==0 for n in range(4) for kap in (kap1,kap4))
pfm_ok = all(pf_check_minus(n,kap)==0 for n in range(4) for kap in (kap1,kap4))
onesided = all(sp.simplify((t**(n+1)-1) + sum(t**j for j in range(n+1))*(1-t))==0 for n in range(4))
print("partial-fraction identities: pf", pf_ok, " pf-minus", pfm_ok, " (t^{n+1}-1)/(1-t)", onesided)
# ---- IBP integrand identities (machine check) ----
def ibp_check(n, typ):
    if typ=='W':
        C = t**(n+1)*w1*w4/(n+1)
        R = -kap1*t**(n+1)*w4/u1 - kap4*t**(n+1)*w1/u4
    elif typ=='M1':
        C = t**(n+1)*L0*w1*w4/(n+1)
        R = -t**n*w1*w4 - kap1*t**(n+1)*L0*w4/u1 - kap4*t**(n+1)*L0*w1/u4
    elif typ=='M2':
        C = (t**(n+1)-1)*L1t*w1*w4/(n+1)
        R = sum(t**j for j in range(n+1))*w1*w4 - kap1*(t**(n+1)-1)*L1t*w4/u1 - kap4*(t**(n+1)-1)*L1t*w1/u4
    elif typ=='M3':
        C = t**(n+1)*w1**2*w4/(n+1)
        R = -2*kap1*t**(n+1)*w1*w4/u1 - kap4*t**(n+1)*w1**2/u4
    elif typ=='M4':
        C = t**(n+1)*w1*w4**2/(n+1)
        R = -kap1*t**(n+1)*w4**2/u1 - 2*kap4*t**(n+1)*w1*w4/u4
    integrand = {'W': t**n*w1*w4, 'M1': t**n*L0*w1*w4, 'M2': t**n*L1t*w1*w4,
                 'M3': t**n*w1**2*w4, 'M4': t**n*w1*w4**2}[typ]
    lhs = (n+1)*integrand
    rhs = (n+1)*sp.diff(C,t) + R
    return sp.simplify(sp.expand(lhs - rhs)) == 0
ibp_ok = all(ibp_check(n,typ) for n in range(4) for typ in ['W','M1','M2','M3','M4'])
print("IBP integrand identities (20):", ibp_ok)
# ---- boundary terms ----
bnd = {'W': (w1*w4).subs(t,1), 'M3': (w1**2*w4).subs(t,1), 'M4': (w1*w4**2).subs(t,1)}
print("boundaries: W =", sp.simplify(bnd['W'].subs(sp.log(phi), E.lf)),
      " M3 =", sp.simplify(bnd['M3'].subs(sp.log(phi), E.lf)),
      " M4 =", sp.simplify(bnd['M4'].subs(sp.log(phi), E.lf)),
      " (M1, M2 boundaries are 0: t^{n+1} L0 -> 0 at 0 and L0(1)=0; (t^{n+1}-1)L1 -> 0 at 1, L1(0)=0)")
