
# step_06.py : genuinely broken family r = q^c (c=1 self-dual). u = 1 + 1/((1+eta)^c - 1)
# W^c(eta) = 2(1+eta)^{-c} int (v(1-v))^{-1+(1+eta)^{-c}} (v^{(1+eta)(1-(1+eta)^{-c})}+(1-v)^{...})^{1/(1+eta)} dv
import sympy as sp
import av9core as C
import time, pickle
t0=time.time()
eta, c = sp.symbols('eta c')
N = 5
def ser(expr):
    return sp.series(expr, eta, 0, N+1).removeO().expand()
alpha  = ser(-1 + (1+eta)**(-c))
beta_  = ser((1+eta)*(1-(1+eta)**(-c)))
gamma  = ser((1+eta)**(-1))
P      = 2*(1+eta)**(-c)
lnI, e = C.general_coeffs(N, alpha, beta_, gamma, P)
tab = C.beta_int_table(N)
h2, e2 = C.selfdual_coeffs(N)
c_ref = {n: C.integrate_poly(4*e2[n], tab) for n in range(1, N+1)}
d = {}
ok=True
for n in range(1, N+1):
    d[n] = sp.simplify(C.integrate_poly(e[n], tab))
    chk = sp.simplify(d[n].subs(c,1) - c_ref[n])
    if chk!=0: ok=False
    print(f"d_{n}(c) done t={time.time()-t0:.1f}  c=1 check: {chk}", flush=True)
print("CONSISTENCY:", ok)
for n in range(1, N+1):
    print(f"\nd_{n}(c) =", sp.factor(d[n]))
for n in [2,3,4,5]:
    z3c = sp.simplify(sp.expand(d[n]).coeff(sp.zeta(3),1))
    print(f"[zeta(3)] d_{n}(c) =", sp.factor(z3c))
pickle.dump({n: str(d[n]) for n in d}, open('scratchpad/d_rqc_family.pkl','wb'))
