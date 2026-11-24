
# step_05.py : broken symmetry family W_c(eta), FIXED machinery. c=1 must reproduce self-dual c_n.
import sympy as sp
import av9core as C
import time, pickle
t0=time.time()
eta, c = sp.symbols('eta c')
N = 6
alpha  = sp.series(-eta/(eta+c), eta, 0, N+1).removeO()
beta_  = sp.series(eta*(1+eta)/(eta+c), eta, 0, N+1).removeO()
gamma  = sp.series(1/(1+eta), eta, 0, N+1).removeO()
P      = 2*c/(eta+c)
lnI, e = C.general_coeffs(N, alpha, beta_, gamma, P)
tab = C.beta_int_table(N)
# self-dual reference
h2, e2 = C.selfdual_coeffs(N)
c_ref = {n: C.integrate_poly(4*e2[n], tab) for n in range(1, N+1)}
d = {}
ok = True
for n in range(1, N+1):
    d[n] = sp.simplify(C.integrate_poly(e[n], tab))
    chk = sp.simplify(d[n].subs(c,1) - c_ref[n])
    if chk != 0: ok=False
    print(f"d_{n}(c) done t={time.time()-t0:.1f}  d_{n}(1)-c_{n} = {chk}", flush=True)
print("CONSISTENCY c=1 vs self-dual:", ok)
for n in range(1, N+1):
    print(f"\nd_{n}(c) =", sp.factor(d[n]))
for n in [3,4,5,6]:
    z3c = sp.simplify(sp.expand(d[n]).coeff(sp.zeta(3),1))
    print(f"[zeta(3)] d_{n}(c) =", sp.factor(z3c))
z5c = sp.simplify(sp.expand(d[5]).coeff(sp.zeta(5),1)); print(f"[zeta(5)] d_5(c) =", sp.factor(z5c))
z56 = sp.simplify(sp.expand(d[6]).coeff(sp.zeta(5),1)); print(f"[zeta(5)] d_6(c) =", sp.factor(z56))
pickle.dump({n: str(d[n]) for n in d}, open('scratchpad/d_cfamily.pkl','wb'))
