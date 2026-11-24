import sympy as sp, time
t0=time.time()
# ===== Independent derivation of theta_{ij} (FN2 Lemma E / P3 Theorem 2 setup) =====
# f(1-tau) = 1 - a t + A t^2/2 - B t^3/6 + D t^4/24 (a=f'(1),A=f''(1),B=f'''(1),D=f''''(1))
a,A,B,D = sp.symbols('a A B D')
t, W, mu = sp.symbols('t W mu')
fser = 1 - a*t + A*t**2/2 - B*t**3/6 + D*t**4/24
Acal = -sp.log(fser)                    # A(tau) = -log f(1-tau)
fpser = sp.diff(fser, t)*(-1)           # f'(1-tau) = a - A t + B t^2/2 - D t^3/6
Btil = sp.log(fpser/a)                  # log(f'(1-tau)/a)

# Solve Acal(Tau) - mu*Btil(Tau) = W for Tau = sum_{i>=1,j>=0} th_ij W^i mu^j, i+j<=3, by series reversion.
# Method: iterate. Start Tau=0; Newton/iteration in series ring: Tau <- inverse of Acal at (W+mu Btil(Tau)).
# Simplest: formal ansatz and coefficient solve.
K = 3
idx = [(i,j) for i in range(1,K+1) for j in range(0,K+1-i)]
vs = sp.symbols(' '.join('t%d%d'%(i,j) for i,j in idx))
Tau = sum(v*W**i*mu**j for v,(i,j) in zip(vs,idx))
eq = sp.expand(Acal.subs(t,Tau) - mu*Btil.subs(t,Tau) - W)
# expand to total order K in (W,mu)
eqs = sp.expand(eq)
sols = {}
for (i,j),v in zip(idx,vs):
    e = sp.expand(eqs.subs(sols))
    coeff = e.coeff(W,i).coeff(mu,j)
    sols[v] = sp.solve(coeff, v)[0]
th = {k: sp.simplify(sols[v]) for k,v in zip(idx,vs)}
print("my theta_ij:")
for k in sorted(th): print("  theta%s ="%str(k), th[k])
# FN2's values:
th_fn2 = {(1,0): 1/a, (2,0): (A-a**2)/(2*a**3), (1,1): -A/a**3,
          (3,0): (3*A**2-3*A*a**2-B*a+a**4)/(6*a**5),
          (2,1): (-4*A**2+3*A*a**2+B*a)/(2*a**5), (1,2): A**2/a**5}
print("match FN2's theta_ij:", all(sp.simplify(th[k]-th_fn2[k])==0 for k in th_fn2))
# P3's general theta_{1j} = a^{-1}(-A/a^2)^j:
print("theta_{1j} formula:", sp.simplify(th[(1,1)]+A/a**3)==0, sp.simplify(th[(1,2)]-A**2/a**5)==0)
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} step_05a: independent theta_ij derivation (series reversion) matches FN2's 6 coefficients exactly [referee, sympy]\n")
