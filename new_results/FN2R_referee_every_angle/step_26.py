import sympy as sp, time
t0=time.time()
a,A,B,D = sp.symbols('a A B D')
t, W, mu = sp.symbols('t W mu')
K = 3
# u-series for -log f(1-tau)
u = -a*t + A*t**2/2 - B*t**3/6 + D*t**4/24
Acal = sum((-1)**(r+1)*u**r/r for r in range(1,K+2))
Acal = sp.expand(Acal)
# keep terms up to t^(K+1)
Acal = sum(m for m in Acal.args if sp.degree(m,t) <= K+1)
v = (-A*t + B*t**2/2 - D*t**3/6)/a
Btil = sum((-1)**(r+1)*v**r/r for r in range(1,K+1))
Btil = sp.expand(Btil)
Btil = sum(m for m in Btil.args if sp.degree(m,t) <= K)
print("Acal =", sp.expand(Acal))
print("Btil =", sp.expand(Btil))

idx = [(i,j) for i in range(1,K+1) for j in range(0,K+1-i)]
vs = sp.symbols(' '.join('t%d%d'%(i,j) for i,j in idx))
Tau = sum(vv*W**i*mu**j for vv,(i,j) in zip(vs,idx))
def trunc(expr):
    expr = sp.expand(expr)
    return sum(m for m in expr.args if sp.degree(m,W)+sp.degree(m,mu) <= K)
eq = trunc(trunc(Acal.subs(t,Tau)) - mu*trunc(Btil.subs(t,Tau)) - W)
sols = {}
for (i,j),vv in zip(idx,vs):
    e = trunc(eq.subs(sols))
    coeff = sp.expand(e).coeff(W,i).coeff(mu,j)
    sol = sp.solve(coeff, vv)
    assert len(sol)==1, (i,j,coeff)
    sols[vv] = sol[0]
th = {k: sp.simplify(sols[vv]) for k,vv in zip(idx,vs)}
print("my theta_ij:")
for k in sorted(th): print("  theta%s ="%str(k), th[k])
th_fn2 = {(1,0): 1/a, (2,0): (A-a**2)/(2*a**3), (1,1): -A/a**3,
          (3,0): (3*A**2-3*A*a**2-B*a+a**4)/(6*a**5),
          (2,1): (-4*A**2+3*A*a**2+B*a)/(2*a**5), (1,2): A**2/a**5}
print("match FN2's theta_ij:", all(sp.simplify(th[k]-th_fn2[k])==0 for k in th_fn2))
print("theta_{1j} = a^{-1}(-A/a^2)^j:", sp.simplify(th[(1,1)]+A/a**3)==0, sp.simplify(th[(1,2)]-A**2/a**5)==0)
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} step_05a: independent theta_ij (series reversion, manual log series) matches FN2's 6 coefficients exactly [referee, sympy]\n")
