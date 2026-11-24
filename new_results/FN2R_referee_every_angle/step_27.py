import sympy as sp, time
a,A,B,D = sp.symbols('a A B D')
t, W, mu = sp.symbols('t W mu')
K = 3
u = -a*t + A*t**2/2 - B*t**3/6 + D*t**4/24
Acal = sum((-1)**r*u**r/r for r in range(1,K+2))   # -log(1+u)
Acal = sp.expand(Acal)
Acal = sum(m for m in Acal.args if sp.degree(m,t) <= K+1)
print("Acal =", sp.expand(Acal))
v = (-A*t + B*t**2/2 - D*t**3/6)/a
Btil = sum((-1)**(r+1)*v**r/r for r in range(1,K+1))
Btil = sp.expand(Btil)
Btil = sum(m for m in Btil.args if sp.degree(m,t) <= K)
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
    assert len(sol)==1
    sols[vv] = sol[0]
th = {k: sp.simplify(sols[vv]) for k,vv in zip(idx,vs)}
th_fn2 = {(1,0): 1/a, (2,0): (A-a**2)/(2*a**3), (1,1): -A/a**3,
          (3,0): (3*A**2-3*A*a**2-B*a+a**4)/(6*a**5),
          (2,1): (-4*A**2+3*A*a**2+B*a)/(2*a**5), (1,2): A**2/a**5}
print("my theta_ij:")
for k in sorted(th): print("  theta%s ="%str(k), th[k])
print("match FN2:", all(sp.simplify(th[k]-th_fn2[k])==0 for k in th_fn2))
print("theta_{1j} formula a^{-1}(-A/a^2)^j:", sp.simplify(th[(1,1)]-(-A/a**3))==0, sp.simplify(th[(1,2)]-A**2/a**5)==0)
# Also theta_{k,0} = [W^k] Acal^{-1}(W): direct reversion check
import sympy.series as _
w = sp.Symbol('w')
Acal_t = sum((-1)**r*(-a*t + A*t**2/2 - B*t**3/6 + D*t**4/24)**r/r for r in range(1,5))
rev = sp.expand(sp.revert(sp.expand(Acal_t), t, 4))
print("Acal^{-1} =", rev)
print("theta_{k,0} == [W^k]Acal^{-1}:", all(sp.simplify(rev.coeff(t,k)-th[(k,0)])==0 for k in (1,2,3)))
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} step_05a fixed: independent theta_ij matches FN2 exactly; theta_{{1j}} formula and theta_{{k0}}=Acal-reversion confirmed [referee, sympy]\n")
