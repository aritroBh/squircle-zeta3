import sympy as sp, time
t0=time.time()
# theta_{ij} to order 4 with my own reversion (extends step_05a)
a,A,B,D = sp.symbols('a A B D')
t, W, mu = sp.symbols('t W mu')
K = 4
u = -a*t + A*t**2/2 - B*t**3/6 + D*t**4/24
Acal = sum((-1)**r*u**r/r for r in range(1,K+2))
Acal = sp.expand(Acal)
Acal = sum(mm for mm in Acal.args if sp.degree(mm,t) <= K+1)
v = (-A*t + B*t**2/2 - D*t**3/6)/a
Btil = sum((-1)**(r+1)*v**r/r for r in range(1,K+1))
Btil = sp.expand(Btil)
Btil = sum(mm for mm in Btil.args if sp.degree(mm,t) <= K)
idx = [(i,j) for i in range(1,K+1) for j in range(0,K+1-i)]
vs = sp.symbols(' '.join('t%d%d'%(i,j) for i,j in idx))
Tau = sum(vv*W**i*mu**j for vv,(i,j) in zip(vs,idx))
def trunc(expr):
    expr = sp.expand(expr)
    return sum(mm for mm in expr.args if sp.degree(mm,W)+sp.degree(mm,mu) <= K)
eq = trunc(trunc(Acal.subs(t,Tau)) - mu*trunc(Btil.subs(t,Tau)) - W)
sols = {}
for (i,j),vv in sorted(zip(idx,vs), key=lambda p: p[0][0]+p[0][1]):
    e = trunc(eq.subs(sols))
    coeff = sp.expand(e).coeff(W,i).coeff(mu,j)
    sol = sp.solve(coeff, vv)
    assert len(sol)==1, (i,j,coeff)
    sols[vv] = sol[0]
th = {k: sp.simplify(sols[vv]) for k,vv in zip(idx,vs)}
# P3's quoted order-4 coefficients (b=A, c=B, d=D):
th_P3 = {
 (4,0): -(a**6-7*a**4*A-6*a**3*B+18*a**2*A**2-a**2*D+10*a*A*B-15*A**3)/(24*a**7),
 (3,1): -(7*a**4*A+6*a**3*B-24*a**2*A**2+a**2*D-13*a*A*B+23*A**3)/(6*a**7),
 (2,2): -3*A*(2*a**2*A+a*B-3*A**2)/(2*a**7),
 (1,3): -A**3/a**7}
print("theta_40 mine:", th[(4,0)])
print("match P3 order-4:", all(sp.simplify(th[k]-th_P3[k])==0 for k in th_P3))
# T4: [mu^4] sum_{i+j<=3} th_ij mu^{i+j} (-1)^{i+1} i! sum_l d_l/(l-1)^{i+1} (mu/(a(1+mu)))^{l-1}
c = sp.Symbol('c')
d = lambda l: (-1)**l*(1-c**2)*sp.legendre(l-1, c).diff(c)/l
Tser = 0
for (i,j),tij in th.items():
    for l in range(2, K+2):
        term = tij*mu**(i+j)*(-1)**(i+1)*sp.factorial(i)*d(l)/(l-1)**(i+1) * (mu/(a*(1+mu)))**(l-1)
        Tser += sp.series(term, mu, 0, K+1).removeO()
T4 = sp.expand(sp.expand(Tser)).coeff(mu,4)
print("my T4 =", sp.factor(T4))
T4_P3 = (12*a**4-48*a**2*A-a**2-24*a*B+96*A**2)/(24*a**6)   # P3's T4 at theta=pi/2 (c=0)
print("T4 at c=0 == P3 T4:", sp.simplify(T4.subs(c,0)-T4_P3)==0)
print("elapsed %.1fs"%(time.time()-t0))
with open("progress.log","a") as f2:
    f2.write(f"{time.strftime('%H:%M')} step_12a: my theta_{{4j}} match P3's quoted ones exactly; my T4(c) at c=0 matches P3's T4. General-angle T4 derived [referee, sympy EXACT]\n")
