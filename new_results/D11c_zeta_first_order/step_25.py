step_06 = r'''# step_06.py -- D11c(3): geometric zeta(3) formulas, 60-70 digit verification
import mpmath as mp
import sympy as sp
mp.mp.dps = 70
z3 = mp.zeta(3); pi = mp.pi

# (i) perimeter third derivative: P(p) = 4(1/p) G(u), u=1/p-1, G(u)=int sqrt((1-t)^{2u}+t^{2u}) dt
# P^(3)(1)/sqrt2 = int_0^1 [ -1/2 ln^3 S - 3 lnS l^2 - 9 ln^2 S - 18 l^2 - 24 lnS ] dt
I3 = mp.quad(lambda t: (-mp.mpf('0.5')*mp.log(t*(1-t))**3 - 3*mp.log(t*(1-t))*mp.log(t/(1-t))**2
                        - 9*mp.log(t*(1-t))**2 - 18*mp.log(t/(1-t))**2 - 24*mp.log(t*(1-t))),
             [0, mp.mpf('0.5'), 1])
print("(i) zeta(3) = pi^2/3 + P^(3)(1)/(6 sqrt2):")
print(mp.nstr(pi**2/3 + I3/6, 55)); print("error:", mp.nstr(abs(pi**2/3+I3/6-z3), 6))

# (ii) simplified J form
J = mp.quad(lambda t: mp.log(t)**2*mp.log(1-t), [0,1])
print("\n(ii) zeta(3) = 3 - pi^2/6 + (1/2) int_0^1 ln^2(t) ln(1-t) dt:")
print(mp.nstr(3-pi**2/6+J/2, 55)); print("error:", mp.nstr(abs(3-pi**2/6+J/2-z3), 6))

# (iii) J as series and classical evaluation
Jser = mp.nsum(lambda k: 1/(k*(k+1)**3), [1, mp.inf])
print("\n(iii) J =", mp.nstr(J, 30), "= -2 sum 1/(k(k+1)^3) =", mp.nstr(-2*Jser, 30))
print("      2 zeta(3) - 6 + pi^2/3 =", mp.nstr(2*z3-6+pi**2/3, 30))

# (iv) double-integral form
DI = mp.quad(lambda t,s: t*mp.log(t)**2/(1-t+t*s) if t>0 else mp.mpf(0), [0,1],[0,1])
print("\n(iv) zeta(3) = 3 - pi^2/6 - (1/2) int int t ln^2 t/(1-t+ts) ds dt:")
print(mp.nstr(3-pi**2/6-DI/2, 55)); print("error:", mp.nstr(abs(3-pi**2/6-DI/2-z3), 6))

# (v) symmetric form
SI = mp.quad(lambda t: mp.log(t)*mp.log(1-t)*mp.log(t*(1-t)) if 0<t<1 else mp.mpf(0), [0,mp.mpf('0.5'),1])
print("\n(v) zeta(3) = 3 - pi^2/6 + (1/4) int ln t ln(1-t) ln(t(1-t)) dt:")
print(mp.nstr(3-pi**2/6+SI/4, 55)); print("error:", mp.nstr(abs(3-pi**2/6+SI/4-z3), 6))

# (vi) independent check: 9-point finite difference of direct arclength quadrature
def perim(p):
    def f(th):
        c, s = mp.cos(th), mp.sin(th)
        dxd = (2/p)*c**(2/p-1)*(-s); dyd = (2/p)*s**(2/p-1)*c
        return mp.sqrt(dxd**2+dyd**2)
    return 4*mp.quad(f, [0, mp.pi/2])
N = 4
B = sp.Matrix(2*N+1, 2*N+1, lambda m,k: sp.Integer(1) if (k-N)==0 and m==0 else (sp.Integer(0) if (k-N)==0 else sp.Integer(k-N)**m))
bs = sp.Matrix([0,0,0,6,0,0,0,0,0]); cs = B.inv()*bs
h = mp.mpf('0.004')
fd = sum(sp.Rational(cs[k])*perim(1+(k-N)*h) for k in range(2*N+1) if cs[k]!=0)/h**3
print("\n(vi) P^(3)(1) 9-point FD:", mp.nstr(fd, 20),
      " vs 6 sqrt2 (zeta3-pi^2/3):", mp.nstr(6*mp.sqrt(2)*(z3-pi**2/3), 20))

# (vii) classical integrals for comparison
print("\n(vii) (1/2) int ln^2 t/(1-t) dt =", mp.nstr(mp.quad(lambda t: mp.log(t)**2/(1-t),[0,1])/2, 30))
print("      Beukers (1/2) int int -ln(xy)/(1-xy) =", mp.nstr(mp.quad(lambda x,y: -mp.log(x*y)/(1-x*y) if x*y>0 else mp.mpf(0),[0,1],[0,1])/2, 30))
'''
open('step_06.py','w').write(step_06)

step_07 = r'''# step_07.py -- D11c(1)+(2): generating function of kappa_n: V-recurrence, asymptotics, Borel sum
import mpmath as mp
mp.mp.dps = 50
R = lambda t: mp.sqrt((1+mp.exp(4*t))/2)
V = {n: mp.diff(R, 0, n) for n in range(0, 20)}

# (a) recurrence from R R' = 2 R^2 - 1
ok = all(abs(sum(mp.binomial(n,k)*V[k]*V[n+1-k] for k in range(0,n+1))
             - (2*sum(mp.binomial(n,k)*V[k]*V[n-k] for k in range(0,n+1)) - (1 if n==0 else 0))) < mp.mpf('1e-35')
         for n in range(0, 15))
print("(a) recurrence sum_k C(n,k) V_k V_{n+1-k} = 2 sum_k C(n,k) V_k V_{n-k} - delta_{n0}:", ok)

# (b) asymptotics of kappa_n
kappa = {3: mp.mpf(1), 5: mp.mpf('-0.75'), 7: mp.mpf('8.0625'), 9: -mp.mpf(22805)/192,
         11: mp.mpf(805173)/256, 13: -mp.mpf(131893755)/1024, 15: mp.mpf(154940467501)/20480}
print("(b) kappa_{2m+1} ~ -4 (-1)^m (2m)! (2/pi)^{2m+1} (2m+1)^{-3/2}:")
for m in range(1, 8):
    n = 2*m+1
    pred = -4*(-1)**m*mp.factorial(2*m)*(2/mp.pi)**n*mp.mpf(n)**(-1.5)
    print(f"  n={n:>2}: kappa={mp.nstr(kappa[n],12)} pred={mp.nstr(pred,12)} ratio={mp.nstr(kappa[n]/pred,8)}")

# (c) Borel sum of G(x)=sum kappa_n x^n at x=0.1 vs optimal truncation
def S(x): return mp.quad(lambda u: (mp.exp(u)-1)/u if u!=0 else mp.mpf(1), [0, x])
def Rodd(v): return (R(v)-R(-v))/2
def Phi(xi):
    E_ = (S(xi)-S(-xi))/2 - xi
    integ = mp.quad(lambda v: Rodd(v)/v if v!=0 else mp.mpf(1), [0, xi/2])
    return -4*E_ - 4*xi + 8*integ
x = mp.mpf('0.1')
Gb = mp.quad(lambda t: mp.e**(-t)*Phi(x*t), [0, mp.inf])
part = 0; best = None; bestterm = mp.mpf('1e9')
for n in range(3, 60, 2):
    kn = 4*(-1)**n/n*(1-2**(1-n)*mp.diff(R,0,n)); part += kn*x**n
    if abs(kn*x**n) < bestterm: bestterm = abs(kn*x**n); best = (n, part)
print("\n(c) Borel sum G(0.1) =", mp.nstr(Gb, 25))
print("    optimal truncation n=", best[0], ":", mp.nstr(best[1], 25), " least term", mp.nstr(bestterm, 6))

# (d) closed GF identity at small x: G(x) = 4x - 4 atanh x + 8 sum_{odd n>=3} V_n (x/2)^n / n
x2 = mp.mpf('0.001')
direct = 4*x2 - 4*mp.atanh(x2) + 8*sum(mp.diff(R,0,n)*(x2/2)**n/n for n in range(3, 30, 2))
via_k = sum((4*(-1)**n/n*(1-2**(1-n)*mp.diff(R,0,n)))*x2**n for n in range(3, 30, 2))
print("\n(d) G(0.001) direct:", mp.nstr(direct, 30))
print("    G(0.001) via kappa series:", mp.nstr(via_k, 30))
'''
open('step_07.py','w').write(step_07)
print("step_06.py, step_07.py written OK")
import os
print([f for f in sorted(os.listdir('.')) if f.startswith('step')])
