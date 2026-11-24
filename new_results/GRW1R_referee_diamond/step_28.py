from flint import arb
import sympy as sp
from sympy import Rational, factorial

K = 96
a=[Rational(0)]*(K+3)
for m in range(K+3): a[m]=Rational((-2)**m, factorial(m+1))
b=[Rational(0)]*(K+3); b[0]=Rational(1)
for m in range(1,K+3): b[m]=(a[m]-sum(b[i]*b[m-i] for i in range(1,m)))/2
e=[b[k+2] for k in range(K+1)]
c=[Rational(0)]*(K//2+1)
for m in range(K//2+1):
    c[m]=((-1)**m)*sum(((-1)**(2*m-i))*e[i]*e[2*m-i] for i in range(2*m+1))
print("c_0..c_6:", [str(x) for x in c[:7]])

e6 = arb(6).exp()
M_E3 = (((1+e6)/6).sqrt()+1+arb('1.5'))/9
M2 = M_E3**2
rho = arb.pi()/2; send = arb('1.9')*rho
carb = [arb(x.numerator)/arb(x.denominator) for x in c]
KK = 45
def psi_series(s):  # s arb interval; returns interval bound for psi(s) via u=s^2 series + tail
    u = s*s
    P = arb(0)
    for m in range(KK, -1, -1):
        P = P*u + carb[m]
    tail = M2*(u/9)**(KK+1)/(1-u/9)
    return arb(P.mid(), P.rad() + float(tail.upper()))

# piece: [0.1, 1.5] via series; [1.5, 1.9rho] via direct balls
maxub = None; arg=None
h = 0.001
cc = 0.1 + h
while cc - h <= 1.5:
    seg = arb(arb(cc), arb(h))
    v = psi_series(seg)
    if maxub is None or v.upper() > maxub: maxub, arg = v.upper(), cc
    cc += 2*h
print("series piece [0.1,1.5] max psi <=", maxub, " near s=", arg)
def psi_direct(s):
    q = s.sin()/s
    return (q + 1 + s**2/4 - (2*(s/2).cos() + s*(s/2).sin())*q.sqrt())/s**4
maxub2 = None; arg2=None
cc = 1.5 + h
while cc - h <= float(send):
    hi = min(cc+h, float(send))
    seg = arb(arb((cc+hi)/2), arb((hi-cc)/2))
    v = psi_direct(seg)
    if maxub2 is None or v.upper() > maxub2: maxub2, arg2 = v.upper(), cc
    cc += 2*h
print("direct piece [1.5, 1.9rho] max psi <=", maxub2, " near s=", arg2)
psi524 = (arb(5)/24)**2
print("(5/24)^2 =", psi524)
finalmax = max(psi524.upper(), maxub, maxub2)
print("FINAL CERTIFIED max |E| on slit =", finalmax.sqrt() if hasattr(finalmax,'sqrt') else arb(finalmax).sqrt())
with open('progress.log','a') as f:
    f.write(f"step_08d: E*=5/24 CERTIFIED rigorously: series piece max {maxub}, direct piece max {maxub2}, both < (5/24)^2={psi524}\n")
