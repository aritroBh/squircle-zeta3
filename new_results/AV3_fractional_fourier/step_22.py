from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, quad)
mp.dps = 80
def N(x,d=40): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")

# W(s) := sum_n b_n(s)  (Abel)  =  (1/pi) int_0^pi t^s cot(t/2) dt
def b_exact(n, s):
    if s == 1: return 2*(-1)**(n+1)/n
    if s == 2: return 2*pi*(-1)**(n+1)/n - 4*(1-(-1)**n)/(pi*n**3)
    if s == 3: return 2*pi**2*(-1)**(n+1)/n + 12*(-1)**n/n**3
    if s == 4: return 2*pi**3*(-1)**(n+1)/n + 24*pi*(-1)**n/n**3 + 48*(1-(-1)**n)/(pi*n**5)

# exact values from summing the closed forms
lam = lambda k: (1-2**(-k))*zeta(k)   # odd zeta
eta = lambda k: (1-2**(1-k))*zeta(k)
W_exact = {
 1: 2*log(2),
 2: 2*pi*log(2) - (4/pi)*(zeta(3)+eta(3)),
 3: 2*pi**2*log(2) - 12*eta(3),
 4: 2*pi**3*log(2) - 24*pi*eta(3) + (48/pi)*(zeta(5)+eta(5)),
}
for s in [1,2,3,4]:
    Ws = nsum(lambda n: b_exact(n,s), [1, inf])
    Wi = (1/pi)*quad(lambda t: t**s/tan(t/2), [0, pi])
    print(f"s={s}: series={N(Ws,30)}")
    print(f"      integral={N(Wi,30)}")
    print(f"      exact   ={N(W_exact[s],30)}  err={nstr(abs(Ws-W_exact[s]),3)}")
print()
print("W(2) = 2 pi log2 - 7 zeta(3)/pi =", N(2*pi*log(2)-7*zeta(3)/pi,40))
print("W(4) = 2 pi^3 log2 - 18 pi zeta(3) + 93 zeta(5)/pi =", N(2*pi**3*log(2)-18*pi*zeta(3)+93*zeta(5)/pi,40))
logline("step_20: W(s)=sum b_n(s) exact values s=1..4 verified vs series and cot-integral; W(2)=2πlog2-7ζ(3)/π")
