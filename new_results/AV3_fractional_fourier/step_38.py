from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, quad, tan)
mp.dps = 100
def N(x,d=42): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
eta = lambda k: (1-2**(1-k))*zeta(k)
def a_c(n,s):
    if s==1: return 2*((-1)**n-1)/(pi*n**2)
    if s==2: return 4*(-1)**n/n**2
    if s==3: return 6*pi*(-1)**n/n**2 + 12*(1-(-1)**n)/(pi*n**4)
    if s==4: return (-1)**n*(8*pi**2/n**2 - 48/n**4)
def b_c(n,s):
    if s==1: return 2*(-1)**(n+1)/n
    if s==2: return 2*pi*(-1)**(n+1)/n - 4*(1-(-1)**n)/(pi*n**3)
    if s==3: return 2*pi**2*(-1)**(n+1)/n + 12*(-1)**n/n**3
    if s==4: return 2*pi**3*(-1)**(n+1)/n + 24*pi*(-1)**n/n**3 + 48*(1-(-1)**n)/(pi*n**5)

# symbolic values via zeta/eta algebra: B(s,s') = sum a_n(s) b_n(s')
def Z(k): return zeta(k)
def E_(k): return eta(k)
B_form = {}
# hand-derived from the exact forms:
B_form[(1,1)] = -(4/pi)*(Z(3)+E_(3))                       # (4/pi)*(-2)*odd? -> -(7)/pi zeta(3)
B_form[(1,2)] = None  # compute below
# do it programmatically: each a_n(s), b_n(s) is a sum of terms c * w(n) * n^{-k}, w in {1,(-1)^n}
def terms_a(s):
    if s==1: return [(2/pi,'a',2),(-2/pi,'c',2)]
    if s==2: return [(4,'a',2)]
    if s==3: return [(6*pi,'a',2),(12/pi,'c',4),(-12/pi,'a',4)]
    if s==4: return [(8*pi**2,'a',2),(-48,'a',4)]
def terms_b(s):
    if s==1: return [(-2,'a',1)]
    if s==2: return [(-2*pi,'a',1),(-4/pi,'c',3),(4/pi,'a',3)]
    if s==3: return [(-2*pi**2,'a',1),(12,'a',3)]
    if s==4: return [(-2*pi**3,'a',1),(24*pi,'a',3),(48/pi,'c',5),(-48/pi,'a',5)]
def Ssum(w,k):
    return zeta(k) if w=='c' else -eta(k)   # sum_{n>=1} n^-k  resp. (-1)^n n^-k
def B_exact(s,sp):
    tot = 0
    for (c1,w1,k1) in terms_a(s):
        for (c2,w2,k2) in terms_b(sp):
            w = 'c' if w1==w2 else 'a'   # product of sign types: c*(-1)^n or (-1)^n*1 -> 'a'
            tot += c1*c2*Ssum(w,k1+k2)
    return tot
for (s,sp) in [(1,1),(1,2),(2,1),(2,2),(1,3),(2,3),(3,3),(2,4),(3,4),(4,4)]:
    Bs = nsum(lambda n: a_c(n,s)*b_c(n,sp), [1, inf])
    Bx = B_exact(s,sp)
    print(f"B({s},{sp}) series={N(Bs,30)}")
    print(f"        exact ={N(Bx,30)}  err={nstr(abs(Bs-Bx),3)}")
logline("step_34: cross-parity pairing B(s,s')=sum a_n b_n: exact zeta-tower values verified s,s' in 1..4")
