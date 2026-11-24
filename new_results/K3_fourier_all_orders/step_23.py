import mpmath as mp
mp.mp.dps = 50
pi = mp.pi; s2 = mp.sqrt(2); z = mp.zeta

phi  = lambda t: mp.exp(mp.sin(t) + mp.cos(3*t)/2)
d2phi = lambda t: mp.diff(phi, t, 2)
d4phi = lambda t: mp.diff(phi, t, 4)
def h(tau, t):
    q = 1/tau
    return (abs(mp.cos(t))**q + abs(mp.sin(t))**q)**tau
def h1(t): return max(abs(mp.cos(t)), abs(mp.sin(t)))
def lhs(tau):
    tot = mp.mpf(0)
    for j in range(4):
        tj = pi/4 + j*pi/2
        pts = [tj-pi/4, tj-8*tau, tj-2*tau, tj, tj+2*tau, tj+8*tau, tj+pi/4]
        tot += mp.quad(lambda t: (h(tau,t)-h1(t))*phi(t), pts)
    return tot
S0 = mp.fsum(phi(pi/4 + j*pi/2) for j in range(4))
S2 = mp.fsum(d2phi(pi/4 + j*pi/2) for j in range(4))
S4 = mp.fsum(d4phi(pi/4 + j*pi/2) for j in range(4))
Ss = {0: S0, 1: S2, 2: S4}
Dt = {
 (2,0): z(2), (3,0): z(3),
 (4,0): -mp.mpf(33)/16*z(4), (4,1): mp.mpf(7)/2*z(4),
 (5,0): -mp.mpf(3)/4*z(5)-mp.mpf(7)/4*z(2)*z(3), (5,1): 2*z(5)+2*z(2)*z(3),
 (6,0): mp.mpf(4039)/256*z(6)-mp.mpf(7)/8*z(3)**2, (6,1): -mp.mpf(733)/16*z(6)+z(3)**2, (6,2): mp.mpf(93)/2*z(6),
 (7,0): mp.mpf(129)/16*z(7)+mp.mpf(141)/16*z(2)*z(5)+mp.mpf(1071)/160*z(2)**2*z(3),
 (7,1): -23*z(7)-25*z(2)*z(5)-mp.mpf(181)/10*z(2)**2*z(3),
 (7,2): 24*z(7)+24*z(2)*z(5)+mp.mpf(84)/5*z(2)**2*z(3),
}
def comb(tau, order):
    return mp.fsum((pi/4)*(s2/pi)*mp.mpf(2)**(-2*m)/mp.factorial(2*m)*Ss[m]*Dt[(n,m)]*tau**n
                   for m in range(4) for n in range(2*m+2, order+1))
def predlim(N):
    n = N+1
    return (pi/4)*(s2/pi)*mp.fsum(mp.mpf(2)**(-2*m)/mp.factorial(2*m)*Ss[m]*Dt[(n,m)] for m in range((n-2)//2+1))

taus = [mp.mpf('0.02'), mp.mpf('0.01'), mp.mpf('0.005')]
print("N: residuals(.02,.01,.005), 2-step Richardson limit, predicted limit, |diff|")
for N in [3, 4, 5, 6]:
    rs = [(lhs(t) - comb(t, N))/t**(N+1) for t in taus]
    L1 = 2*rs[1]-rs[0]; L2 = 2*rs[2]-rs[1]; L3 = (4*L2-L1)/3
    print(f"N={N}: {[mp.nstr(r,9) for r in rs]} extrap={mp.nstr(L3,10)} pred={mp.nstr(predlim(N),10)} diff={mp.nstr(abs(L3-predlim(N)),3)}")
