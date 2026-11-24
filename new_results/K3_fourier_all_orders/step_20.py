import mpmath as mp
mp.mp.dps = 50
pi = mp.pi; s2 = mp.sqrt(2); z = mp.zeta

phi  = lambda t: mp.exp(mp.sin(t) + mp.cos(3*t)/2)
def d2phi(t): return mp.diff(phi, t, 2)
def d4phi(t): return mp.diff(phi, t, 4)
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

Dt = {
 (2,0): [z(2)], (3,0): [z(3)],
 (4,0): [-sp_rat if False else -mp.mpf(33)/16*z(4)], (4,1): [mp.mpf(7)/2*z(4)],
 (5,0): [-mp.mpf(3)/4*z(5)-mp.mpf(7)/4*z(2)*z(3)], (5,1): [2*z(5)+2*z(2)*z(3)],
 (6,0): [mp.mpf(4039)/256*z(6)-mp.mpf(7)/8*z(3)**2], (6,1): [-mp.mpf(733)/16*z(6)+z(3)**2], (6,2): [mp.mpf(93)/2*z(6)],
}
def Dseries(tau, m):
    return s2/pi*mp.fsum(Dt[(n,m)][0]*tau**n for n in range(2*m+2, 7))

S0 = mp.fsum(phi(pi/4 + j*pi/2) for j in range(4))
S2 = mp.fsum(d2phi(pi/4 + j*pi/2) for j in range(4))
S4 = mp.fsum(d4phi(pi/4 + j*pi/2) for j in range(4))
print("S0 =", mp.nstr(S0, 10), " S2 =", mp.nstr(S2, 10), " S4 =", mp.nstr(S4, 10))

def comb(tau, order):  # include terms with tau^n, n <= order
    tot = mp.mpf(0)
    for m in range(3):
        if 2*m+2 <= order:
            S = [S0, S2, S4][m]
            tot += (pi/4) * mp.mpf(2)**(-2*m)/mp.factorial(2*m) * S * Dseries(tau, m)
    return tot

print("\nresiduals after subtracting comb through order N:")
for tau in [mp.mpf('0.04'), mp.mpf('0.02'), mp.mpf('0.01'), mp.mpf('0.005')]:
    L = lhs(tau)
    r3 = (L - comb(tau, 3))/tau**4
    r4 = (L - comb(tau, 4))/tau**5
    r5 = (L - comb(tau, 5))/tau**6
    r6 = (L - comb(tau, 6))/tau**7
    print(f"tau={mp.nstr(tau,4)}: /t^4 {mp.nstr(r3,8)} | /t^5 {mp.nstr(r4,8)} | /t^6 {mp.nstr(r5,8)} | /t^7 {mp.nstr(r6,8)}")

# predicted tau^7 term: m=0 D7,0 + m=1 D7,1 (+m=2 D7,2? no, that starts at tau^6... D_{7,2}? m=2 starts n=6; n=7 has m=0,1,2? r<= (n-2)/2=2.5 -> m=0,1,2... wait for n=7: m up to (7-2)//2 = 2)
# predicted residual/tau^7 limit = (pi/4)[ S0*(sqrt2/pi)*Dt70 + (1/8) S2*(sqrt2/pi)*Dt71 + (1/192) S4*(sqrt2/pi)*Dt72 ]
Dt70 = mp.mpf(129)/16*z(7)+mp.mpf(141)/16*z(2)*z(5)+mp.mpf(1071)/160*z(2)**2*z(3)
Dt71 = -23*z(7)-25*z(2)*z(5)-mp.mpf(181)/10*z(2)**2*z(3)
Dt72 = 24*z(7)+24*z(2)*z(5)+mp.mpf(84)/5*z(2)**2*z(3)
pred7 = (pi/4)*(s2/pi)*(S0*Dt70 + S2*Dt71/8 + S4*Dt72/192)
print("\npredicted limit of residual/tau^7 =", mp.nstr(pred7, 10))
