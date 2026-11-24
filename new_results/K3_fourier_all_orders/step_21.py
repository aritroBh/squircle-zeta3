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
 (2,0): z(2), (3,0): z(3),
 (4,0): -mp.mpf(33)/16*z(4), (4,1): mp.mpf(7)/2*z(4),
 (5,0): -mp.mpf(3)/4*z(5)-mp.mpf(7)/4*z(2)*z(3), (5,1): 2*z(5)+2*z(2)*z(3),
 (6,0): mp.mpf(4039)/256*z(6)-mp.mpf(7)/8*z(3)**2, (6,1): -mp.mpf(733)/16*z(6)+z(3)**2, (6,2): mp.mpf(93)/2*z(6),
}
S0 = mp.fsum(phi(pi/4 + j*pi/2) for j in range(4))
S2 = mp.fsum(d2phi(pi/4 + j*pi/2) for j in range(4))
S4 = mp.fsum(d4phi(pi/4 + j*pi/2) for j in range(4))
Ss = {0: S0, 1: S2, 2: S4}

def comb(tau, order):
    tot = mp.mpf(0)
    for m in range(3):
        for n in range(2*m+2, order+1):
            tot += (pi/4)*(s2/pi)*mp.mpf(2)**(-2*m)/mp.factorial(2*m)*Ss[m]*Dt[(n,m)]*tau**n
    return tot

print("residual/tau^(N+1) after subtracting comb through order N:")
for tau in [mp.mpf('0.04'), mp.mpf('0.02'), mp.mpf('0.01'), mp.mpf('0.005')]:
    L = lhs(tau)
    vals = []
    for N in [3, 4, 5, 6]:
        vals.append(mp.nstr((L - comb(tau, N))/tau**(N+1), 8))
    print(f"tau={mp.nstr(tau,4)}: N=3..6: {vals}")

# predicted limits
pred = {}
for N in [3,4,5,6]:
    n = N+1
    tot = mp.mpf(0)
    for m in range((n-2)//2+1):
        pass
print("\npredicted tau^4 limit = (√2 ζ4/64)(-33 S0 + 7 S2) =", mp.nstr(s2*z(4)/64*(-33*S0+7*S2), 10))
