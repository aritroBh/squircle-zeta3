# Referee item 5: independent end-to-end I(tau) and R3/tau^4 with Richardson (own code).
import sys
from mpmath import mp, mpf, quad, exp, log, log1p, expm1, tan, cos, sin, sqrt, pi, zeta, nstr, quadgl
DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 60
mp.dps = DPS
def bps(tau):
    ks = [0,0.25,0.5,1,1.5,2,3,4,6,8,12,16,24,32,48,64,96,128]
    pts = [tau*mpf(k) for k in ks if tau*k < pi/4]
    return pts + [pi/4]
def I_squircle(tau, method='tanh-sinh'):
    # theta = pi/4 - x, x in [0, pi/4]; m = (cos x + sin x)/sqrt2 ; r = tan(pi/4 - x)
    def d(x):
        m = (cos(x)+sin(x))/sqrt(2)
        if x == pi/4: return mpf(0)
        r = tan(pi/4 - x)
        return m*expm1(tau*log1p(exp(log(r)/tau)))
    v, e = quad(d, bps(tau), error=True, method=method)
    return 4*sqrt(2) + 8*v, 8*e
def I_softmax(tau, method='tanh-sinh'):
    def d(x):
        c, s_ = cos(x), sin(x)
        m = c + s_
        L = tau*log1p(exp(-2*s_/tau) + exp(-2*c/tau) + exp(-2*m/tau))
        return m*expm1(L)
    v, e = quad(d, bps(tau), error=True, method=method)
    return 8 + 8*v, 8*e
z = zeta
cases = {
 'squircle': (I_squircle, 4*sqrt(2), sqrt(2)*z(2), sqrt(2)*z(3), -33*sqrt(2)/16*z(4)),
 'softmax' : (I_softmax, mpf(8), pi**2/3, 2*z(3), z(4)/4),
}
invtaus = [20,40,80,160,320,640,1280,2560]
for name,(fn,c0,c2,c3,c4) in cases.items():
    print(f"=== {name}  (dps {DPS})  claimed c4 = {nstr(c4,25)}")
    taus=[]; R=[]
    for it in invtaus:
        tau = mpf(1)/it
        I, err = fn(tau)
        Igl, _ = fn(tau, 'gauss-legendre')
        r3 = (I - c0 - c2*tau**2 - c3*tau**3)/tau**4
        taus.append(tau); R.append(r3)
        print(f"1/tau={it:5d}  I={nstr(I,50)}  quad_err_est={nstr(err,2)}  |TS-GL|={nstr(abs(I-Igl),2)}")
        print(f"            R3/tau^4={nstr(r3,30)}  (R3/tau^4-c4)/tau={nstr((r3-c4)/tau,15)}")
    # Richardson (Neville to tau=0, polynomial in tau)
    def neville(xs, ys):
        P = list(ys); n=len(xs)
        for k in range(1,n):
            for i in range(n-k):
                P[i] = (xs[i+k]*P[i] - xs[i]*P[i+1])/(xs[i+k]-xs[i])  # extrapolate to 0
        return P[0]
    req = [invtaus.index(k) for k in (20,80,320,1280)]
    for label, idx in (("required 4 pts 1/tau=20,80,320,1280", req),
                       ("last 4 pts 1/tau=320..2560", list(range(4,8))),
                       ("last 5 pts 1/tau=160..2560", list(range(3,8))),
                       ("last 6 pts 1/tau=80..2560", list(range(2,8))),
                       ("all 8 pts", list(range(8)))):
        ext = neville([taus[i] for i in idx],[R[i] for i in idx])
        print(f"  Richardson [{label}]: c4 ~ {nstr(ext,25)}  ext-c4 = {nstr(ext-c4,3)}")
    # c5 estimate
    S5 = [(R[i]-c4)/taus[i] for i in range(8)]
    ext5 = neville(taus[3:], S5[3:])
    print(f"  c5 estimate (Richardson on (R3/tau^4-c4)/tau, last 5): {nstr(ext5,15)}")
