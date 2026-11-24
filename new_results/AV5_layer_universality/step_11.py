import mpmath
from mpmath import mp, mpf, sqrt, pi, log, cos, sin, cosh, tanh, sech, matrix, lu_solve
mp.dps = 60

# numerical G(tau,s) = log(h/m) at diamond end, h=(c^q+s^q)^(1/q), m=max
def Gnum(tau, s):
    th = pi/4 + tau*s
    q = 1/tau
    h = (cos(th)**q + sin(th)**q)**(1/q)
    m = max(cos(th), sin(th))
    return log(h/m)

def g1f(s): return log(2*cosh(s)) - s
def g3f(s): return mpf(2)/3*s**3*(tanh(s)-1)
def g5f(s): return mpf(2)/3*s**5*(tanh(s)-1) + mpf(2)/9*s**6*sech(s)**2
def g7f(s): return (-mpf(8)/81*s**9*tanh(s)*sech(s)**2 + mpf(4)/9*s**8*sech(s)**2
                    + mpf(244)/315*s**7*(tanh(s)-1))

for s0 in [mpf('1.0'), mpf('0.7')]:
    taus = [mpf('0.004')/2**j for j in range(5)]
    Gs = [Gnum(t, s0) for t in taus]
    # G = g1 t + g3 t^3 + g5 t^5 + g7 t^7 + g9 t^9 ; divide by t: fit odd poly
    n=5
    A = matrix(n,n); b = matrix(n,1)
    for i,t in enumerate(taus):
        for j in range(n): A[i,j] = t**(2*j)
        b[i] = Gs[i]/t
    c = lu_solve(A,b)
    print(f"s={s0}: g1 {mp.nstr(c[0],20)} vs {mp.nstr(g1f(s0),20)}")
    print(f"        g3 {mp.nstr(c[1],20)} vs {mp.nstr(g3f(s0),20)}")
    print(f"        g5 {mp.nstr(c[2],20)} vs {mp.nstr(g5f(s0),20)}")
    print(f"        g7 {mp.nstr(c[3],20)} vs {mp.nstr(g7f(s0),20)}")
