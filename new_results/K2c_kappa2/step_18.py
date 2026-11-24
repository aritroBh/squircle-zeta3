
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi, gamma
mp.dps = 70

def saddle_xs(N, Z):
    return (-Z + sqrt(Z**2 + 4*N*Z))/2

def Jm_pieces(N, Z, m):
    x1 = saddle_xs(N, Z)
    d2 = N*(1/x1**2 - 1/(x1+Z)**2); ad2 = abs(d2)
    v = exp(-mpc(0,1)*mp.arg(d2)/2)
    w = 1/sqrt(ad2)
    S, T = 12.0*w, 12.0*w
    xa = x1 - S*v; xb = x1 + T*v
    def f(x):
        return exp(-x)*(1+Z/x)**(-N)*(x+Z)**(m+mpf(3)/2)/x
    IA = xa*mp.quadgl(lambda t: f(t*xa), [0, mpf('0.3'), mpf('0.6'), mpf('0.85'), 1])
    IB = v*mp.quadgl(lambda s: f(x1+s*v), [-S, 0, T])
    IC = mp.quadgl(lambda s: f(xb+s), [0, mpf('10'), mpf('40'), mp.inf])
    phis = x1 + N*log(1+Z/x1)
    model = exp(-phis)*sqrt(2*pi/d2)*(x1+Z)**(m+mpf(3)/2)/x1
    return IA, IB, IC, model, xa, xb, S, T

a = mp.pi/2; Z = 1j*a
for N in [1500, 4600, 10000]:
    IA, IB, IC, model, xa, xb, S, T = Jm_pieces(N, Z, 0)
    print('N=%d |x1|=%.1f S=T=%.1f |xa|=%.1f arg(xa)=%.3f' % (N, float(abs(saddle_xs(N,Z))), float(S), float(abs(xa)), float(mp.arg(xa))))
    print('  IA/model =', mp.nstr(IA/model, 8))
    print('  IB/model =', mp.nstr(IB/model, 8))
    print('  IC/model =', mp.nstr(IC/model, 8))
    print('  sum/model=', mp.nstr((IA+IB+IC)/model, 8))
