# Planner: transcribe the real-Li2 theorem from verify_proof.md by hand and compare with direct quadrature.
import mpmath as mp
mp.mp.dps = 40
L2 = lambda x: mp.re(mp.polylog(2, x))
lg = mp.log
def M1_thm(a):
    c, s = mp.cos(a/2), mp.sin(a/2)
    Xm, Xp = (1+s-c)/(1-c), (1+s+c)/(1+c)
    U = (-L2(-1/s) + L2(-s/(1-c)) + L2(-s/(1+c)) - L2(-(1+s)/c) - L2(-s*(1+s)/c**2) + L2(-s*(1+s)/(c*(1+c)))
         + c*L2(-(1+s-c)/c) + c*L2(-(1+s-c)/(c*(1-c))) - 2*c*L2(-(1+s-c)/(2*c))
         - 2*L2(1/(1+s)) + L2(c**2/(1+s)) + L2(c/(1+s)) + L2(s/(1+s)) - L2(c*(1-c)/(s*(1+s)))
         - c*L2((1-c)/(1+s-c)) + c*L2((1+c)/(1+s+c)) + c*L2(c/(1+s+c)) + c*L2(c*(1+c)/(1+s+c)) - 2*c*L2(2*c/(1+s+c)))
    V = ((1-c)*L2(-1/c) + L2(1-1/c**2) - (1+c)*L2(-(1-c)/c) + 2*c*L2(-(1-c)/(2*c)) - L2(c**2) - (1+c)*L2(c)
         + (1-c)*L2(c/(1+c)) + 2*c*L2(2*c/(1+c)))
    Lp = (lg(1+s)*lg(c/(1+s)) + lg(c/(1+c))*lg((1+s)/s) + c*lg(4)*lg(Xm) + c/2*lg(Xm)**2
          - c/2*lg(Xp)*lg((1+c)*(1+s+c)))
    return mp.zeta(2) + U + V + Lp
def M1_quad(a):
    ca = mp.cos(a)
    def f(g):
        S0 = 1 - 2*ca*g + g*g; f0 = mp.sqrt(S0)/((1+g)*g); l = mp.log(g)
        return f0*(-mp.log(1+g) + l + (-l*g*g + ca*l*g)/S0) - l/g
    return mp.quad(f, [0, mp.mpf(1)/4, mp.mpf(1)/2, 1])
for a in [mp.mpf('0.05'), mp.pi/5, mp.pi/2, 2*mp.pi/3, mp.mpf('2.5'), mp.mpf('3.1')]:
    t, q = M1_thm(a), M1_quad(a)
    print(mp.nstr(a, 6), mp.nstr(t, 22), "diff", mp.nstr(abs(t-q), 3))
