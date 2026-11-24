# Referee item 1: S(p) from the definition vs 8(1-tau)^2 T(tau).
# Route A (definition, radial chart): x = r(u) u on the unit sphere, F = sum |x_i|^p,
#   dA = r^3 |grad F(x)| / p dsigma   (since u.n = p/(r |grad F|) for F homogeneous of degree p)
#   |grad F|/p = (sum x_i^{2(p-1)})^{1/2}.  8 octants.
# Route B (proof's start): T = int_Delta (t1 t2 t3)^{-tau} (sum t_i^{2 tau})^{1/2} dt1 dt2,
#   derived by z=(1-x^p-y^p)^{1/p}, dA = |grad F|/|F_z| dx dy, t_i = x_i^p.
import mpmath as mp
mp.mp.dps = 20

def S_radial(p):
    p = mp.mpf(p)
    def integrand(th, ph):
        u = [mp.sin(th)*mp.cos(ph), mp.sin(th)*mp.sin(ph), mp.cos(th)]
        r = (sum(ui**p for ui in u))**(-1/p)
        x = [r*ui for ui in u]
        g = mp.sqrt(sum(xi**(2*(p-1)) for xi in x))
        return r**3*g*mp.sin(th)
    return 8*mp.quad(integrand, [0, mp.pi/2], [0, mp.pi/2])

def S_T(tau):
    tau = mp.mpf(tau)
    def f(t1, t2):
        t3 = 1-t1-t2
        if t3 <= 0 or t1 <= 0 or t2 <= 0:
            return mp.mpf(0)
        return (t1*t2*t3)**(-tau)*mp.sqrt(t1**(2*tau)+t2**(2*tau)+t3**(2*tau))
    T = mp.quad(lambda t1: mp.quad(lambda t2: f(t1, t2), [0, 1-t1]), [0, 1])
    return 8*(1-tau)**2*T, T

if __name__ == "__main__":
    for tau in ['0.5', '0.2', '0.02']:
        tau = mp.mpf(tau); p = 1/(1-tau)
        A = S_radial(p)
        B, T = S_T(tau)
        print(f"tau={mp.nstr(tau,4)} p={mp.nstr(p,8)}  S_radial={mp.nstr(A,16)}  8(1-tau)^2T={mp.nstr(B,16)}  rel.diff={mp.nstr(abs(A-B)/A,3)}")
    print("sphere check 4*pi =", mp.nstr(4*mp.pi, 16))
    print("S(1)=8*T(0)=8*sqrt(3)/2=4sqrt3 =", mp.nstr(4*mp.sqrt(3), 16))
