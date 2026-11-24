# LAY1: numerical sanity checks for the Layer Theorem (mpmath, dps 50).
import mpmath as mp, sys
mp.mp.dps = 50
z2, z3, z4 = mp.zeta(2), mp.zeta(3), mp.zeta(4)
s2 = mp.sqrt(2)

def pts(a, b, tau, at_left):
    # breakpoints geometric in units of tau, clustered at the layer end
    p = [mp.mpf(0)]
    k = mp.mpf(1)
    while k*tau < (b-a):
        p.append(k*tau); k *= 2
    p.append(b-a)
    return [a+x for x in p] if at_left else [b-x for x in reversed(p)]

def P_squircle(tau):
    q = 1/tau
    f = lambda t: mp.cos(t)*(1+mp.exp(q*mp.log(mp.tan(t))))**tau if t > 0 else mp.cos(t)
    I, err = mp.quad(f, pts(mp.mpf(0), mp.pi/4, tau, False), error=True)
    return 8*I, 8*err

def P_softmax(tau):
    def f(t):
        c, s = mp.cos(t), mp.sin(t)
        return (c+s)*(1+mp.exp(-2*s/tau)+mp.exp(-2*c/tau)+mp.exp(-2*(c+s)/tau))**tau
    I, err = mp.quad(f, pts(mp.mpf(0), mp.pi/4, tau, True), error=True)
    return 8*I, 8*err

cases = {
  'squircle': (P_squircle, 4*s2, s2*z2, s2*z3, -33*s2*z4/16),
  'softmax' : (P_softmax, mp.mpf(8), mp.pi**2/3, 2*z3, z4/4),
}
taus = [mp.mpf(1)/n for n in (10, 20, 40, 80, 160, 320, 640, 1280, 2560)]
for name, (P, c0, c2, c3, c4) in cases.items():
    print(f"== {name}: c0={mp.nstr(c0,20)} c2={mp.nstr(c2,20)} c3={mp.nstr(c3,20)} predicted c4={mp.nstr(c4,20)}")
    print(" 1/tau | quad err | R3/tau^4 = (P-c0-c2 t^2-c3 t^3)/t^4 | R4/tau^5 = (P-...-c4 t^4)/t^5")
    for t in taus:
        v, e = P(t)
        R3 = (v - c0 - c2*t**2 - c3*t**3)/t**4
        R4 = (v - c0 - c2*t**2 - c3*t**3 - c4*t**4)/t**5
        print(f" {int(1/t):5d} | {mp.nstr(e,3):8s} | {mp.nstr(R3,25):32s} | {mp.nstr(R4,15)}")
        sys.stdout.flush()
