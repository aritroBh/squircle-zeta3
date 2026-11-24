from mpmath import mp, mpf, quad, nstr, pi, sin, cos, legenp
mp.dps = 55
def Jfull_trig(s, alpha):  # int_0^{pi-alpha} sin^{s-1}u sin^{-s}(alpha+u) du
    return quad(lambda u: sin(u)**(s-1)*sin(alpha+u)**(-s), [0, pi-alpha])
def Jcand(s, ca): return pi*legenp(s-1, 0, -ca, type=2)/sin(pi*s)
worst = 0
for alpha in [pi/2, pi/3, 2*pi/3, pi/6, 3*pi/4, pi/5]:
    ca = cos(alpha)
    for s in [mpf('0.3'), mpf('0.5'), mpf('0.7'), mpf('0.25'), mpf('0.9')]:
        a = Jfull_trig(s, alpha); b = Jcand(s, ca)
        d = abs(a-b)
        worst = max(worst, d)
        print(f"alpha={nstr(alpha,6)}, s={nstr(s,4)}: |diff|={nstr(d,3)}")
print("WORST:", nstr(worst,3))
