from mpmath import mp, mpf, quad, sqrt, nstr, pi, sin, cos, gamma, legenp
mp.dps = 50
def Qf(w, ca): return sqrt(1 - 2*w*ca + w**2)
def Eint(s, ca): return quad(lambda w: (w**(1+s) + w**(-s))/Qf(w,ca)**3, [0, 1])
def Ecand(s, ca, alpha):
    return gamma(1-s)*gamma(2+s)*legenp(-1-s, -1, -ca, type=2)/sin(alpha)
for alpha in [pi/2, pi/3, 2*pi/3, pi/6, 5*pi/6, pi/4]:
    ca = cos(alpha)
    for s in [mpf('0.3'), mpf('-0.2'), mpf('0.55'), mpf('-0.75')]:
        a = Eint(s, ca); b = Ecand(s, ca, alpha)
        print(f"alpha={nstr(alpha,6)}, s={nstr(s,5)}: |E-cand|={nstr(abs(a-b),3)}")
