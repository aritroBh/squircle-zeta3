from mpmath import mp, mpf, quad, sqrt, nstr, pi, sin, cos, legenp
mp.dps = 50
def Qf(w, ca): return sqrt(1 - 2*w*ca + w**2)
def Jfull(s, ca):  # int_0^inf w^{s-1} Q^{-1} dw, 0<Re s<1
    return quad(lambda w: w**(s-1)/Qf(w,ca), [0, 1]) + quad(lambda w: w**(s-1)/Qf(w,ca), [1, mpf('inf')])
def Jcand(s, ca): return pi*legenp(s-1, 0, -ca, type=2)/sin(pi*s)
for alpha in [pi/2, pi/3, 2*pi/3, pi/6, 3*pi/4]:
    ca = cos(alpha)
    for s in [mpf('0.3'), mpf('0.5'), mpf('0.7'), mpf('0.25')]:
        a = Jfull(s, ca); b = Jcand(s, ca)
        print(f"alpha={nstr(alpha,6)}, s={nstr(s,4)}: |Jfull - pi P_{{s-1}}(-ca)/sin|={nstr(abs(a-b),3)}")
