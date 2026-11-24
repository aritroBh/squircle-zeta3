from mpmath import mp, mpf, quad, sqrt, nstr, pi, sin, cos, beta, legenp
mp.dps = 50
def Qf(w, ca): return sqrt(1 - 2*w*ca + w**2)
def Iint(mu, ca):
    return quad(lambda w: w**(mu-1)/Qf(w,ca)**3, [0, 1]) + quad(lambda w: w**(mu-1)/Qf(w,ca)**3, [1, mpf('inf')])
def Icand(mu, ca, alpha):
    return (2/sin(alpha))*beta(mu, 3-mu)*legenp(mu-2, -1, -ca, type=2)
for alpha in [pi/2, pi/3, 2*pi/3, pi/6, 5*pi/6]:
    ca = cos(alpha)
    for mu in [mpf('0.5'), mpf('1.3'), mpf('2.2')]:
        a = Iint(mu, ca); b = Icand(mu, ca, alpha)
        print(f"alpha={nstr(alpha,6)}, mu={nstr(mu,4)}: |I-cand|={nstr(abs(a-b),3)}")
