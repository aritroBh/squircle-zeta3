from mpmath import mp, mpf, quad, sqrt, nstr, pi, sin
mp.dps = 50
def Qf(w, ca): return sqrt(1 - 2*w*ca + w**2)
def Ea(s, ca): return quad(lambda w: (w**(1+s) + w**(-s))/Qf(w,ca)**3, [0, 1])
def WW(s, ca): return quad(lambda w: w**(-s)*(1+w)/Qf(w,ca)**3, [0, 1]) + quad(lambda w: w**(-s)*(1+w)/Qf(w,ca)**3, [1, mpf('inf')])
def GG(s, ca):
    g1 = lambda w: (1-w)/(1+w)*((1+w)/Qf(w,ca) - 1)/w
    return quad(lambda w: w**(-s)*g1(w), [0, 1]) + quad(lambda w: w**(-s)*g1(w), [1, mpf('inf')])
def e(x): return nstr(abs(x), 3)
for ca in [mpf('0'), mpf('0.28'), mpf('-0.5'), mpf('0.96'), mpf('-0.875')]:
    for s in [mpf('0.3'), mpf('-0.2'), mpf('0.7'), mpf('-0.6')]:
        lhs = GG(s,ca); rhs = 2*pi/sin(pi*s) - (1-ca)*WW(s,ca)/s
        print(f"ca={nstr(ca,6)}, s={nstr(s,4)}: |GG - [2pi/sin - (1-ca)WW/s]|={e(lhs-rhs)}")
