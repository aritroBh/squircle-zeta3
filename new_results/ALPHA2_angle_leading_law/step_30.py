from mpmath import mp, mpf, quad, sqrt, nstr
mp.dps = 50
def Qf(w, ca): return sqrt(1 - 2*w*ca + w**2)
def Ea(s, ca): return quad(lambda w: (w**(1+s) + w**(-s))/Qf(w,ca)**3, [0, 1])
def WW(s, ca): return quad(lambda w: w**(-s)*(1+w)/Qf(w,ca)**3, [0, 1]) + quad(lambda w: w**(-s)*(1+w)/Qf(w,ca)**3, [1, mpf('inf')])
def JJ(s, ca): return quad(lambda w: (w**(s-1) + w**(-s))/Qf(w,ca), [0, 1])  # symmetrized, -1<Re s<1
def e(x): return nstr(abs(x), 3)
for ca in [mpf('0'), mpf('0.28'), mpf('-0.5'), mpf('0.96')]:
    for s in [mpf('0.3'), mpf('-0.2'), mpf('0.7')]:
        ww = WW(s,ca)
        jj = JJ(s,ca)-JJ(-s,ca); rhs3 = (1-ca)*ww/s
        print(f"ca={nstr(ca,5)}, s={nstr(s,4)}: |dJ-(1-ca)WW/s|={e(jj-rhs3)}")
