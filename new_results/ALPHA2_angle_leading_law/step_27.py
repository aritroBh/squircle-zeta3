from mpmath import mp, mpf, quad, sqrt, cos, pi
mp.mp.dps = 50

def Qf(w, ca): return sqrt(1 - 2*w*ca + w**2)
def Ea(s, ca):
    return quad(lambda w: (w**(1+s) + w**(-s))/Qf(w,ca)**3, [0, 1])
def WW(s, ca):
    return quad(lambda w: w**(-s)*(1+w)/Qf(w,ca)**3, [0, 1]) + quad(lambda w: w**(-s)*(1+w)/Qf(w,ca)**3, [1, mpf('inf')])
def FF(s, ca):
    h = lambda w: ((1+w)/Qf(w,ca) - 1)/w
    return quad(lambda w: w**(-s)*h(w), [0, 1]) + quad(lambda w: w**(-s)*h(w), [1, mpf('inf')])
def JJ(s, ca):
    return quad(lambda w: w**(s-1)/Qf(w,ca), [0, 1]) + quad(lambda w: w**(s-1)/Qf(w,ca), [1, mpf('inf')])

for ca in [mpf('0'), mpf('0.28'), mpf('-0.5')]:
    for s in [mpf('0.3'), mpf('-0.2')]:
        ww = WW(s,ca); ee = Ea(s,ca)+Ea(-s,ca); ff = FF(s,ca); rhs2 = (1+ca)*(Ea(s,ca)-Ea(-s,ca))/s
        jj = JJ(s,ca)-JJ(-s,ca); rhs3 = (1-ca)*ww/s
        print(f"ca={ca}, s={s}: |WW-(E+E-)|={abs(ww-ee):.2e}  |FF-(1+ca)(E-E-)/s|={abs(ff-rhs2):.2e}  |dJ-(1-ca)WW/s|={abs(jj-rhs3):.2e}")
