from mpmath import mp, mpf, nstr, pi, acos, gamma, hyper, diffs, factorial, polylog, zeta, log
mp.dps = 60
def Jc(c):
    c = mpf(c)
    return log(1-c**2)*(polylog(2,c**2)-zeta(2)) + 2*polylog(3,1-c**2) - 2*zeta(3) + log(c**2)*log(1-c**2)**2
def e3(ca):
    def Eafun(s):
        return gamma(1-s)*gamma(2+s)*hyper([1+s, -s],[2], (1+ca)/2)/(1-ca)
    return list(diffs(Eafun, 0, 3))[-1]/factorial(3)
for c in ['1/2','3/5','1/3','4/5','2/3','5/7','1/4','6/7']:
    cm = mpf(c); ca = 2*cm**2-1
    lhs = e3(ca); rhs = -Jc(cm)/(2*cm**2)
    print(f"c={c}: e_3={nstr(lhs,25)}  -J/(2c^2)={nstr(rhs,25)}  |diff|={nstr(abs(lhs-rhs),3)}")
