from mpmath import mp, mpf, nstr, pi, sin, cos, gamma, hyper, diffs, factorial, log, sqrt
mp.dps = 60
def dn(f, x, n): return list(diffs(f, x, n))[-1]/factorial(n)
def Eafun(s, ca, alpha):
    c2 = (1+ca)/2
    return gamma(1-s)*gamma(2+s)*hyper([1+s, -s],[2], c2)/(1-ca)
def Esq(s): return gamma(1+s/2)*gamma(mpf('0.5')-s/2)/sqrt(pi)
esq = [dn(Esq, 0, n) for n in range(9)]
ea = [dn(lambda s: Eafun(s, mpf('0'), pi/2), 0, n) for n in range(9)]
print("w : squircle-E  vs  Ea-closed-form  |diff|")
for w in range(9):
    print(w, nstr(esq[w], 15), nstr(ea[w], 15), nstr(abs(esq[w]-ea[w]),3))
for alpha in [pi/3, 2*pi/3, pi/6, pi/4]:
    ca = cos(alpha)
    e1 = dn(lambda s: Eafun(s, ca, alpha), 0, 1)
    e1f = -log(sin(alpha/2))/cos(alpha/2)**2
    print(f"alpha={nstr(alpha,6)}: e_1={nstr(e1,20)} closed={nstr(e1f,20)} |diff|={nstr(abs(e1-e1f),3)}")
