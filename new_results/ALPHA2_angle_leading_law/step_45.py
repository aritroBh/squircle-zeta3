from mpmath import mp, mpf, nstr, pi, sin, cos, gamma, hyper, diffs, factorial, log
mp.dps = 60
def Eafun(s, ca, alpha):
    c2 = (1+ca)/2
    return gamma(1-s)*gamma(2+s)*hyper([1+s, -s],[2], c2)/(1-ca)
def ecoeffs(ca, alpha, nmax=9):
    f = lambda s: Eafun(s, ca, alpha)
    return [diffs(f, 0, n)/factorial(n) for n in range(nmax)]
# squircle check: e_w should be [y^w] Gamma(1+y/2)Gamma(1/2-y/2)/sqrt(pi)
from mpmath import sqrt
def Esq(s): return gamma(1+s/2)*gamma(mpf('0.5')-s/2)/sqrt(pi)
esq = [diffs(Esq, 0, n)/factorial(n) for n in range(9)]
ea = ecoeffs(mpf('0'), pi/2)
print("w : squircle-E Taylor vs Ea(0)-closed-form")
for w in range(9):
    print(w, nstr(esq[w], 15), nstr(ea[w], 15), nstr(abs(esq[w]-ea[w]),3))
# e_1 closed form check at several angles
for alpha in [pi/3, 2*pi/3, pi/6, pi/4]:
    ca = cos(alpha)
    e1 = ecoeffs(ca, alpha, 2)[1]
    e1f = -log(sin(alpha/2))/cos(alpha/2)**2
    print(f"alpha={nstr(alpha,6)}: e_1={nstr(e1,20)} closed={nstr(e1f,20)} |diff|={nstr(abs(e1-e1f),3)}")
