import mpmath as mp
mp.mp.dps = 90
sq2 = mp.sqrt(2)
x1 = 1-sq2; x2 = 1-1/sq2
kap1, kap2 = -x1, -x2
rho = x1/x2
N = 300
def Linst(s1,s2,a,b):
    tot = mp.mpf(0)
    for n in range(2,N):
        an = a**n/n**s1
        for m in range(1,n):
            tot += an*(b**m)/m**s2
    return tot
l1 = -mp.log(1-x1)
Li21 = mp.polylog(2,x1); Li22 = mp.polylog(2,x2); Li2p = mp.polylog(2,x1*x2)
L31 = mp.polylog(3,x1); L32 = mp.polylog(3,x2); L3p = mp.polylog(3,x1*x2)
def A(x): return mp.nsum(lambda n: mp.harmonic(n)*x**n/n**2, [1,mp.inf], tol=mp.mpf('1e-100'))
A1, A2 = A(x1), A(x2)
Li21a = Linst(2,1,x2,rho); Li21b = Linst(2,1,x1,1/rho)
Li12a = Linst(1,2,x2,rho); Li12b = Linst(1,2,x1,1/rho)
def H1(a,b):
    tot = mp.mpf(0)
    for n in range(2,N):
        an = a**n/n
        for m in range(1,n):
            tot += an*(b**m)*mp.harmonic(m)/m
    return tot
H1a = H1(x2,rho); H1b = H1(x1,1/rho)
I21 = mp.quad(lambda t: mp.log(1+kap1*t)**2*mp.log(1+kap2*t), [0,1])
I12 = mp.quad(lambda t: mp.log(1+kap1*t)*mp.log(1+kap2*t)**2, [0,1])
names = ['Li21(x2,r)','Li21(x1,1/r)','Li12(x2,r)','Li12(x1,1/r)','H1(x2,r)','H1(x1,1/r)',
         'A(x1)','A(x2)','Li3(x1)','Li3(x2)','Li3(x1x2)','l1^3','l1 Li2(x1)','l1 Li2(x2)','l1 Li2(x1x2)','1']
basis = [Li21a,Li21b,Li12a,Li12b,H1a,H1b,A1,A2,L31,L32,L3p,l1**3,l1*Li21,l1*Li22,l1*Li2p,mp.mpf(1)]
print("l1 =", mp.nstr(l1,30))
for tag, I in [('I21',I21),('I12',I12)]:
    rel = mp.pslq([I]+basis, tol=mp.mpf('1e-70'), maxcoeff=10**9, maxsteps=20000)
    if rel is None:
        print(tag, ": PSLQ FAILED")
        continue
    resid = mp.fsum(mp.mpf(c)*v for c,v in zip(rel, [I]+basis))
    print(f"\n{tag} = {mp.nstr(I,30)}")
    print(f"  PSLQ residual: {mp.nstr(abs(resid),3)}")
    for nm, c in zip(['I']+names, rel):
        if c != 0: print("   ", nm, ":", c)
