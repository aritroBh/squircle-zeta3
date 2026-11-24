
import mpmath as mp
import pickle
mp.mp.dps = 110
N = pickle.load(open('/tmp/Nfinal.pkl','rb'))['N']
M = pickle.load(open('/tmp/Mvals.pkl','rb'))['M']
pi=mp.pi

# ratios and normalized growth for M_n
mp.mp.dps = 30
print("n : M_n                 |M_n|^(1/n)      M_{n+1}/M_n")
for n in range(10, 40):
    print(f"{n:2d}  {mp.nstr(M[n],12):>15}  {mp.nstr(abs(M[n])**(1/n),8)}   {mp.nstr(M[n+1]/M[n],10)}")

# optimal truncation check at p = 1+eta, eta=0.05 (dual: h = eta/(1+eta) = 1/21)
mp.mp.dps = 110
def pi_cosh_p(p):
    p=mp.mpf(p); h=1/p
    f = lambda u: mp.cosh(u/2)**(-1-h) * mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f,[0,mp.inf])
for eta in [mp.mpf('0.05'), mp.mpf('0.02'), mp.mpf('0.1')]:
    p = 1+eta
    true = pi_cosh_p(p)/4
    # partial sums of M series
    s = mp.mpf(0); best=None; be=mp.mpf(1)
    for n in range(0,42):
        s += M[n]*eta**n
        err = abs(s-true)
        if err<be: be=err; best=n
    print(f"eta={eta}: pi_p/4={mp.nstr(true,25)}  best trunc at n={best}, err={mp.nstr(be,3)}")
