import pickle
K1d="research/08_theory/K1d_engine_fixed"
Jtab = pickle.load(open(K1d+'/Jtab_t.pkl','rb'))
ks = sorted(Jtab.keys())
print("all keys (n,a,b,c):", ks)
# verify a couple numerically
import mpmath as mp
mp.mp.dps=40
r2=mp.sqrt(2); k1=r2-1; x1v=-k1
import sympy as sp
L1,L2,L3,AA,z2,x = sp.symbols('L1 L2 L3 AA z2 x')
def atomvals(xv,N=3000):
    L1v=-mp.log(1-xv); L2v=mp.polylog(2,xv); L3v=mp.polylog(3,xv)
    s=mp.mpf(0); Hn=mp.mpf(0)
    for nn in range(1,N+1):
        Hn+=1/mp.mpf(nn); s+=xv**nn*Hn/nn**2
    return {L1:L1v,L2:L2v,L3:L3v,AA:s,z2:mp.pi**2/6,x:xv}
sub = atomvals(x1v)
ok=True
for key in ks:
    n,a,b,c = key
    f=lambda t: t**n*mp.log(t)**a*mp.log(1-t)**b*mp.log(1+k1*t)**c
    Jn = mp.quad(f,[0,1])
    Js = mp.mpf(sp.N(Jtab[key].subs(sub, simultaneous=True), 60))
    d = abs(Jn-Js)
    if d > mp.mpf('1e-30'): ok=False; print("MISMATCH", key, mp.nstr(Jn,20), mp.nstr(Js,20), mp.nstr(d,3))
print("Jtab check all keys at x1:", "OK" if ok else "FAIL")
