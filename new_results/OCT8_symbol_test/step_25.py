import mpmath as mp, pickle, re
from fractions import Fraction as Fr
mp.mp.dps = 420
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator

d = pickle.load(open("symbol_identities.pkl","rb"))
G_names, D_names, sol_out = d["G_names"], d["D_names"], d["sol_out"]
Lam_oct7 = pickle.load(open("../OCT7_generic_octagon/lam_lll_450.pkl","rb"))["Lam"]
Lam = {k: mp.mpf(v) for k, v in Lam_oct7.items()}

def Li21_series(x, y, N=1400):
    x = M(x); y = M(y); s = mp.mpf(0); yn = y
    for n in range(2, N):
        yn *= y
        t = mp.mpf(0); xm = x
        for m in range(1, n):
            t += xm/m**2; xm *= x
        s += t*yn/n
        if abs(yn) < mp.mpf(10)**(-mp.mp.dps-10): break
    return s
L1 = Li21_series(-Fr(1,3), Fr(1,3)); L3 = Li21_series(Fr(1,4), Fr(2,5))
print("cross-check Lam1:", mp.nstr(abs(L1-Lam["Lam1"]),3), " Lam3:", mp.nstr(abs(L3-Lam["Lam3"]),3))

def Li21_PV(x,y,eps=mp.mpf('1e-6')):
    x=M(x); y=M(y); t0=1/y
    res = -mp.polylog(2,x)
    f=lambda t: y*mp.polylog(2,x*y*t)/(1-y*t)
    def g(t):  # regularized
        return f(t) - res/(t-t0)
    a = mp.quad(f,[0,t0-eps]) + mp.quad(f,[t0+eps,1])
    b = mp.quad(g,[t0-eps,t0]) + mp.quad(g,[t0,t0+eps])
    return a + b + res*mp.log(eps/eps)
L2pv = Li21_PV(-Fr(1,9), 3); L4pv = Li21_PV(Fr(1,10), Fr(5,2))
print("PV Lam2 vs OCT7:", mp.nstr(abs(L2pv-Lam["Lam2"]),3))
print("PV Lam4 vs OCT7:", mp.nstr(abs(L4pv-Lam["Lam4"]),3))
pickle.dump({"Lam450":{k:mp.nstr(v,450) for k,v in Lam.items()}}, open("lam_vals.pkl","wb"))
