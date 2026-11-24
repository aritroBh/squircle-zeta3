import mpmath as mp, pickle, re
from fractions import Fraction as Fr
mp.mp.dps = 420

# ---- load symbol identities ----
d = pickle.load(open("symbol_identities.pkl","rb"))
G_names, D_names, sol_out = d["G_names"], d["D_names"], d["sol_out"]

def parse_Li3(name):
    z = re.match(r"Li3\((.*)\)", name).group(1)
    return Fr(z)
def parse_prod(name):
    # '-Li2(w)*logp'  or  'loga*logb*logc'
    if name.startswith("-Li2("):
        m = re.match(r"-Li2\((.*)\)\*log(\d+)", name)
        return ('Li2log', Fr(m.group(1)), int(m.group(2)))
    ps = [int(x[3:]) for x in name.split("*")]
    return ('logs', ps)

# ---- numeric values ----
Lam_oct7 = pickle.load(open("../OCT7_generic_octagon/lam_lll_450.pkl","rb"))["Lam"]
Lam = {k: mp.mpf(v) for k, v in Lam_oct7.items()}

# independent cross-checks
def Li21_series(x, y, N=1200):
    x = mp.mpf(x); y = mp.mpf(y); s = mp.mpf(0); yn = y
    for n in range(2, N):
        yn *= y
        t = mp.mpf(0); xm = x
        for m in range(1, n):
            t += xm/m**2; xm *= x
        s += t*yn/n
        if abs(yn) < mp.mpf(10)**(-mp.mp.dps-20): break
    return s
L1 = Li21_series(Fr(-1,3), Fr(1,3))
L3 = Li21_series(Fr(1,4), Fr(2,5))
print("cross-check Lam1:", mp.nstr(abs(L1-Lam["Lam1"]),3), " Lam3:", mp.nstr(abs(L3-Lam["Lam3"]),3))
# PV check for Lam2, Lam4 (pole at t=1/y in (0,1))
def Li21_PV(x,y):
    x=mp.mpf(x); y=mp.mpf(y); t0=mp.mpf(1)/y
    f=lambda t: y*mp.polylog(2,x*y*t)/(1-y*t)
    a = mp.quad(f,[0,t0-mp.mpf('1e-8')]) + mp.quad(f,[t0+mp.mpf('1e-8'),1])
    b = (mp.quad(lambda t: f(t)-(-mp.polylog(2,x))/(t-t0),[t0-mp.mpf('1e-8'),t0]) +
         mp.quad(lambda t: f(t)-(-mp.polylog(2,x))/(t-t0),[t0,t0+mp.mpf('1e-8')]) )
    res = -mp.polylog(2,x)  # residue of f at t0
    # principal value: a + b + res*log((1-t0+... )) -> use symmetric log ratio
    return a + b + res*mp.log((mp.mpf('1e-8'))/(mp.mpf('1e-8')))
L2pv = Li21_PV(Fr(-1,9), Fr(3,1))
L4pv = Li21_PV(Fr(1,10), Fr(5,2))
print("PV Lam2 vs OCT7:", mp.nstr(abs(L2pv-Lam["Lam2"]),3), " Lam4:", mp.nstr(abs(L4pv-Lam["Lam4"]),3))
