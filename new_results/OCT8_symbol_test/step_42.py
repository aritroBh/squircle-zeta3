import mpmath as mp, pickle, re
from fractions import Fraction as Fr
mp.mp.dps = 450
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator
d = pickle.load(open("symbol_identities.pkl","rb"))   # convention-A extraction (S_Li21)
G_names, D_names, sol_out = d["G_names"], d["D_names"], d["sol_out"]
logv = {p: mp.log(p) for p in (2,3,5)}
z3, z2 = mp.zeta(3), mp.zeta(2)

def RHS(name):
    c = [Fr(x) for x in sol_out[name]]
    tot = mp.mpf(0)
    for i, gn in enumerate(G_names):
        if c[i]: tot += M(c[i])*mp.polylog(3, M(Fr(gn[4:-1])))
    for j, dn in enumerate(D_names):
        cf = c[len(G_names)+j]
        if not cf: continue
        m = re.match(r"-Li2\((.*)\)\*log(\d+)", dn)
        w, p = Fr(m.group(1)), int(m.group(2))
        tot += M(cf)*(-mp.polylog(2, M(w)))*logv[p]
    return tot

def Li21A(x, y, N=2000):
    x = M(x); y = M(y); s = mp.mpf(0); yn = y
    for n in range(2, N):
        yn *= y
        t = mp.mpf(0); xm = x
        for m in range(1, n):
            t += xm/m**2; xm *= x
        s += t*yn/n
        if abs(yn) < mp.mpf(10)**(-mp.mp.dps-10): break
    return s

targetsA = {"Lam1": (Fr(-1,3), Fr(1,3)), "Lam3": (Fr(1,4), Fr(2,5))}
for name,(x,y) in targetsA.items():
    lhs = Li21A(x,y)
    rem = lhs - RHS(name)
    basis = [z3, z2*logv[2], z2*logv[3], z2*logv[5]]
    rel = mp.pslq([rem]+basis, tol=mp.mpf('1e-400'), maxcoeff=10**60, maxsteps=1000000)
    print(name, "convA: lhs =", mp.nstr(lhs,12), " |rem| =", mp.nstr(abs(rem),3), " PSLQ:", rel)
