import sys, os, pickle, json
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
from flint import acb, ctx
import dilog_exact as DE

mp.mp.dps = 1250
ctx.prec = 4200

tails3 = pickle.load(open("m3_tails_exact.pkl", "rb"))
tails4 = pickle.load(open("m4_tails_exact.pkl", "rb"))

base_fl = {"lam": acb(2).log(), "L": (acb(2).sqrt()+1).log(), "pi2": acb.pi()**2, "z3": acb(3).zeta()}
def pval_fl(P):
    tot = acb(0)
    for k, c in P.d.items():
        x = acb(c.p.numerator)/c.p.denominator + acb(2).sqrt()*acb(c.q.numerator)/c.q.denominator
        for a in k:
            if isinstance(a, str):
                x *= base_fl[a]
            else:
                z = DE.float_exact(a[1])
                x *= acb(mp.nstr(z, 1200)).polylog(acb(2 if a[0]=="Li2" else 3))
        tot += x
    return tot

W3_k3 = pval_fl(tails3["W3p"])
W3_k4 = pval_fl(tails4["W3p"])

# rem:k4 claim: W3^(4) = W3'(k3) - z3/8 - pi2*lam/24 - lam^3/12
claim = W3_k3 - base_fl["z3"]/8 - base_fl["pi2"]*base_fl["lam"]/24 - base_fl["lam"]**3/12
print("W3^(4) (pickle, flint)        =", W3_k4.real.mid().str(80, radius=False, more=True))
print("claim (W3'k3 - z3/8 - ...)    =", claim.real.mid().str(80, radius=False, more=True))
d1 = W3_k4 - claim
print("pickle - claim =", d1.real.mid().str(8, radius=False, more=True), " rad", d1.real.rad().str(5))

# quadrature W3^(4) high precision
r2 = mp.sqrt(2); aN = 1+r2; bN = r2-1
omega = json.load(open("omega4.json"))
terms = {tuple(map(int, k.split(","))): [(z, acb(2).sqrt()*c[1] + acb(c[0])) for z, c in v] for k, v in omega.items()}
pole = {"0": acb(0), "a": acb(2).sqrt()+1, "-a": -acb(2).sqrt()-1, "b": acb(2).sqrt()-1, "-b": 1-acb(2).sqrt()}
def piece_acb(i, j):
    ts = terms[(i, j)]
    def f(v, analytic=False):
        l1 = ((v-bN)*(v+aN)/(4*v)).log(analytic=True)
        l2 = ((aN-v)*(v+bN)/(4*v)).log(analytic=True)
        s = acb(0)
        for z, c in ts:
            s += c/(v - pole[z])
        return l1**i * l2**j * s
    aNacb = acb(2).sqrt()+1; bNacb = acb(2).sqrt()-1
    return acb.integral(f, bNacb, 1, aNacb)
q = piece_acb(1,1) + piece_acb(2,0)
print("W3^(4) (acb quadrature)       =", q.real.mid().str(80, radius=False, more=True))
d2 = q - W3_k4
print("acb-quad - pickle =", d2.real.mid().str(8, radius=False, more=True), " rad", d2.real.rad().str(5))
