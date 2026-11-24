
# step_06d.py -- full M3(2pi/3) from R3 + Psi + Phi vs series value; W4 accounting at 1/sqrt2
import sys; sys.path.insert(0, '.')
import mpmath as mp, json, re
from hpl_num import hpl
mp.mp.dps = 120
def load(f):
    data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/" + f))
    out = []
    for s, c in data:
        m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
        if m:
            out.append((int(m.group(1)), tuple(int(t) for t in m.group(2).split(',')), mp.mpf(c)))
        else:
            out.append((None, None, mp.mpf(c)))   # constant term
    return out
R3 = load("hpl_fit_k3_(0,0,0).json")
PSI = load("hpl_fit_k3_(1,0,0).json")
def ev_terms(terms, c0):
    tot = mp.mpf(0)
    for e, w, cf in terms:
        if w is None:
            tot += cf
        else:
            tot += cf * (c0 ** e) * hpl(w, c0)
    return tot
def Phi(c):
    return mp.mpf(1)/2 - (1 - c**2)/(2*c) * mp.atanh(c)
def M3(c):
    return -mp.zeta(3)*Phi(c) + mp.pi**2/6 * ev_terms(PSI, c) + ev_terms(R3, c)
c = mp.mpf(1)/2
val = M3(c)
print("M3(2pi/3) via R3+Psi+Phi =", mp.nstr(val, 70))
print("series value               = 0.01071736293228946541840957629122489585571897573846433576015561...")
print("diff =", mp.nstr(val - mp.mpf('0.01071736293228946541840957629122489585571897573846433576015561'), 5))

# weight-4 part of M3 at c: Omega - zeta3*Phi_wt1 + zeta2*Psi_wt2
def Omega(c0):
    return mp.fsum(cf * (c0 ** e) * hpl(w, c0) for e, w, cf in R3 if w is not None and len(w) == 4)
def Psi_wt2(c0):
    return mp.fsum(cf * (c0 ** e) * hpl(w, c0) for e, w, cf in PSI if w is not None and len(w) == 2)
def Phi_wt1(c):
    return -(1 - c**2)/(2*c) * mp.atanh(c)
for cc, name in ((mp.mpf(1)/2, "2pi/3"), (1/mp.sqrt(2), "pi/2")):
    w4part = Omega(cc) - mp.zeta(3)*Phi_wt1(cc) + mp.pi**2/6*Psi_wt2(cc)
    print(f"weight-4 part of M3 at {name} =", mp.nstr(w4part, 50))
print("paper W4 = 0.79578505983224081132")
