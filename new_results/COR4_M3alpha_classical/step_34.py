import sys, json, re, time
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
mp.mp.dps = 60
from hpl_series import H_series
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
R3 = json.load(open(G1 + "hpl_fit_k3_(0,0,0).json"))
def parse(name, coef):
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
    num, den = (coef.split("/") + ["1"])[:2]
    return e, w, mp.mpf(int(num))/int(den)
c0 = mp.sqrt(2)/2
F = sum(cf * c0**e * H_series(w, c0) for n, c in R3 for e, w, cf in [parse(n, c)] if len(w) == 4)
sq2 = mp.sqrt(2); b = sq2 - 1
L = mp.log(1+sq2); lam = mp.log(2)
li4 = lambda x: mp.polylog(4, x); li3 = lambda x: mp.polylog(3, x); li2 = lambda x: mp.polylog(2, x)
z3 = mp.zeta(3); pi2 = mp.pi**2; pi4 = mp.pi**4
Li4part = ((1-573*sq2)*li4(b) - (1+573*sq2)*li4(-b) + 56*sq2*li4(-b**2)
           + 17*sq2*li4(1/sq2) - (mp.mpf(37)/4 + 17*sq2/16)*li4(mp.mpf(1)/2)
           + (10+86*sq2)*li4(1-1/sq2) - (10-86*sq2)*li4(2-sq2))
prods = (L*((40+7*sq2)*li3(b) - (40-7*sq2)*li3(-b) - 72*li3(1/sq2) + (mp.mpf(63)/8-3*sq2/16)*z3)
         + lam*((6+13*sq2)*li3(b) + (6-13*sq2)*li3(-b) - 47*sq2/2*li3(1/sq2) + (329*sq2/128-mp.mpf(3)/8)*z3)
         + 30*li2(b)**2 + li2(b)*((30+7*sq2/2)*L**2 - 24*lam*L - 19*sq2/8*lam**2 + (sq2/24-mp.mpf(19)/4)*pi2)
         + (mp.mpf(103)/384 - 4963*sq2/5760)*pi4
         + pi2*((mp.mpf(67)/24+185*sq2/96)*L**2 + (41*sq2/12-mp.mpf(35)/12)*lam*L + (mp.mpf(5)/12-1105*sq2/384)*lam**2)
         + (mp.mpf(13)/24-25*sq2/24)*L**4 - (mp.mpf(73)/6+29*sq2/3)*L**3*lam + (mp.mpf(63)/8+67*sq2/16)*L**2*lam**2
         + (mp.mpf(41)/24+17*sq2/6)*L*lam**3 + (137*sq2/192-mp.mpf(61)/192)*lam**4)
W4 = Li4part + prods
print("W4 =", mp.nstr(W4, 45))
print("exp= 0.795785059832240811320359819629726762546527413871")
D = F - W4
print("F  =", mp.nstr(F, 30))
print("D = F - W4 =", mp.nstr(D, 30))
