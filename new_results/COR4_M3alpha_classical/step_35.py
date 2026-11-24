import sys, json, re
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
F = sum(cf * c0**e * H_series(w, c0) for n, c in R3 for e, w, cf in [parse(n, c)] if len(w) == 4)
D = F - W4

# explicit product/lower-weight expression for D:
Phi = mp.mpf(1)/2 - (1-c0**2)/(2*c0)*mp.atanh(c0)
Psi = (mp.mpf(1)/2 + (mp.mpf(3)/4 - 1/(4*c0) - 1/(2*c0**2))*H_series((1,),c0)
       + (-mp.mpf(3)/4 - 1/(4*c0) + 1/(2*c0**2))*H_series((-1,),c0)
       + ((1-c0**2)/(4*c0))*(H_series((0,1),c0)+H_series((0,-1),c0)+H_series((1,-1),c0)-H_series((-1,1),c0)))
R3low = sum(cf * c0**e * H_series(w, c0) for n, c in R3 for e, w, cf in [parse(n, c)] if len(w) < 4)
W3p = ((7*sq2/32-mp.mpf(3)/16)*z3 + (1+sq2)*li3(b) + (1-sq2)*li3(-b) - 2*sq2*li3(1/sq2)
       + (2*L-sq2*lam)*li2(b)
       + (sq2/12-mp.mpf(1)/8)*pi2*L - (1+sq2)/24*pi2*lam + (1+sq2)/24*lam**3 + sq2/4*lam**2*L
       - (1+sq2)/2*lam*L**2 + (5-sq2)/6*L**3)
T = (pi2 + 6*lam**2)/24
Dexp = z3*Phi - mp.zeta(2)*Psi - R3low + W3p + T
print("D          =", mp.nstr(D, 35))
print("D explicit =", mp.nstr(Dexp, 35))
print("|D - Dexp| =", mp.nstr(abs(D - Dexp), 3))
# now: D must be decomposable (weight<=3 + zeta(2)*weight2). Check that D is a combination of products:
# simplest: D - (z3*Phi - zeta2*Psi1 - R3low + W3p + T) should equal zeta(2)*Psi2(c0)
Psi2 = ((1-c0**2)/(4*c0))*(H_series((0,1),c0)+H_series((0,-1),c0)+H_series((1,-1),c0)-H_series((-1,1),c0))
print("|D - (W4products part)|: residual test below")
