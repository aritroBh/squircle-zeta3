import sys, json, re
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
mp.mp.dps = 80
from hpl_series import H_series
from fractions import Fraction as Fr

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
R3 = json.load(open(G1 + "hpl_fit_k3_(0,0,0).json"))
def parse(name, coef):
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
    num, den = (coef.split("/") + ["1"])[:2]
    return e, w, mp.mpf(int(num))/int(den)

c0 = mp.sqrt(2)/2
F = mp.mpf(0)
for n, c in R3:
    e, w, cf = parse(n, c)
    if len(w) == 4:
        F += cf * c0**e * H_series(w, c0)
print("F(1/sqrt2) =", mp.nstr(F, 40))

sq2 = mp.sqrt(2)
b = sq2 - 1
L = mp.log(2); lam = mp.log(1 + sq2)
li4 = lambda x: mp.polylog(4, x)
li3 = lambda x: mp.polylog(3, x)
li2 = lambda x: mp.polylog(2, x)
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
print("W4(paper formula) =", mp.nstr(W4, 40))
M3_paper = mp.mpf('0.022770521491786577456730509742327800250183691729217')
print("M3 - W3' - T vs W4 consistency: computing W3', T ...")
# W3' and T from thm:M3
W3p = ((7*sq2/32-mp.mpf(3)/16)*z3 + (1+sq2)*li3(b) + (1-sq2)*li3(-b)
       + (mp.mpf(0)) )
# need full W3' formula; instead check M3 formula from G1 at c0:
# M3 = -z3*Phi + z2*Psi + R3
Phi = mp.mpf(1)/2 - (1-c0**2)/(2*c0)*mp.atanh(c0)
def H(w, x): return H_series(w, x)
Psi = (mp.mpf(1)/2 + (mp.mpf(3)/4 - 1/(4*c0) - 1/(2*c0**2))*H((1,),c0)
       + (-mp.mpf(3)/4 - 1/(4*c0) + 1/(2*c0**2))*H((-1,),c0)
       + ((1-c0**2)/(4*c0))*(H((0,1),c0)+H((0,-1),c0)+H((1,-1),c0)-H((-1,1),c0)))
R3val = mp.mpf(0)
for n, c in R3:
    e, w, cf = parse(n, c)
    R3val += cf * c0**e * (mp.mpf(1) if len(w)==0 else H(w, c0))
M3_G1 = -z3*Phi + mp.zeta(2)*Psi + R3val
print("M3 (G1 formula) =", mp.nstr(M3_G1, 40))
print("|M3_G1 - M3_paper| =", mp.nstr(abs(M3_G1 - M3_paper), 3))
D = F - Li4part
print("D = F - Li4part =", mp.nstr(D, 30))
# W4 should equal M3 - W3' - T; and F - Li4part should be a combo of products
