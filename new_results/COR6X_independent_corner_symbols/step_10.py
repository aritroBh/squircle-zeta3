
# step_10.py -- task (1): confirm Omega(1/2) = five-Li4 + P at 200 digits (own code)
import sys; sys.path.insert(0, '<project>/research/08_theory/COR6X_independent_corner_symbols')
import mpmath as mp, json, re
from hpl_num import hpl
mp.mp.dps = 200
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[0], parse_term(s)[1], mp.mpf(c)) for s, c in data if len(parse_term(s)[1]) == 4]
def Omega(c0):
    c0 = mp.mpf(c0)
    return mp.fsum(c * (c0 ** e) * hpl(w, c0) for e, w, c in TERMS)
O = Omega(mp.mpf(1)/2)
L4 = lambda z: mp.polylog(4, z)
C = (-21*L4(-2) - mp.mpf(65)/8*L4(mp.mpf(1)/3) + mp.mpf(13)/2*L4(mp.mpf(2)/3)
     + mp.mpf(137)/6*L4(mp.mpf(1)/2) + mp.mpf(39)/64*L4(mp.mpf(1)/9))
l2, l3 = mp.log(2), mp.log(3)
z3 = mp.zeta(3)
L32 = mp.polylog(3, mp.mpf(1)/3); L3m = mp.polylog(3, -mp.mpf(1)/2)
L23 = mp.polylog(2, mp.mpf(1)/3)
p2 = mp.pi**2
P = (mp.mpf(15)/2*L32*l2 + mp.mpf(37)/2*L3m*l2 - mp.mpf(5)/4*L3m*l3
     + mp.mpf(119)/16*z3*l2 - mp.mpf(37)/32*z3*l3
     + L23**2 - mp.mpf(23)/4*L23*l2**2 + mp.mpf(7)/4*L23*l2*l3 + L23*l3**2
     - mp.mpf(1)/6*p2*L23
     - mp.mpf(455)/144*p2*l2**2 + mp.mpf(9)/8*p2*l2*l3 - mp.mpf(7)/32*p2*l3**2
     - mp.mpf(949)/2160*mp.pi**4
     - mp.mpf(211)/144*l2**4 + mp.mpf(101)/24*l2**3*l3 - mp.mpf(15)/8*l2**2*l3**2
     - mp.mpf(35)/24*l2*l3**3 + mp.mpf(113)/192*l3**4)
D = O - C
print("Omega(1/2)              =", mp.nstr(O, 50))
print("Omega - fiveLi4         =", mp.nstr(D, 50))
print("P (COR5 products)       =", mp.nstr(P, 50))
print("difference              =", mp.nstr(D - P, 6))
