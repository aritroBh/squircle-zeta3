
import json, sys, re
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
import mpmath as mp
mp.mp.dps = 60
from hpl_series import H_series, Li22
from fractions import Fraction as Fr
G1 = '<project>/research/10_frontier/G1_zeta3_any_angle/'
R3 = json.load(open(G1+'hpl_fit_k3_(0,0,0).json'))
def combo(cc, upto=None):
    tot = mp.mpf(0)
    for name, coef in R3:
        m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
        e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
        if upto and len(w) != upto: continue
        tot += mp.mpf(Fr(coef).numerator)/Fr(coef).denominator * cc**e * H_series(w, cc)
    return tot
cc = mp.mpf('0.5')
# R3(0.5) via M3(2pi/3) = -zeta3 Phi + zeta2 Psi + R3
Psi = json.load(open(G1+'hpl_fit_k3_(1,0,0).json'))
sys.path.insert(0, G1)
import hpl_eval as he
phi = he.Phi(cc)
psi_val = mp.mpf(0)
for name, coef in Psi:
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
    # use classical dict evaluation via m3_hpl_check's H: fall back to H_series for len<=3 via classical_dict
    psi_val += mp.mpf(Fr(coef).numerator)/Fr(coef).denominator * cc**e * H_series(w, cc)
m3 = he.stored('2pi_3', 3)
R3_check = m3 + mp.zeta(3)*phi - mp.zeta(2)*psi_val
print("R3 from M3 decomposition:", mp.nstr(R3_check, 30))
print("R3 from 74 HPL terms:    ", mp.nstr(combo(cc), 30))
print("Omega (w4 part):         ", mp.nstr(combo(cc, 4), 30))
print("w2 part:", mp.nstr(combo(cc, 2), 30), " w3 part:", mp.nstr(combo(cc, 3), 30))
print("w2+w3+w4 - R3:", mp.nstr(combo(cc,2)+combo(cc,3)+combo(cc,4)-combo(cc), 5))
