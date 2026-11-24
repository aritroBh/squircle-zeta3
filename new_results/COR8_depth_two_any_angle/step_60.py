
import json, sys, re
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
import mpmath as mp
mp.mp.dps = 60
from hpl_series import H_series
from fractions import Fraction as Fr
G1 = '<project>/research/10_frontier/G1_zeta3_any_angle/'
R3 = json.load(open(G1+'hpl_fit_k3_(0,0,0).json'))
def combo(cc, want=None):
    tot = mp.mpf(0)
    for name, coef in R3:
        m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
        e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(","))
        if want and len(w) != want: continue
        tot += mp.mpf(Fr(coef).numerator)/Fr(coef).denominator * cc**e * H_series(w, cc)
    return tot
cc = mp.mpf('0.5')
phi = mp.mpf('0.5') - ((1-cc**2)/(2*cc))*mp.atanh(cc)
def H(w, x): return H_series(w, x)
psi_val = (mp.mpf('0.5') + (mp.mpf('0.75') - 1/(4*cc) - 1/(2*cc**2))*H((1,),cc)
           + (-mp.mpf('0.75') - 1/(4*cc) + 1/(2*cc**2))*H((-1,),cc)
           + ((1-cc**2)/(4*cc))*(H((0,1),cc)+H((0,-1),cc)+H((1,-1),cc)-H((-1,1),cc)))
m3 = mp.mpf('0.010717362932289465418409576291224895855718975738464335760155612417841103619782418171396503742199126038436380620836337865304662784078999406937965739269009749615547401879802207222417678105864671523471037376024873555429679637766699958870900647844995948532283743592952407677723163070778394289744442351066844742744206712464681')
R3_check = m3 + mp.zeta(3)*phi - mp.zeta(2)*psi_val
print("R3 from M3 decomposition:", mp.nstr(R3_check, 40))
print("R3 from 74 HPL terms:    ", mp.nstr(combo(cc), 40))
print("Omega (w4 part):         ", mp.nstr(combo(cc, 4), 40))
