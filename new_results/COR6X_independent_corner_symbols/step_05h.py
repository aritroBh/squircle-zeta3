
import sys; sys.path.insert(0, '.')
import mpmath as mp, json, re
from hpl_num import hpl
mp.mp.dps = 120
data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
def parse_term(s):
    m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
    return int(m.group(1)), tuple(int(t) for t in m.group(2).split(','))
TERMS = [(parse_term(s)[0], parse_term(s)[1], mp.mpf(c)) for s, c in data if len(parse_term(s)[1]) == 4]
def Omega(c0):
    c0 = mp.mpf(c0)
    return mp.fsum(c * (c0 ** e) * hpl(w, c0) for e, w, c in TERMS)
O = Omega(1/mp.sqrt(2))
print("Omega(1/sqrt2) =", mp.nstr(O, 60))
print("paper W4       = 0.79578505983224081132...")
print("difference     =", mp.nstr(O - mp.mpf('0.79578505983224081132'), 5))
