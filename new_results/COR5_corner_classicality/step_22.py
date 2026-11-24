import os, json
base = "<project>/research"
g1 = os.path.join(base, "10_frontier/G1_zeta3_any_angle")
R3 = json.load(open(os.path.join(g1, "hpl_fit_k3_(0,0,0).json")))
Psi = json.load(open(os.path.join(g1, "hpl_fit_k3_(1,0,0).json")))
print("R3 terms:", len(R3), "Psi terms:", len(Psi))
print("R3[:8]:", json.dumps(R3[:8]))
print()
print("Psi:", json.dumps(Psi))
# weights of R3 terms
import re
def parse(name):
    m = re.match(r'c\^(-?\d+) H\[(.*)\]', name)
    if name == '1': return 0, ()
    if name == 'c': return 1, ()
    return int(m.group(1)), tuple(int(z) for z in m.group(2).split(','))
from collections import Counter
cnt = Counter(len(parse(n)[1]) for n, c in R3)
print("R3 weight counts:", cnt)
cnt2 = Counter(parse(n)[0] for n, c in R3)
print("R3 c-exponent counts:", cnt2)
