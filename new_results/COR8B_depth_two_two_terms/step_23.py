
import json
d = json.load(open('../../10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
import re
def weight(term):
    m = re.search(r'H\[(.*)\]', term[0])
    return len(m.group(1).split(','))
for t in d:
    print(weight(t[0]), t[0], t[1])
