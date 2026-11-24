base = "<project>/research/10_frontier/G1_zeta3_any_angle/"
import json
R3 = json.load(open(base+'hpl_fit_k3_(0,0,0).json'))
print(len(R3))
for t in R3[:10]: print(t)
# also check hpl_eval.py for the conventions
with open(base+"hpl_eval.py") as f:
    txt = f.read()
print(txt[:3000])
