import json, pickle
G2 = "<project>/research/08_theory/G2_M4alpha/"
import re
for fn in ["hpl_fit_k4_(0,0,0).json", "hpl_fit_k4_(1,0,0).json"]:
    d = json.load(open(G2+fn))
    exps = set(); wl = {}
    kinds = set()
    for nm, v in d:
        if nm in ('1','c'): kinds.add(nm); continue
        m = re.match(r'c\^(-?\d+) (He?)\[(.*)\]', nm)
        e = int(m.group(1)); w = tuple(int(t) for t in m.group(2).split(','))
        exps.add(e); wl[len(w)] = wl.get(len(w), 0) + 1; kinds.add(m.group(2))
    print(fn, "| exponents:", sorted(exps), "| words by weight:", wl, "| kinds:", kinds)
# mu4_zsum.pkl (the stored output of step_11)
z = pickle.load(open(G2 + "mu4_zsum.pkl", "rb"))
print("mu4_zsum.pkl parts/words:", {p: len(E) for p, E in z.items()})
# which Z-indices appear
idx = set()
for p, E in z.items():
    for w, c in E:
        for a in w: idx.add(a)
print("Z-indices used:", sorted(idx), "depths:", sorted(set(len(w) for p,E in z.items() for w,c in E)))
