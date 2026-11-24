
import json, re
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
def load_fit(name):
    D = json.load(open(G1 + "/" + name)); out = []
    for key, val in D:
        m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
        if key == "1": e, w = 0, ()
        elif key == "c": e, w = 1, ()
        else:
            e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
        out.append((e, w))
    return out
WORDS = sorted({w for name in ["hpl_fit_k3_(0,0,0).json", "hpl_fit_k3_(1,0,0).json"]
                for e, w in load_fit(name) if w}, key=lambda x: (len(x), x))
print(len(WORDS), "words")
print("ending with 0:", [w for w in WORDS if w[-1] == 0])
print("starting with 0:", len([w for w in WORDS if w[0] == 0]))
print(WORDS)
