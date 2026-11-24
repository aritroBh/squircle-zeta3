import sys, json, pickle, re
sys.path.insert(0, ".")
from fractions import Fraction as Fr

base = "<project>/research"
R3 = json.load(open(base + "/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
results = pickle.load(open("word_classes.pkl", "rb"))

def parse(name, coef):
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", name)
    e = int(m.group(1))
    w = tuple(int(z) for z in m.group(2).split(","))
    num, den = (coef.split("/") + ["1"])[:2]
    return e, w, Fr(int(num), int(den))

tot_num = [Fr(0)]*3
tot_pipe = [Fr(0)]*3
nfail = 0
for name, coef in R3:
    e, w, cf = parse(name, coef)
    if len(w) != 4: continue
    res = results[w]
    cval = Fr(2) ** (-e)
    if res[0] == "FAIL":
        nfail += 1
        print("FAIL word:", w)
        continue
    nc, pc, h = res
    for j in range(3):
        tot_num[j] += cf * cval * nc[j]
        tot_pipe[j] += cf * cval * pc[j]
print("words that failed PSLQ:", nfail)
print("vF from NUMERIC classes:", tuple(tot_num), "  (expect (55/2, 13/2, 0))")
print("vF from pipeline classes:", tuple(tot_pipe), "  (was (5/2, 1/2, 0))")
