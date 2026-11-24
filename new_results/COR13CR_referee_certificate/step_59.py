import sys, json
sys.dont_write_bytecode = True
base = '<project>/research/08_theory'
sys.path.insert(0, base + '/HEXC1_hexagon_certificate')
import cmzv6
from fractions import Fraction as Fr

def tag_families(n, lower):
    fam = {}
    def emit(row, f):
        r = cmzv6.primitive(row)
        if r and r not in fam: fam[r] = f
    if n == 1:
        for r in cmzv6.weight1_rows():
            rr = cmzv6.primitive(r)
            if rr: fam[rr] = 'weight1_branch'
    else:
        for a in range(1, n):
            b = n - a
            for u in cmzv6.words(a, False):
                u_adm = u[0] != 0
                if u_adm and a > b: continue
                for v in cmzv6.words(b):
                    if u_adm and a == b and u > v: continue
                    emit(cmzv6.ds(u, v), 'double_shuffle')
    for d in (2, 3, 6):
        for row in cmzv6.distribution(n, d): emit(row, f'distribution_d{d}')
    for k, basis in lower.items():
        for row in basis:
            for w in cmzv6.words(n - k):
                emit(cmzv6.mul(dict(row), {w: 1}), f'lift_w{k}')
    return fam

lower = {}
rows1, _ = cmzv6.make_relations(1, lower); lower[1] = rows1
rows2, _ = cmzv6.make_relations(2, lower); lower[2] = rows2
fam = tag_families(2, lower)
certE = json.load(open(base + '/COR13C_third_corner_certificate/cert13_euler.json'))
cert23 = json.load(open(base + '/COR13C_third_corner_certificate/cert13_w23.json'))
from collections import Counter
usageE = Counter(fam[rows2[i]] for i in certE["prow"])
usage2 = Counter(fam[rows2[i]] for i in cert23["2"]["prow"])
print("euler cert family usage:", dict(usageE))
print("w2 cert family usage:", dict(usage2))
yE = [Fr(int(a), int(b)) for a, b in certE["y"]]
print("euler max height digits:", max(max(len(str(abs(v.numerator))), len(str(v.denominator))) for v in yE),
      "nonzero:", sum(1 for v in yE if v != 0), "of", len(yE))
with open('progress.log','a') as f:
    f.write("referee step(verifier): verify_cor13cr.py fresh run 580s, ALL PASS (A exact + B numeric)\n")
print("done")
