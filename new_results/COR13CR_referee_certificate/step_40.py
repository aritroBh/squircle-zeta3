import json
from fractions import Fraction as Fr
base = '<project>/research/08_theory'
cert4 = json.load(open(base + '/COR13C_third_corner_certificate/cert13_w4.json'))
y = [Fr(s) for s in cert4["y"]]
hn = max(len(str(abs(yy.numerator))) for yy in y)
hd = max(len(str(yy.denominator)) for yy in y)
hmax = max(hn, hd)
import statistics
digs = sorted(max(len(str(abs(yy.numerator))), len(str(yy.denominator))) for yy in y)
print("max numerator digits:", hn, "max denominator digits:", hd)
print("median of max-digits:", digs[len(digs)//2])
nz = [yy for yy in y if yy != 0]
print("nonzero multipliers:", len(nz), "of", len(y))
# same for w3
cert23 = json.load(open(base + '/COR13C_third_corner_certificate/cert13_w23.json'))
for k in ('2','3'):
    yy = [Fr(int(a), int(b)) for a, b in cert23[k]["y"]]
    print(f"w{k}: max height digits:", max(max(len(str(abs(v.numerator))), len(str(v.denominator))) for v in yy),
          "nonzero:", sum(1 for v in yy if v != 0), "of", len(yy))
