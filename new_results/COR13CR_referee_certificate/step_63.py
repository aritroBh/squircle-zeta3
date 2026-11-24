import json
from fractions import Fraction as Fr
base = '<project>/research/08_theory'
D = base + '/COR13C_third_corner_certificate'
# final consistency digest of numbers quoted in the report
c23 = json.load(open(D+'/cert13_w23.json')); c4 = json.load(open(D+'/cert13_w4.json')); ce = json.load(open(D+'/cert13_euler.json'))
print("cert sizes:", len(c23['2']['prow']), len(c23['3']['prow']), len(c4['prow']), "euler", len(ce['prow']))
y4 = [Fr(s) for s in c4['y']]
print("w4: nonzero", sum(1 for v in y4 if v), "max height digits",
      max(max(len(str(abs(v.numerator))), len(str(v.denominator))) for v in y4))
import statistics
hs = sorted(max(len(str(abs(v.numerator))), len(str(v.denominator))) for v in y4)
print("w4 median height digits:", hs[len(hs)//2])
for k in ('2','3'):
    y = [Fr(int(a),int(b)) for a,b in c23[k]['y']]
    print(f"w{k}: nonzero {sum(1 for v in y if v)}/{len(y)}, max height",
          max(max(len(str(abs(v.numerator))), len(str(v.denominator))) for v in y))
y = [Fr(int(a),int(b)) for a,b in ce['y']]
print("euler: nonzero", sum(1 for v in y if v), "max height", max(max(len(str(abs(v.numerator))), len(str(v.denominator))) for v in y))
print("cols:", len(c23['2']['cols']), len(c23['3']['cols']), len(c4['cols']))
