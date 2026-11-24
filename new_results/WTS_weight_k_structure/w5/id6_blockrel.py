# Q-linear relations among the 10 block integrals and sqrt2 * blocks (20 numbers), LLL exact Gram, D digits.
import sys, json
sys.dont_write_bytecode = True
from id_lib import *
D = int(sys.argv[1])
setprec(1450)
print("flint ctx.prec =", ctx.prec, "bits; D =", D)
B = json.load(open("results/blocks_D1500_m10.json"))["blocks"]
sq2 = arb(2).sqrt()
names = list(B.keys())
vals = [arb(B[n]["mid"]) for n in names]
allv = vals + [v * sq2 for v in vals]
alln = names + ["sq2*" + n for n in names]
# also include 1-dim sanity: sizes
best, top = lll_relation(allv, D)
print("best row maxcoef", best[2], "last", best[1])
for row in top:
    print({n: c for n, c in zip(alln, row[:-1]) if c != 0}, "last", row[-1])
