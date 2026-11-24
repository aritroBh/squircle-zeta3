import json, math
from fractions import Fraction as F
M4E = '<project>/research/08_theory/M4e_w5_lchi'
open(M4E+'/progress.log','a').write(
 "step_28-31: setup; copied M4d scripts+jsons. L(chi8,5) to 1650 digits (mpmath Hurwitz), "
 "flint acb cross-check diff 4.2e-1656; 3-form check at 60d OK. Saved Lchi8_5.json\n")
sym = json.load(open(M4E+'/nice_symbol.json'))
kernel = [[F(c) for c in v] for v in sym['kernel']]
print("kernel vectors:", len(kernel))
def prim(vec):
    from math import gcd, lcm
    D = lcm(*[c.denominator for c in vec])
    ints = [int(c*D) for c in vec]
    g = 0
    for i in ints: g = gcd(g, abs(i))
    return [i//g for i in ints]
K = [prim(v) for v in kernel]
for k in K[:8]:
    print("max|entry|=", max(map(abs,k)))
import flint
t0=__import__('time').monotonic()
red = flint.fmpz_mat(K).lll(gram='exact')
Kr = [[int(red[i,j]) for j in range(red.ncols())] for i in range(red.nrows())]
print("LLL kernel reduce:", __import__('time').monotonic()-t0, "s")
lens = sorted(max(abs(i) for i in r) for r in Kr)
print("reduced kernel max-entries:", lens)
