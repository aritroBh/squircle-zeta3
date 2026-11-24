import json, mpmath as mp
from fractions import Fraction as Fr
BASE = "<project>/research/08_theory"
GOLD = BASE + "/GOLD5_golden_corner_w5/"
G2 = BASE + "/G2_M4alpha/"
DPS = 1450
mp.mp.dps = DPS + 40
flat = json.load(open(G2 + 'flat_exact_N640_K4.json'))
KEY = {'(0, 0, 0)': '1', '(1, 0, 0)': 'z2', '(0, 1, 0)': 'z3', '(0, 0, 1)': 'z4', '(2, 0, 0)': 'z2^2'}
mu = []
for nv in range(1, 641):
    mu.append({KEY[kk]: Fr(v) for kk, v in flat['4'][str(nv)].items()})
def flat_val(c):
    z2 = mp.pi**2/6; z3 = mp.zeta(3); z4 = mp.pi**4/90
    vals = {'1': mp.mpf(1), 'z2': z2, 'z3': z3, 'z4': z4, 'z2^2': z2**2}
    c2 = c*c; cp = c2; tot = mp.mpf(0)
    for nv in range(640):
        for kk, v in mu[nv].items():
            tot += mp.mpf(v.numerator)/v.denominator*vals[kk]*cp
        cp *= c2
    return tot
out = json.load(open(GOLD + 'values_golden.json'))
r5 = mp.sqrt(5)
for tag, c in [('gol1', (r5-1)/2), ('gol3', r5-2)]:
    M4 = mp.mpf(out[tag]['M4'])
    fv = flat_val(c)
    d = abs(M4 - fv)
    tail = abs(c*c)**641 * 100
    print(f"{tag}: |B4 - flat_N640(block 4)| = {mp.nstr(d, 3)}  (tail bound {mp.nstr(tail, 3)})")
    out[tag]['flat_diff'] = mp.nstr(d, DPS)
    out[tag]['flat_tail'] = mp.nstr(tail, DPS)
    out[tag]['flat_val_N640'] = mp.nstr(fv, DPS)
json.dump(out, open(GOLD + 'values_golden.json', 'w'), indent=1)
print("values_golden.json updated")
with open(GOLD + 'progress.log', 'a') as f:
    f.write("step_02 flat check corrected (block '4' only): " +
            "; ".join(f"{t} |B4-flat|={mp.nstr(abs(mp.mpf(out[t]['M4'])-flat_val((r5-1)/2 if t=='gol1' else r5-2)),3)}" for t in ('gol1','gol3')) + "\n")
