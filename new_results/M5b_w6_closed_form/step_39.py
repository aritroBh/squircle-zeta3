import sys, os, json
sys.path.insert(0,'<project>/research/08_theory/M5b_w6_closed_form')
os.chdir('<project>/research/08_theory/M5b_w6_closed_form')
import mpmath as mp
mp.mp.dps = 200
r2 = mp.sqrt(2); aN, bN = 1+r2, r2-1
om = json.load(open('omega5.json'))
pole = {"0": mp.mpf(0), "a": aN, "-a": -aN, "b": bN, "-b": -bN}
def piecequad(i, j, lo, hi):
    ts = om[f"{i},{j}"]
    def f(v):
        l1 = mp.log((v-bN)*(v+aN)/(4*v)); l2 = mp.log((aN-v)*(v+bN)/(4*v))
        return l1**i*l2**j*mp.fsum((mp.mpf(c[0])+r2*mp.mpf(c[1]))/(v-pole[z]) for z,c in ts)
    return mp.quad(f, [lo, (lo+hi)/2, hi])
parts = {}
for name, lo, hi in [('left', bN, 1), ('mid', 1, aN-1), ('right', aN-1, aN)]:
    s = mp.mpf(0)
    for ij in [(3,2),(4,1),(5,0)]:
        s += piecequad(*ij, lo, hi)
    parts[name] = s
    print(name, mp.nstr(s, 40))
print('total', mp.nstr(sum(parts.values()), 40))
# now the acb pieces from checkpoint
d = json.load(open('w6_600_partial.json'))
for k in ['left','mid','right']:
    print('acb', k, d[k][0][:50])
