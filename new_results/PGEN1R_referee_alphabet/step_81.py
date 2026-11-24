from fractions import Fraction as Fr
import pickle
def wedge(uu, vv):
    out = {}
    for a, ca in uu.items():
        for bb, cb in vv.items():
            if a == bb: continue
            key, sg = ((a,bb),1) if a < bb else ((bb,a),-1)
            out[key] = out.get(key, Fr(0)) + sg*ca*cb
    return {k:c for k,c in out.items() if c}
def delta_li3(vz, vomz):
    w = wedge(vz, vomz)
    out = {}
    for ab, cab in w.items():
        for l3, c3 in vz.items():
            out[(ab, l3)] = out.get((ab, l3), Fr(0)) + cab*c3
    return {k:c for k,c in out.items() if c}
F = Fr
sf = pickle.load(open('../OCT9d_octagon_short/short_form.pkl','rb'))
G8 = [F(x) for x in sf['G8']]; B0 = [F(x) for x in sf['B0']]
sol0 = [F(x) for x in sf['sol0']]; sols123 = {e:[F(x) for x in v] for e,v in sf['sols123'].items()}
def vec_of(z):
    z = abs(F(z)); out = {}
    for p in (2,3,5,7):
        while z.numerator % p == 0: out[p] = out.get(p,0)+1; z /= p
        while z.denominator % p == 0: out[p] = out.get(p,0)-1; z *= p
    assert z == 1
    return {k:v for k,v in out.items() if v}
def delta_edge(args, coefs):
    acc = {}
    for z, c in zip(args, coefs):
        for k,v in delta_li3(vec_of(z), vec_of(1-z)).items():
            acc[k] = acc.get(k, Fr(0)) + c*v
    return {k:v for k,v in acc.items() if v}
mine = {1: delta_edge(G8, sols123[1][:8]), 2: delta_edge(G8, sols123[2][:8]),
        3: delta_edge(G8, sols123[3][:8]), 0: delta_edge(B0, sol0[:20])}
pg = pickle.load(open('../PGEN1_polygon_alphabet/delta_oct9d.pkl','rb'))
for e in [0,1,2,3]:
    print(f'edge {e}: mine == PGEN1:', mine[e]==pg[e], f' ({len(mine[e])} terms)')
    if mine[e]!=pg[e]:
        print('  mine:', mine[e])
# wedge^3 checks (edge 0 has 4 letters: real check)
def wedge3_image(acc):
    out = {}
    for ((a,b),l3), c in acc.items():
        if len({a,b,l3})<3: continue
        trip = (a,b,l3); srt = tuple(sorted(trip))
        inv = (trip[0]>trip[1]) + (trip[0]>trip[2]) + (trip[1]>trip[2])
        out[srt] = out.get(srt, Fr(0)) + c*(-1 if inv%2 else 1)
    return {k:v for k,v in out.items() if v}
for e in [0,1,2,3]:
    print(f'edge {e}: wedge^3(delta) =', wedge3_image(mine[e]))
open('progress.log','a').write('step: generic octagon delta per edge == PGEN1 exactly (0,1,2,3); wedge3=0\n')
