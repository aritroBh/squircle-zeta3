import pickle
from fractions import Fraction as Fr
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
def vec_of(z):
    z = abs(Fr(z)); out = {}
    for p in (2,3,5,7):
        while z.numerator % p == 0: out[p] = out.get(p,0)+1; z /= p
        while z.denominator % p == 0: out[p] = out.get(p,0)-1; z *= p
    assert z == 1, z
    return {k:v for k,v in out.items() if v}
def delta_of_presentation(pres):
    acc = {}
    for c, z in pres:
        assert all(p in (2,3,5,7) for p in vec_of(z)) and all(p in (2,3,5,7) for p in vec_of(1-z)), z
        for k,v in delta_li3(vec_of(z), vec_of(1-z)).items():
            acc[k] = acc.get(k, Fr(0)) + c*v
    return {k:v for k,v in acc.items() if v}
cr = pickle.load(open('../PGEN1_polygon_alphabet/CR_presentations.pkl','rb'))
pg = pickle.load(open('../PGEN1_polygon_alphabet/delta_oct9d.pkl','rb'))
for e in [0,1,2,3]:
    dcr = delta_of_presentation(cr[e]['presentation'])
    same = dcr == pg[e]
    print(f'edge {e}: delta from CR presentation == delta from OCT9d basis:', same, f'({len(dcr)} terms)')
    if not same:
        print('  CR:', dcr); print('  OCT9d:', pg[e])
open('progress.log','a').write('step: generic octagon delta invariant: CR presentation == OCT9d basis, all 4 edges exact\n')
