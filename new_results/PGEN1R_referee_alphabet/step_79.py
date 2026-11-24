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
def F(s): return Fr(s)
sf = pickle.load(open('../OCT9d_octagon_short/short_form.pkl','rb'))
G8 = [F(x) for x in sf['G8']]
B0 = [F(x) for x in sf['B0']]
sol0 = [F(x) for x in sf['sol0']]
sols123 = {e: [F(x) for x in v] for e,v in sf['sols123'].items()}
print('G8:', G8)
print('B0:', B0)
print('sol0[:20]:', sol0[:20])
print('sols123[2][:8]:', sols123[2][:8])
print('sols123[3][:8]:', sols123[3][:8])
def vec_of(z):
    z = abs(F(z)); out = {}
    for p in (2,3,5,7):
        while z.numerator % p == 0: out[p] = out.get(p,0)+1; z /= p
        while z.denominator % p == 0: out[p] = out.get(p,0)-1; z *= p
    assert z == 1, z
    return {k:v for k,v in out.items() if v}
def delta_edge(args, coefs):
    acc = {}
    for z, c in zip(args, coefs):
        dd = delta_li3(vec_of(z), vec_of(1-z))
        for k,v in dd.items():
            acc[k] = acc.get(k, Fr(0)) + c*v
    return {k:v for k,v in acc.items() if v}
# edges 1-3: coefficients = first 8 of sols123[e]
d1 = delta_edge(G8, sols123[1][:8])
d2 = delta_edge(G8, sols123[2][:8])
d3 = delta_edge(G8, sols123[3][:8])
d0 = delta_edge(B0, sol0[:20])
print('edge 1 delta terms:', len(d1), ' != 0:', bool(d1))
print('edge 2 delta terms:', len(d2), ' != 0:', bool(d2))
print('edge 3 delta terms:', len(d3), ' != 0:', bool(d3))
print('edge 0 delta terms:', len(d0), ' != 0:', bool(d0))
pg = pickle.load(open('../PGEN1_polygon_alphabet/delta_oct9d.pkl','rb'))
print('PGEN1 delta_oct9d keys:', list(pg.keys()) if isinstance(pg, dict) else type(pg))
