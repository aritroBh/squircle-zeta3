
import sys, pickle, itertools, time
sys.path.insert(0,'.')
from fractions import Fraction as Fr
S = pickle.load(open('syms_omega_dt.pkl','rb'))
S_diff = S['S_diff']
B = pickle.load(open('solve81_data.pkl','rb')) if False else None

E = [(1,0,0),(0,1,0),(0,0,1)]
def expand(sym):
    out = {}
    for k,c0 in sym.items():
        opts = [[(i,vi) for i,vi in enumerate(letter) if vi] for letter in k]
        for combo in itertools.product(*opts):
            mult = Fr(c0)
            for i,vi in combo: mult *= vi
            nk = tuple(i for i,vi in combo)
            out[nk] = out.get(nk, Fr(0)) + mult
        out = {k:v for k,v in out.items() if v}
    return out
def vec81(sym):
    e = expand(sym)
    v = [Fr(0)]*81
    for k,c0 in e.items():
        v[k[0]*27 + k[1]*9 + k[2]*3 + k[3]] = c0
    return tuple(v)
def target_vec(e):
    out = {}
    for k, poly in S_diff.items():
        if any(letter[0] != 0 for letter in k): continue
        c0 = poly.get(e, Fr(0))
        if c0:
            nk = tuple(letter[1:] for letter in k)   # project to (c,p,q)
            out[nk] = out.get(nk, Fr(0)) + c0
    return vec81(out)
T = {e: target_vec(e) for e in (-1,0,1)}
print("target component counts:", {e: sum(1 for x in v if x) for e,v in T.items()})
pickle.dump(T, open('targets81.pkl','wb'))
