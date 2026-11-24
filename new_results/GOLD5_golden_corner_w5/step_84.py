import sys, itertools, json
GOLD = '<project>/research/08_theory/GOLD5_golden_corner_w5/'
sys.path.insert(0, GOLD)
from fractions import Fraction as Fr
import gold5_lib as G
Fraction = Fr
li5cl = {eval(k): tuple((Fr(q[0]), Fr(q[1])) for q in v)
         for k, v in json.load(open(GOLD + 'li5_classes_golden.json')).items()}
nz = [z for z, cl in li5cl.items() if any(any(q) for q in cl)]
adata = json.load(open(GOLD + 'atom_cobrackets_golden.json'))
def q12(v):
    out = []
    for q in v: out += [q[0], q[1]]
    return tuple(out)
def nfmul(a, b):
    return (a[0]*b[0] + 5*a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def nfadd(a, b):
    return (a[0]+b[0], a[1]+b[1])
def nfscale(q, a):
    return (q*a[0], q*a[1])
SQ5v = (Fr(0), Fr(1))
def cols_z(z):
    cl = li5cl[z]
    return [q12(cl), q12(tuple(nfmul(SQ5v, q) for q in cl))]

# class(Omega_5) gol1 (computed before); recompute to be safe
P = G.golden_point('gol1')
S5 = P.omega5()
cl_om = G.quotient5(G.ft_expand_sym_nf(S5, 2))
print("class(Omega5) gol1 =", [G.nf_str(q) for q in cl_om])

# atom symbols' classes: compute for all atoms used in the 710 solutions
keys = [k for k, v in adata.items() if any(Fr(c) for q in v['cvec'] for c in q)]
import pickle, os
if os.path.exists(GOLD + 'atom_classes_golden.json'):
    atomcl = json.load(open(GOLD + 'atom_classes_golden.json'))
else:
    atomcl = {}
    F = P.F
    for k in keys:
        x, y = eval(adata[k]['xy'][0]), eval(adata[k]['xy'][1])
        S_at = G.li32_sym_nf(P, x, y)
        cl = G.quotient5(G.ft_expand_sym_nf(S_at, 2))
        atomcl[k] = [[str(c) for c in q] for q in cl]
    json.dump(atomcl, open(GOLD + 'atom_classes_golden.json', 'w'), indent=1)
atomcl = {k: tuple((Fr(q[0]), Fr(q[1])) for q in v) for k, v in atomcl.items()}
print(f"{len(atomcl)} atom classes")
pickle_ok = True
