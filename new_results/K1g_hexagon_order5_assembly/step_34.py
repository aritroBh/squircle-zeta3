import os, pickle, time
import sympy as sp
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
with open(f'{WD}/build_blocks.pkl','rb') as f:
    bb = pickle.load(f)
J5=bb['J5']; JLOC_n=bb['JLOC_n']; atoms=bb['atoms']; z2=bb['z2']; z3=bb['z3']; z5=bb['z5']; x=bb['x']
k1f = '<project>/research/08_theory/K1f_hexagon_order5'
with open(os.path.join(k1f,'E5_polys.pkl'),'rb') as f: E5 = pickle.load(f)[5]

basis_syms = set(atoms.values()) | {z2,z3,z5}
def decompose(J):
    """J (sympy in atoms,x,z's) -> list of (monomial_key, rational_coeff).
    monomial_key = (tuple sorted (atomname,exp)), x_exp, (z2e,z3e,z5e))."""
    out=[]
    for term in sp.Add.make_args(sp.expand(J)):
        pd = term.as_powers_dict()
        xe = pd.get(x,0)
        at=[]; ze=[0,0,0]; rat=sp.Integer(1)
        for s,e in pd.items():
            if s==x: continue
            elif s==z2: ze[0]=e
            elif s==z3: ze[1]=e
            elif s==z5: ze[2]=e
            elif s in basis_syms: at.append((str(s),e))
            else: rat *= s**e
        at=tuple(sorted(at))
        out.append(((at,xe,tuple(ze)), rat))
    return out

# Test decompose on a J5 entry
dc = decompose(J5[(0,0,1,1)])
print("decompose J5[(0,0,1,1)]:")
for mk,rc in dc: print("   ", mk, "*", rc)
