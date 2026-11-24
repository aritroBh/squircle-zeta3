# Copied from research/08_theory/K1g_hexagon_order5_assembly/step_35.py (our notes) with the working directory taken from the environment variable G9WORK.
# Only paths and the module name of the e-sum engine were changed.
import os
import os, pickle, time
import sympy as sp
WD = os.environ['G9WORK']
with open(f'{WD}/build_blocks.pkl','rb') as f: bb = pickle.load(f)
J5=bb['J5']; JLOC_n=bb['JLOC_n']; atoms=bb['atoms']; z2=bb['z2']; z3=bb['z3']; z5=bb['z5']; x=bb['x']
k1f = os.environ['G9WORK']
with open(os.path.join(k1f,'E5_polys.pkl'),'rb') as f: E5 = pickle.load(f)[5]
basis_syms = set(atoms.values()) | {z2,z3,z5}
def decompose(J):
    out=[]
    for term in sp.Add.make_args(sp.expand(J)):
        pd = term.as_powers_dict()
        xe = pd.get(x,0); at=[]; ze=[0,0,0]; rat=sp.Integer(1)
        for s,e in pd.items():
            if s==x: continue
            elif s==z2: ze[0]=e
            elif s==z3: ze[1]=e
            elif s==z5: ze[2]=e
            elif s in basis_syms: at.append((str(s),e))
            else: rat *= s**e
        out.append(((tuple(sorted(at)),xe,tuple(ze)), rat))
    return out

t0=time.time()
# coefficient accumulator: monomial -> edge-data coefficient (sympy), for the SUM (before 120x)
coef_acc = {}
# Non-local
for key,co in E5.items():
    n,a,b,c = key
    if c>=1:
        for mk,rat in decompose(J5[key]):
            if mk in coef_acc: coef_acc[mk] = coef_acc[mk] + co*rat
            else: coef_acc[mk] = co*rat
    else:
        for mk,rat in decompose(JLOC_n[(n,a,b)]):
            if mk in coef_acc: coef_acc[mk] = coef_acc[mk] + co*rat
            else: coef_acc[mk] = co*rat
print("assembly accumulated, %.1f s, %d distinct monomials"%(time.time()-t0, len(coef_acc)))
import pickle as pk
with open(f'{WD}/coef_acc_raw.pkl','wb') as f: pk.dump(coef_acc,f)
print("saved coef_acc_raw.pkl")
# summarize monomials
npoles = {}
for mk in coef_acc:
    xe = mk[1]
    npoles[xe] = npoles.get(xe,0)+1
print("x-power distribution:", dict(sorted(npoles.items())))
