import os, importlib.util
import mpmath as mp
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
spec = importlib.util.spec_from_file_location("atom_eval", os.path.join(WD,'atom_eval.py'))
ae = importlib.util.module_from_spec(spec); spec.loader.exec_module(ae)
base = '<project>/research/08_theory'
spec2 = importlib.util.spec_from_file_location("mzv", os.path.join(base,'E5c_nonlocal_z5','mzv.py'))
mzv = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(mzv)
mp.mp.dps=120
B = mzv.basis(dps=120)
labels=[b[0] for b in B]
atom_names = ['H1','Li2','Li3','H011','Li4','H0011','H0111','Li5','H00011','H00101','H00111','H01011','H01111']
import sympy as sp
atom_mzv = {}
print("PSLQ of atoms at x=-1 into alternating-MZV basis:")
for nm in atom_names:
    v = ae.atom_eval(nm, -1)
    rel = mzv.pslq_relation(v, B, tol=mp.mpf(10)**(-95), maxcoeff=10**9)
    if rel is None:
        print(f"  {nm}: PSLQ FAILED"); continue
    # rel: [coeff_X, coeff_b0, ...]; X = -sum(coeff_bi/coeff_X * bi)
    cX = rel[0]
    terms=[]
    expr=sp.Integer(0)
    for i,ci in enumerate(rel[1:]):
        if ci!=0:
            cf = sp.Rational(-ci, cX)
            terms.append((labels[i], cf))
            expr += cf*sp.Symbol(labels[i].replace('(','_').replace(')','').replace('/','_').replace(' ','').replace('^','p'))
    atom_mzv[nm]=terms
    print(f"  {nm}(-1) = " + " + ".join(f"({cf})*{lab}" for lab,cf in terms))
import pickle
with open(f'{WD}/atom_mzv_m1.pkl','wb') as f: pickle.dump(atom_mzv,f)
