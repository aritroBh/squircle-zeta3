import os, pickle, time
import sympy as sp, mpmath as mp
import importlib.util
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
spec = importlib.util.spec_from_file_location("atom_eval", os.path.join(WD,'atom_eval.py'))
ae = importlib.util.module_from_spec(spec); spec.loader.exec_module(ae)
with open(f'{WD}/build_blocks.pkl','rb') as f: bb = pickle.load(f)
J5=bb['J5']; JLOC_n=bb['JLOC_n']
k1f = '<project>/research/08_theory/K1f_hexagon_order5'
with open(os.path.join(k1f,'E5_polys.pkl'),'rb') as f: E5 = pickle.load(f)[5]
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
x = sp.Symbol('x')
atom_names = ['H1','Li2','Li3','H011','Li4','H0011','H0111','Li5','H00011','H00101','H00111','H01011','H01111']
atoms = {nm: sp.Symbol(nm) for nm in atom_names}
z2,z3,z5 = sp.symbols('z2 z3 z5')
all_syms = [A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB]

def eval_Earc_numeric(ed, dps=80):
    """ed: dict of mpf values A0,A1,B,kap,L2,q0,V02,w0,w1,wo (weights as given)."""
    mp.mp.dps=dps
    xv = -ed['kap']
    # log values (unweighted): lA0=log(A0/w0) etc
    lA0v = mp.log(ed['A0']/ed['w0']); lA1v=mp.log(ed['A1']/ed['w1']); lBv=mp.log(ed['B']/ed['wo'])
    valmap = {A0:ed['A0'],A1:ed['A1'],B:ed['B'],kap:ed['kap'],L2:ed['L2'],q0:ed['q0'],V02:ed['V02'],
              lA0:lA0v,lA1:lA1v,lB:lBv}
    z2v,z3v,z5v = mp.zeta(2),mp.zeta(3),mp.zeta(5)
    if xv!=0:
        avals={nm:ae.atom_eval(nm,xv) for nm in atom_names}
    total=mp.mpf(0)
    for key,co in E5.items():
        n,a,b,c=key
        # numeric coefficient via evalf with substitutions
        cov = co.subs({s: mp.mpf(valmap[s]) for s in all_syms})
        covn = mp.mpf(sp.N(cov, dps))
        if covn==0: continue
        if c>=1:
            if xv==0: continue
            Je = J5[key].subs({atoms[nm]:avals[nm] for nm in atom_names} | {x:xv,z2:z2v,z3:z3v,z5:z5v})
            Jv = mp.mpf(sp.N(Je, dps))
        else:
            Jv = mp.mpf(sp.N(JLOC_n[(n,a,b)].subs({z2:z2v,z3:z3v,z5:z5v}), dps))
        total += covn*Jv
    return 120*total

with open(f'{WD}/golden_edges.pkl','rb') as f: GD = pickle.load(f)
edges=GD['edges']
t0=time.time()
tot=mp.mpf(0)
for ed in edges:
    edn = {kk: (mp.mpf(sp.N(ed[kk],90)) if kk not in ('w0','w1','wo') else mp.mpf(ed[kk])) for kk in ['A0','A1','B','kap','L2','q0','V02','w0','w1','wo']}
    ev = eval_Earc_numeric(edn)
    print(f"edge {ed['e']} (r={sp.nsimplify(ed['r'])}): E_arc = {mp.nstr(ev,30)}")
    tot += ev
P5_formula = 2*tot
print("\nFormula P^(5)(1) = 2*sum =", mp.nstr(P5_formula, 30))
print("Direct  P^(5)(1)          = -78.3585998739646...")
print("diff =", mp.nstr(P5_formula - mp.mpf('-78.3585998739646'), 5))
print("time %.1f s"%(time.time()-t0))
