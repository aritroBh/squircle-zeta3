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
def eval_Earc_full(valmap, dps=60):
    mp.mp.dps=dps
    xv = -valmap[kap]
    z2v,z3v,z5v = mp.zeta(2),mp.zeta(3),mp.zeta(5)
    avals={nm:ae.atom_eval(nm,xv) for nm in atom_names} if xv!=0 else None
    total=mp.mpf(0)
    for key,co in E5.items():
        n,a,b,c=key
        covn = mp.mpf(sp.N(co.subs({s: mp.mpf(valmap[s]) for s in valmap}), dps))
        if covn==0: continue
        if c>=1:
            if xv==0: continue
            Jv = mp.mpf(sp.N(J5[key].subs({atoms[nm]:avals[nm] for nm in atom_names}|{x:xv,z2:z2v,z3:z3v,z5:z5v}), dps))
        else:
            Jv = mp.mpf(sp.N(JLOC_n[(n,a,b)].subs({z2:z2v,z3:z3v,z5:z5v}), dps))
        total += covn*Jv
    return 120*total
# synthetic edge: B=1/2, L2=2, q0=-1, V02=1, weights=1; A1=1-B, A0=1-B-B*kap, logs=log(term coeffs)
Bv=mp.mpf(1)/2
print("E_arc(kappa) as kappa->0 (synthetic edge, should converge, NOT blow up):")
for j in [1,2,3,4,5]:
    kv = mp.mpf(10)**(-j)
    A1v=1-Bv; A0v=1-Bv-Bv*kv
    vm={A0:A0v,A1:A1v,B:Bv,kap:kv,L2:mp.mpf(2),q0:mp.mpf(-1),V02:mp.mpf(1),
        lA0:mp.log(A0v),lA1:mp.log(A1v),lB:mp.log(Bv)}
    ev=eval_Earc_full(vm)
    print(f"  kap=1e-{j}: E_arc = {mp.nstr(ev,25)}")
# kappa=0 local value
A1v=1-Bv; A0v=1-Bv
vm0={A0:A0v,A1:A1v,B:Bv,kap:mp.mpf(0),L2:mp.mpf(2),q0:mp.mpf(-1),V02:mp.mpf(1),
     lA0:mp.log(A0v),lA1:mp.log(A1v),lB:mp.log(Bv)}
e0=eval_Earc_full(vm0)
print(f"  kap=0 (local): E_arc = {mp.nstr(e0,25)}")
