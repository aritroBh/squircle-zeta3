import os, pickle, time
import sympy as sp, mpmath as mp
import importlib.util
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
spec = importlib.util.spec_from_file_location("atom_eval", os.path.join(WD,'atom_eval.py'))
ae = importlib.util.module_from_spec(spec); spec.loader.exec_module(ae)
mp.mp.dps=80
with open(f'{WD}/coef_acc_raw.pkl','rb') as f: coef_acc = pickle.load(f)

A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
# Build a monomial->value function
atom_names = ['H1','Li2','Li3','H011','Li4','H0011','H0111','Li5',
              'H00011','H00101','H00111','H01011','H01111']

def eval_Earc(A0v,A1v,Bv,kapv,L2v,q0v,V02v,lA0v,lA1v,lBv, dps=80):
    mp.mp.dps = dps
    xv = -mp.mpf(kapv)   # x = -kap = 1-r
    # atom values
    avals = {nm: ae.atom_eval(nm, xv) for nm in atom_names}
    z2v=mp.zeta(2); z3v=mp.zeta(3); z5v=mp.zeta(5)
    subd = {A0:sp.nsimplify(A0v),A1:sp.nsimplify(A1v),B:sp.nsimplify(Bv),kap:sp.nsimplify(kapv),
            L2:sp.nsimplify(L2v),q0:sp.nsimplify(q0v),V02:sp.nsimplify(V02v)}
    # logs: keep as mp numbers, substitute into coefficient via evalf
    logd = {lA0: lA0v, lA1: lA1v, lB: lBv}
    total = mp.mpf(0)
    for (at,xe,ze),co in coef_acc.items():
        # numeric coefficient
        cov = co.subs(subd)
        covn = mp.mpf(sp.N(cov.subs({lA0:sp.Float(lA0v),lA1:sp.Float(lA1v),lB:sp.Float(lBv)}, 60)))
        if covn==0: continue
        # monomial value
        mv = mp.mpf(1)
        for nm,e in at: mv *= avals[nm]**e
        mv *= xv**xe
        mv *= z2v**ze[0]*z3v**ze[1]*z5v**ze[2]
        total += covn*mv
    return 120*total

# Test: b=0 edge 1, r=2
t0=time.time()
val = eval_Earc(A0v='1/2',A1v='3/4',Bv='1/4',kapv=1,L2v=5,q0v=1,V02v=1,
                lA0v=mp.log(mp.mpf(1)/2), lA1v=mp.log(mp.mpf(1)/4), lBv=mp.log(mp.mpf(1)/4))
print("b=0 edge 1 (r=2), K1f convention logs:")
print("E_arc =", mp.nstr(val, 40))
print("expected (E5e conv B) = -852.2552563361562115016653047398929502035")
print("diff =", mp.nstr(val - mp.mpf('-852.2552563361562115016653047398929502035'),5))
print("time %.1f s"%(time.time()-t0))
