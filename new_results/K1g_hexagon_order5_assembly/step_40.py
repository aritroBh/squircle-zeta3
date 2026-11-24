import os, pickle, time
import sympy as sp, mpmath as mp
import importlib.util
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
spec = importlib.util.spec_from_file_location("atom_eval", os.path.join(WD,'atom_eval.py'))
ae = importlib.util.module_from_spec(spec); spec.loader.exec_module(ae)
with open(f'{WD}/coef_acc_raw.pkl','rb') as f: coef_acc = pickle.load(f)
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
atom_names = ['H1','Li2','Li3','H011','Li4','H0011','H0111','Li5',
              'H00011','H00101','H00111','H01011','H01111']

def eval_Earc(edge, dps=80):
    """edge: dict with sympy-exact A0,A1,B,kap,L2,q0,V02 (Rational) and r; weights w0,w1,wo."""
    mp.mp.dps = dps
    kapv = edge['kap']
    xv = -mp.mpf(sp.N(kapv, dps+10))
    avals = {nm: ae.atom_eval(nm, xv) for nm in atom_names}
    z2v=mp.zeta(2); z3v=mp.zeta(3); z5v=mp.zeta(5)
    # exact sympy substitution for rational data
    rd = {A0:edge['A0'],A1:edge['A1'],B:edge['B'],kap:edge['kap'],
          L2:edge['L2'],q0:edge['q0'],V02:edge['V02'],
          lA0:sp.log(edge['A0']/edge['w0']), lA1:sp.log(edge['A1']/edge['w1']),
          lB:sp.log(edge['B']/edge['wo'])}
    total = mp.mpf(0)
    for (at,xe,ze),co in coef_acc.items():
        cov = co.subs(rd)
        covn = mp.mpf(sp.N(cov, dps))
        if covn==0: continue
        mv = mp.mpf(1)
        for nm,e in at: mv *= avals[nm]**e
        mv *= xv**xe
        mv *= z2v**ze[0]*z3v**ze[1]*z5v**ze[2]
        total += covn*mv
    return 120*total

b0e1 = {'A0':sp.Rational(1,2),'A1':sp.Rational(3,4),'B':sp.Rational(1,4),'kap':sp.Integer(1),
        'L2':sp.Integer(5),'q0':sp.Integer(1),'V02':sp.Integer(1),'w0':sp.Integer(1),'w1':sp.Integer(3),'wo':sp.Integer(1)}
t0=time.time()
val = eval_Earc(b0e1)
print("b=0 edge1 (r=2), K1f conv logs, E_arc =")
print(mp.nstr(val, 45))
exp = mp.mpf('-852.2552563361562115016653047398929502035283608983064058739')
print("diff from E5e conv B:", mp.nstr(val-exp, 5))
print("time %.1f s"%(time.time()-t0))
