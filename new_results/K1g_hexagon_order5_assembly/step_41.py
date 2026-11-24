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
atom_names = ['H1','Li2','Li3','H011','Li4','H0011','H0111','Li5',
              'H00011','H00101','H00111','H01011','H01111']
atoms = {nm: sp.Symbol(nm) for nm in atom_names}
z2,z3,z5 = sp.symbols('z2 z3 z5')

def eval_Earc_perkey(edge, dps=80):
    mp.mp.dps=dps
    kapv=edge['kap']; xv = -mp.mpf(sp.N(kapv,dps+10))
    rd = {A0:edge['A0'],A1:edge['A1'],B:edge['B'],kap:edge['kap'],
          L2:edge['L2'],q0:edge['q0'],V02:edge['V02'],
          lA0:sp.log(edge['A0']/edge['w0']), lA1:sp.log(edge['A1']/edge['w1']),
          lB:sp.log(edge['B']/edge['wo'])}
    z2v,z3v,z5v = mp.zeta(2),mp.zeta(3),mp.zeta(5)
    if xv != 0:
        avals = {nm: ae.atom_eval(nm,xv) for nm in atom_names}
        atomsub = dict(avals); atomsub[x]=xv; atomsub[z2]=z2v; atomsub[z3]=z3v; atomsub[z5]=z5v
    total = mp.mpf(0)
    for key,co in E5.items():
        n,a,b,c = key
        covn = mp.mpf(sp.N(co.subs(rd), dps))
        if covn==0: continue
        if c>=1:
            if xv==0: continue  # non-local vanishes at x=0
            Jexpr = J5[key]
            Jv = mp.mpf(sp.N(Jexpr.subs({atoms[nm]:avals[nm] for nm in atom_names} | {x:xv, z2:z2v, z3:z3v, z5:z5v}), dps))
        else:
            Jv = mp.mpf(sp.N(JLOC_n[(n,a,b)].subs({z2:z2v,z3:z3v,z5:z5v}), dps))
        total += covn*Jv
    return 120*total

# Test on edge 3 (kappa=0, local law)
b0e3 = {'A0':sp.Rational(1,4),'A1':sp.Rational(1,4),'B':sp.Rational(3,4),'kap':sp.Integer(0),
        'L2':sp.Integer(2),'q0':sp.Integer(-1),'V02':sp.Integer(1),'w0':sp.Integer(1),'w1':sp.Integer(1),'wo':sp.Integer(3)}
t0=time.time()
v3 = eval_Earc_perkey(b0e3)
print("b=0 edge3 (kappa=0): E_arc =", mp.nstr(v3,45))
print("expected (E5c conv B) 7.484865725290184402298838186623603542263891963")
print("diff:", mp.nstr(v3 - mp.mpf('7.484865725290184402298838186623603542263891963'),5))
print("time %.1f"%(time.time()-t0))
