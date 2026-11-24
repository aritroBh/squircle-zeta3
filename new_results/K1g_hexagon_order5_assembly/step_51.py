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
def eval_Earc_logs(ed, lA0v,lA1v,lBv, dps=80):
    mp.mp.dps=dps
    xv = -ed['kap']
    valmap = {A0:ed['A0'],A1:ed['A1'],B:ed['B'],kap:ed['kap'],L2:ed['L2'],q0:ed['q0'],V02:ed['V02'],
              lA0:lA0v,lA1:lA1v,lB:lBv}
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

# edge 1 data
ed1 = {'A0':mp.mpf(1)/2,'A1':mp.mpf(3)/4,'B':mp.mpf(1)/4,'kap':mp.mpf(1),
       'L2':mp.mpf(5),'q0':mp.mpf(1),'V02':mp.mpf(1)}
# |m_0|=1/4, |m_1|=(3/4)sqrt2, |m_2|=1/4
m0=mp.mpf(1)/4; m1=mp.mpf(3)/4*mp.sqrt(2); m2=mp.mpf(1)/4
# polar logs
plA0=mp.log(ed1['A0']/1); plA1=mp.log(ed1['A1']/3); plB=mp.log(ed1['B']/1)
# conv-B logs = polar - log|m_j|
blA0=plA0-mp.log(m0); blA1=plA1-mp.log(m1); blB=plB-mp.log(m2)
print("polar logs:", mp.nstr(plA0,10), mp.nstr(plA1,10), mp.nstr(plB,10))
print("convB logs:", mp.nstr(blA0,10), mp.nstr(blA1,10), mp.nstr(blB,10))
t0=time.time()
vB = eval_Earc_logs(ed1, blA0,blA1,blB)
print("edge1 conv-B E_arc =", mp.nstr(vB,45))
print("E5e value           = -852.2552563361562115016653047398929502035283608983064058739")
print("diff =", mp.nstr(vB - mp.mpf('-852.2552563361562115016653047398929502035283608983064058739'),5))
print("time %.1f"%(time.time()-t0))
