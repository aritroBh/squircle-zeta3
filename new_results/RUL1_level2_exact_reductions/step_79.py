import pathlib
code = r'''
import sys, time, pickle
sys.dont_write_bytecode = True
import sympy as sp
import mpmath as mp
mp.mp.dps = 80
K1G = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
K1F = '<project>/research/08_theory/K1f_hexagon_order5'
coef_acc = pickle.load(open(f'{K1G}/coef_acc_raw.pkl','rb'))
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
l2, l3 = sp.symbols('l2 l3')
# numeric atom values via Landen (independent)
import importlib.util, os
spec = importlib.util.spec_from_file_location("atom_eval", os.path.join(K1G,'atom_eval.py'))
ae = importlib.util.module_from_spec(spec); spec.loader.exec_module(ae)
atom_names = list(ae.ATOMWORDS.keys())
avals = {nm: ae.atom_eval(nm, mp.mpf(-1)) for nm in atom_names}
z2v,z3v,z5v = mp.zeta(2),mp.zeta(3),mp.zeta(5)
def foldnum(lA0v,lA1v,lBv):
    rd = {A0:mp.mpf(1)/2,A1:mp.mpf(3)/4,B:mp.mpf(1)/4,kap:mp.mpf(1),L2:mp.mpf(5),q0:mp.mpf(1),V02:mp.mpf(1),
          lA0:lA0v,lA1:lA1v,lB:lBv}
    tot = mp.mpf(0)
    for (at,xe,ze),co in coef_acc.items():
        cov = co.subs(rd)
        if cov == 0: continue
        mv = mp.mpf((-1)**xe)
        for nm,e in at: mv *= avals[nm]**e
        mv *= z2v**ze[0]*z3v**ze[1]*z5v**ze[2]
        tot += mp.mpf(sp.N(cov, 30))*mv
    return 120*tot
# conv-B logs
blA0 = mp.log(2); blA1 = -mp.log(2)/2 - mp.log(3); blB = mp.mpf(0)
vB = foldnum(blA0, blA1, blB)
print("numeric conv-B fold =", mp.nstr(vB, 30))
print("E5e value           = -852.2552563361562115016653047398929502035283608983064058739")
print("diff =", mp.nstr(vB - mp.mpf("-852.2552563361562115016653047398929502035283608983064058739"), 3))
# conv-B logs with log3 REPLACED by 0 to see the log3 content:
vB0 = foldnum(blA0, -mp.log(2)/2, blB)
print("log3->0 fold       =", mp.nstr(vB0, 30))
print("log3-content (vB-vB0) =", mp.nstr(vB - vB0, 6))
'''
pathlib.Path('rul1_step05d.py').write_text(code)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice',PY,'rul1_step05d.py'], capture_output=True, text=True, timeout=1200)
print(r.stdout[-1500:]); print(r.stderr[-800:])
