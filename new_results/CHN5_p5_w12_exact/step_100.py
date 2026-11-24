import os, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
with open('progress.log','a') as f:
    f.write(time.strftime('%H:%M:%S ')+"STEP3c trace2 (mpf, a1 1e-7..1e6) DONE: 119 interior pts, ALL QA<0, QC>0 (4 stratum artifacts at a1=1 excluded: a2=a4,a5~1e-61 where QA=QC=0 on stratum {a5=0,a2=a4,a1=a3}). Small-a1 asymptotic EXACT-linear: nQA=QA/S^8 ~ -8.16201e-6*a1 (a1=1e-7..1e-11), nQC ~ +1.606e-4*a1 -> clean boundary approach, no crossing. Large-a1: Gamma exists to a1>=1e9, nQA plateaus ~-1.89e-6,-2.33e-5. min|nQA| interior ~1e-7..1e-6, never 0. => NUMERICAL: A<0, CLASS>0 on all of Gamma+.\n")
    f.write(time.strftime('%H:%M:%S ')+"Cross-checks: 11 independent dense solves (a1 1e-5..3e6) all QA<0,QC>0. Branch counts 1-5/a1 (stratum degeneracy at a1=1). Approaching resultant elimination attempt for exact certificate.\n")
print("logged")
res = r'''
import sys, time, pickle
sys.path.insert(0,'.')
from fractions import Fraction as F
import sympy as sp
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
A=[a1,a2,a3,a4,a5]
def dict_to_sym(p):
    return sp.expand(sum(sp.Integer(c.numerator)/c.denominator*sp.prod(A[i]**e[i] for i in range(5)) for e,c in p.items()))
Q=pickle.load(open('Q8_quotients.pkl','rb'))
QA=dict_to_sym(Q['A']).subs(a3,1)
C4=sp.sympify(open('C4.txt').read()).subs(a3,1)
P7=sp.sympify(open('P7.txt').read()).subs(a3,1)
with open('resultants.log','w') as lg:
    for name,P in [('Res_a5(C4,P7)',P7), ('Res_a5(C4,QA)',QA)]:
        t0=time.time()
        R=sp.resultant(C4, P, a5)
        lg.write(f"{name}: computed [{time.time()-t0:.0f}s], terms={len(sp.Poly(R).terms())}, total deg={sp.Poly(R).total_degree()}\n"); lg.flush()
        pickle.dump(str(R), open(f'resultant_{"P7" if "P7" in name else "QA"}.pkl','wb'))
        t0=time.time()
        fl=sp.factor_list(R)
        lg.write(f"{name}: factor_list [{time.time()-t0:.0f}s] -> {len(fl[1])} factors:\n")
        for fac,mult in fl[1]:
            lg.write(f"   mult {mult}, deg {sp.Poly(fac).total_degree()}, terms {len(sp.Poly(fac).terms())}\n")
        lg.flush()
    lg.write("DONE\n")
print("resultants done")
'''
open(W+'/run_resultants.py','w').write(res)
import subprocess
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_resultants.py'],
                      cwd=W, stdout=open(W+'/run_resultants.out','w'), stderr=subprocess.STDOUT)
open(W+'/run_resultants.pid','w').write(str(proc.pid))
print("launched resultants pid", proc.pid, "at", time.strftime('%H:%M:%S'))
