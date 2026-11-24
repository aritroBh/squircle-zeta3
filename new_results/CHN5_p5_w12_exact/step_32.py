import os
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
ideal = r'''
import sys, time, pickle
sys.path.insert(0,'.')
from fractions import Fraction as F
import sympy as sp
from linalg_tools import ideal_member, monoms_homog, solve_linear
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
A=[a1,a2,a3,a4,a5]
def dict_to_sym(p):
    return sp.expand(sum(sp.Integer(c.numerator)/c.denominator*sp.prod(A[i]**e[i] for i in range(5)) for e,c in p.items()))
Q=pickle.load(open('Q8_quotients.pkl','rb'))
QA=dict_to_sym(Q['A']); QB=dict_to_sym(Q['B']); QC=dict_to_sym(Q['CLASS']); QZ=dict_to_sym(Q['Z'])
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
with open('ideal.log','w') as lg:
    for name,expr in [('B-75/14*A', QB-sp.Rational(75,14)*QA),
                      ('Z+25/3*A',  QZ+sp.Rational(25,3)*QA)]:
        expr=sp.expand(expr)
        if expr==0:
            lg.write(f"{name}: quotient difference is IDENTICALLY 0 (proportionality is a polynomial identity!)\n"); continue
        t0=time.time()
        mem=ideal_member(expr, [(C4,4),(Q7,6)], A, 8)
        if mem is None:
            lg.write(f"{name}: NOT in ideal (C4,Q7') [{time.time()-t0:.0f}s]\n")
        else:
            # verify
            check=sp.expand(expr-(mem[0]*C4+mem[1]*Q7))
            lg.write(f"{name}: IN ideal (C4,Q7')! U deg4 {len(sp.Poly(mem[0]).terms())} terms, V deg2 {len(sp.Poly(mem[1]).terms())} terms, verify={check==0} [{time.time()-t0:.0f}s]\n")
            pickle.dump({'U':str(mem[0]),'V':str(mem[1])}, open(f'ideal_{name.replace("/","_").replace("*","")}.pkl','wb'))
        lg.flush()
    lg.write("DONE\n")
print("ideal done")
'''
open(W+'/run_ideal.py','w').write(ideal)
print("run_ideal.py written")
import time
time.sleep(120)
print(open('eval12_0.log').read())
