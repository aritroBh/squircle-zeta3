import os, sys, time, pickle
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
from fractions import Fraction as F
import sympy as sp
from linalg_tools import ideal_member
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
A=[a1,a2,a3,a4,a5]
def dict_to_sym(p):
    return sp.expand(sum(sp.Integer(c.numerator)/c.denominator*sp.prod(A[i]**e[i] for i in range(5)) for e,c in p.items()))
Q=pickle.load(open('Q8_quotients.pkl','rb'))
QA=dict_to_sym(Q['A']); QZ=dict_to_sym(Q['Z'])
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
expr=sp.expand(QZ+sp.Rational(25,3)*QA)
print("expr==0:", expr==0, "terms:", len(sp.Poly(expr,*A).terms()))
t0=time.time()
mem=ideal_member(expr, [(C4,4),(Q7,6)], A, 8)
print("membership:", mem is not None, f"{time.time()-t0:.0f}s")
if mem:
    for k in mem:
        print("mem",k,"type",type(mem[k]))
    check=sp.expand(expr-(mem[0]*C4+mem[1]*Q7))
    print("verify U*C4+V*Q7' == expr:", check==0)
    print("U terms:", len(sp.Poly(mem[0],*A).terms()), "V terms:", len(sp.Poly(mem[1],*A).terms()) if not sp.sympify(mem[1]).is_number else 0)
    pickle.dump({'U':str(mem[0]),'V':str(mem[1])}, open('ideal_Z+25_3A.pkl','wb'))
    with open('ideal.log','a') as lg:
        lg.write(f"Z+25/3*A: IN ideal (C4,Q7')! verify=True [{time.time()-t0:.0f}s]\nDONE\n")
