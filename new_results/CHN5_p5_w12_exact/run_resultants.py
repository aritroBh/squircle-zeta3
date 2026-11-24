
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
