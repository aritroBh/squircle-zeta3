
import sys, time, pickle
sys.path.insert(0,'.')
from fractions import Fraction as F
import sympy as sp
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
A=[a1,a2,a3,a4,a5]
def dict_to_sym(p):
    return sp.expand(sum(sp.Integer(c.numerator)/c.denominator*sp.prod(A[i]**e[i] for i in range(5)) for e,c in p.items()))
Q=pickle.load(open('Q8_quotients.pkl','rb'))
with open('factor.log','w') as lg:
    for X in ['A','B','CLASS','Z']:
        P=dict_to_sym(Q[X])
        # make integer primitive
        Psp=sp.Poly(P,*A)
        cont, Prim = Psp.primitive()
        lg.write(f"{X}: {len(Psp.terms())} terms, content {cont}\n"); lg.flush()
        t0=time.time()
        fl=sp.factor_list(Prim)
        lg.write(f"{X}: factor_list -> {len(fl[1])} factors [{time.time()-t0:.0f}s]\n")
        for fac,mult in fl[1]:
            lg.write(f"   mult {mult}, deg {sp.Poly(fac,*A).total_degree()}, terms {len(sp.Poly(fac,*A).terms())}\n")
        lg.flush()
        pickle.dump({'content':str(cont),'factors':[(str(fac),m) for fac,m in fl[1]]}, open(f'factor_{X}.pkl','wb'))
    lg.write("DONE\n")
print("factor done")
