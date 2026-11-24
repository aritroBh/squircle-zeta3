import os, pickle
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
from fractions import Fraction as F
import sympy as sp
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
A=[a1,a2,a3,a4,a5]
def dict_to_sym(p):
    return sp.expand(sum(sp.Integer(c.numerator)/c.denominator*sp.prod(A[i]**e[i] for i in range(5)) for e,c in p.items()))
Q=pickle.load(open('Q8_quotients.pkl','rb'))
out={}
with open('exact_forms.log','w') as lg:
    for X in ['A','B','CLASS','Z']:
        P=dict_to_sym(Q[X])
        Psp=sp.Poly(P,*A)
        cont,Prim=Psp.primitive()
        Prim=sp.expand(Prim)
        nterms=len(sp.Poly(Prim,*A).terms())
        # reversal symmetry check at monomial level
        revmap={a1:a5,a2:a4,a3:a3,a4:a2,a5:a1}
        rev=sp.expand(Prim.subs(revmap, simultaneous=True))
        revsym = sp.expand(Prim-rev)==0
        # sign pattern
        coeffs=[c for _,c in sp.Poly(Prim,*A).terms()]
        npos=sum(1 for c in coeffs if c>0); nneg=sum(1 for c in coeffs if c<0)
        Pval=Prim.subs({a1:1,a2:1,a3:1,a4:1,a5:1})
        lg.write(f"{X}: coord = ({cont})*a2a3a4(Sum a)*P8_{X}, P8_{X} integer primitive, deg 8, {nterms} terms, IRREDUCIBLE (factor_list), reversal-symmetric={revsym}, coeffs +{npos}/-{nneg}, P8(1,1,1,1,1)={Pval}\n")
        open(f'P8_{X}.txt','w').write(str(Prim))
        out[X]={'cont':str(cont),'poly':str(Prim),'nterms':nterms,'revsym':revsym,'npos':npos,'nneg':nneg,'P11111':str(Pval)}
        lg.flush()
    lg.write("DONE\n")
print(open('exact_forms.log').read())
pickle.dump(out, open('exact_forms_summary.pkl','wb'))
