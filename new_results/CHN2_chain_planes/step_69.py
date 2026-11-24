
with open('progress.log','a') as f:
    f.write("2026-10-29 C1 engine: engine5.py symbolic P5 engine (A+A_rev single-series, Euler sums->MZV words, DS reduce) W<=10 -> forms10.pkl. Controls ALL PASS: u1=0/u5=0 == CHN1 P4 (5 pts, W<=10), u3=0/u2=0/u4=0 == Gamma-products, reversal symbolic, Cauchy/FFT c2..c8 at (1,2,4,3,5) reldiff ~1e-20. Two self-corrected algebra bugs (p_n subtract a4^(n-1) not a2^(n-1); q_m=(-1)^m(Y^m-X Y^(m-1))).\n")

import pickle
from fractions import Fraction as F
import sympy as sp
forms10 = pickle.load(open('forms10.pkl','rb'))
a1,a2,a3,a4,a5 = sp.symbols('a1 a2 a3 a4 a5')
AV = [a1,a2,a3,a4,a5]
def to_sym(p):
    return sp.sum([sp.Integer(c.numerator)/c.denominator * sp.prod(AV[i]**e[i] for i in range(5)) for e,c in p.items()], sp.Integer(0))

P8 = to_sym(forms10[8][3])   # z(3,5) at W=8
print("P8^(5) terms:", len(forms10[8][3]))
# divisibility checks
for i,v in enumerate(AV):
    print(f"  a{i+1} divides P8:", sp.expand(P8.subs(v,0))==0)
Sig = a1+a2+a3+a4+a5
print("Sigma divides P8:", sp.expand(sp.rem(P8, Sig, *AV))==0 if True else None)
# factor out a2a3a4
Q5 = sp.expand(P8/(a2*a3*a4))
print("Q5 = a2a3a4 quotient is polynomial:", sp.Poly(Q5, *AV).is_Poly if Q5!=0 else None)
print("Q5 terms:", len(sp.Poly(Q5,*AV).terms()))
# check Sigma | Q5
print("Sigma divides Q5:", sp.expand(sp.rem(Q5, Sig, *AV))==0)
