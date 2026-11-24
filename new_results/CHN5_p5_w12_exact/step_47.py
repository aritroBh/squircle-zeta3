import sympy as sp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
import os; os.chdir(W)
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
print("C4:", sp.Poly(C4).total_degree(), "deg;", "a5-deg", sp.Poly(C4,a5).degree(), "; terms", len(sp.Poly(C4).terms()))
print("Q7':", sp.Poly(Q7).total_degree(), "deg;", "a5-deg", sp.Poly(Q7,a5).degree(), "; terms", len(sp.Poly(Q7).terms()))
print("C4 head:", str(C4)[:300])
print("Q7' head:", str(Q7)[:300])
