import os, sys
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
import sympy as sp
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
SA=sp.sympify(open('subres1_A.txt').read())
SB=sp.sympify(open('subres1_B.txt').read())
PA=sp.Poly(SA,a1,a2,a4); PB=sp.Poly(SB,a1,a2,a4)
print("subres1_A: terms", len(PA.terms()), "total degree", sp.Poly(SA).total_degree(), "deg a4", sp.Poly(SA,a4).degree())
print("subres1_B: terms", len(PB.terms()), "total degree", sp.Poly(SB).total_degree(), "deg a4", sp.Poly(SB,a4).degree())
G=sp.sympify([l.split(':',1)[1].strip() for l in open('R_factors.txt') if 'deg 22' in l][0])
PG=sp.Poly(G,a1,a2,a4)
print("G22: terms", len(PG.terms()), "total degree", sp.Poly(G).total_degree(), "deg a4", sp.Poly(G,a4).degree())
# also C4, Q7' degrees
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
print("C4 terms", len(sp.Poly(C4).terms()), "deg", sp.Poly(C4).total_degree(), "deg a5", sp.Poly(C4,a5).degree())
print("Q7' terms", len(sp.Poly(Q7).terms()), "deg", sp.Poly(Q7).total_degree(), "deg a5", sp.Poly(Q7,a5).degree())
P7=sp.sympify(open('P7.txt').read())
print("P7 terms", len(sp.Poly(P7).terms()), "deg", sp.Poly(P7).total_degree(), "deg a5", sp.Poly(P7,a5).degree(), "deg a4", sp.Poly(P7,a4).degree())
