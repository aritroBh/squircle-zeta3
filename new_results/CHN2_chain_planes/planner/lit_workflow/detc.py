import sys, pickle; sys.dont_write_bytecode=True
sys.path.insert(0,'.')
import sympy as sp
a1,a2,a3,a4=sp.symbols('a1 a2 a3 a4')
Ms=pickle.load(open('M357_a.pkl','rb')); M3,M5,M7=Ms[3],Ms[5],Ms[7]
C=sp.expand(M5*M3-M3*M5)
print('c11 =',sp.factor(C[0,0]))
print('c12 =',sp.factor(C[0,1]))
print('c21 =',sp.factor(C[1,0]))
d=sp.factor(sp.expand(C.det()))
print('det[M5,M3] =',d)
# ratio r_k = M_k[0,1]/M_k[1,0]; commutator c11 = M5_12 M3_21 - M3_12 M5_21 ; Pi1,Pi2 <=> r_3 = r_5 ?
r3=sp.factor(M3[0,1]/M3[1,0]); print('r3 = M3_12/M3_21 =',r3)
print('M3_12 - M3_21 =',sp.factor(M3[0,1]-M3[1,0]))
print('M3_11 - M3_22 =',sp.factor(M3[0,0]-M3[1,1]))
